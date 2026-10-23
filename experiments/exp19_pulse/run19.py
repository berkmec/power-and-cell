"""Exp19: cells that multiply after a sudden current change (follow-up of exp18 v2, where an abrupt kill made the common
node jump and new cells nucleate at the nearest free healthy places).
A killspeed: the exp18 target (seed 1, V_hold 3.3 / 5.0 V) is killed with a -300 mV band (30 um margins) that ramps up in
  T_r = 0 ... 30 s; 20 ms steps until T_r + 1 s, then 0.5 s steps for 10 s more.
B rsq: abrupt kill (T_r = 0) at 3.3 V with lateral body sheet resistance R_sq = 1e10, 1e11, 1e12 (fresh growth sweeps).
C step: 4 equal seeds (+50 mV, 10 um) at 75, 225, 375, 525 um; slow sweep (0.02 V/s) vs abrupt V_DD steps from rest to V_s
  (20 ms steps for 1 s, then 0.5 s steps up to 4 s); 10 random fields (seeds 400-409).
usage: python run19.py predict | run [nproc] [killspeed|rsq|step] | analyze   ->  results/exp19/"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "exp18_programmable"))
import regen18 as rg  # noqa: E402  (imports run18: low priority, paths)
import run18  # noqa: E402
from run18 import law, sp, prog_array, W, DY, LAM, R_D, R_L  # noqa: E402

OUT = run18.ROOT / "results" / "exp19"
BASEKW = dict(W=W, dy=DY, R_L=R_L, r_D=R_D, contacts="distributed", r_v=R_D * LAM**2, down=False, sigma_scale=1.0,
              **law.BASE)
TRS = (0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0)
RSQS = (1e10, 1e11, 1e12)
SEEDS4 = [(75.0, 10.0, 0.05), (225.0, 10.0, 0.05), (375.0, 10.0, 0.05), (525.0, 10.0, 0.05)]
VSTEPS = (0.95, 1.0, 1.05, 1.1, 1.2, 1.4)
FIELDS = tuple(range(400, 410))

PREDICTIONS = dict(
    written_before_runs=True,
    A1="births (new cells farther than 10 um from all earlier cells) do not increase with the ramp time T_r",
    A2="T_r = 0 reproduces exp18 v2 (3.3 V: 2 births, 5.0 V: 3 births); T_r = 30 s reproduces exp18 v3 (0 births)",
    A3="the shortest T_r with zero births lies between 0.1 s and 10 s at both drives",
    B1="hypothesis, uncertain: stronger lateral body coupling (smaller R_sq) lets existing fronts absorb the released current "
       "faster, so an abrupt kill gives no more births at 1e10 / 1e11 than at 1e12",
    C1="slow sweep: exactly one of the four +50 mV seeds hosts the first birth in >= 7 of 10 fields (winner takes all)",
    C2="abrupt steps: the number of seeds that become cells does not decrease with V_s; at V_s = 1.05 and 1.1 V all four seeds "
       "become cells in >= 7 of 10 fields with no cell outside the seeds",
    C3="at V_s = 1.4 V cells also appear outside the seeds (background cells with V_up < 1.4 V)",
    design="V_up(local overdrive) from law.constants: +30 mV 1.18 V, +50 mV 1.02 V, +70 mV 0.89 V; seed regions reach about "
           "+50 + 20 mV, background maxima about +20...+40 mV (exp18 births without seeds at V_DD 1.18-1.34 V)",
)


def new_cells(prev, now):
    return [iv for iv in now if all(iv[1] < a - 10 or iv[0] > b + 10 for a, b in prev)]


def evolve(V, prev, vdd, regions_of_t, t_end, t_fine_end, seed, R_sq, log, dt_fine=0.02):
    t = 0.0
    while t < t_end - 1e-9:
        dt = dt_fine if t < t_fine_end - 1e-9 else 0.5
        t += dt
        r = sp.simulate(vdd_up=(vdd, vdd), dvdd=1.0, seed=seed, R_sq=R_sq, t_relax=dt, n_sub=2,
                        prog=prog_array(regions_of_t(t)), V0=V, **BASEKW)
        V = np.array(r["V_final"])
        h = r["hist"][0]
        now = rg.intervals(h)
        log.append(dict(t=round(t, 4), vc=h["vds"], I_uA=h["I_uA"], cells=now, new=new_cells(prev, now)))
        prev = now
    return V


def kill_zone(target):
    return ((target[0] + target[1]) / 2, target[1] - target[0] + 60.0)


def job_killspeed(j):
    g = json.loads((run18.ROOT / "results" / "exp18" / "regen" / f"grow_v{j['vhold']:g}.json").read_text(encoding="utf-8"))
    cells = rg.intervals(g["hist"][-1])
    c, w = kill_zone(cells[-1])
    Tr = j["Tr"]
    log = []
    evolve(np.array(g["V_final"]), cells, j["vhold"], lambda t: [(c, w, -0.30 * (min(1.0, t / Tr) if Tr else 1.0))],
           Tr + 11.0, Tr + 1.0, 1, 1e12, log)
    return dict(cells_before=cells, log=log)


def job_rsq(j):
    g = sp.simulate(vdd_up=(0.6, 3.3), dvdd=0.02, seed=1, R_sq=j["R_sq"], t_relax=1.0, **BASEKW)
    cells = rg.intervals(g["hist"][-1])
    log = []
    if cells:
        c, w = kill_zone(cells[-1])
        evolve(np.array(g["V_final"]), cells, 3.3, lambda t: [(c, w, -0.30)], 11.0, 1.0, 1, j["R_sq"], log)
    return dict(cells_before=cells, log=log)


def job_step(j):
    if j["vs"] is None:                                   # slow sweep, the exp18 protocol
        r = sp.simulate(vdd_up=(0.6, 1.6), dvdd=0.02, seed=j["field"], R_sq=1e12, t_relax=1.0, prog=prog_array(SEEDS4),
                        **BASEKW)
        b = next((h for h in r["hist"] if h["n_fil"] > 0), None)
        return dict(first_birth=rg.intervals(b) if b else [], final=rg.intervals(r["hist"][-1]),
                    first_birth_vdd=b["vdd"] if b else None)
    log = []
    evolve(np.zeros(int(round(W / DY))), [], j["vs"], lambda t: SEEDS4, 4.0, 1.0, j["field"], 1e12, log)
    return dict(final=log[-1]["cells"], log=log)


def run(j):
    t0 = time.time()
    res = {"killspeed": job_killspeed, "rsq": job_rsq, "step": job_step}[j["mode"]](j)
    res.update(job=j, wall=time.time() - t0)
    (OUT / "runs").mkdir(parents=True, exist_ok=True)
    (OUT / "runs" / f"{j['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    n_new = sum(len(x["new"]) for x in res.get("log", []))
    return f"{j['name']}: yeni hücre {n_new} | son {len(res.get('final') or (res['log'][-1]['cells'] if res.get('log') else []))} ({res['wall']:.0f} s)"


def jobs(mode):
    if mode == "killspeed":
        return [dict(mode=mode, name=f"kill_v{vh:g}_Tr{Tr:g}", vhold=vh, Tr=Tr) for vh in (3.3, 5.0) for Tr in TRS]
    if mode == "rsq":
        return [dict(mode=mode, name=f"rsq_{R:.0e}", R_sq=R) for R in RSQS]
    return [dict(mode=mode, name=f"step_f{f}_{'slow' if vs is None else f'v{vs:g}'}", field=f, vs=vs)
            for f in FIELDS for vs in (None,) + VSTEPS]


def on_seed(iv):
    return [k for k, (c, w, _) in enumerate(SEEDS4) if min(iv[1], c + w / 2 + 2) - max(iv[0], c - w / 2 - 2) > 0]


def persistent_new(before, after):
    """cells present at the end that overlap none of the cells before the kill (transient fragments are not counted)"""
    return [iv for iv in after if all(min(iv[1], b) - max(iv[0], a) <= 0 for a, b in before)]


def seed_strengths(field):
    """largest local overdrive (mismatch + programmed offset, mV) inside each seed region: what decides the race"""
    n = int(round(W / DY))
    y = (np.arange(n) + 0.5) * DY
    d = np.random.default_rng(field).normal(0.0, 3.356e-3 / np.sqrt(DY * 0.15), n)
    return [float((d[np.abs(y - c) < w / 2].max() + dv) * 1e3) for c, w, dv in SEEDS4]


def analyze():
    runs = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in (OUT / "runs").glob("*.json")}
    out = {"killspeed": {}, "rsq": {}, "step": {}}
    for vh in (3.3, 5.0):
        for Tr in TRS:
            r = runs.get(f"kill_v{vh:g}_Tr{Tr:g}")
            if r:
                out["killspeed"][f"{vh:g}|{Tr:g}"] = dict(births=sum(len(x["new"]) for x in r["log"]),
                                                         persistent_new=len(persistent_new(r["cells_before"], r["log"][-1]["cells"])),
                                                         vc_max=max(x["vc"] for x in r["log"]),
                                                         n_before=len(r["cells_before"]), n_end=len(r["log"][-1]["cells"]))
    for R in RSQS:
        r = runs.get(f"rsq_{R:.0e}")
        if r:
            out["rsq"][f"{R:.0e}"] = dict(cells_before=r["cells_before"], births=sum(len(x["new"]) for x in r["log"]),
                                          persistent_new=len(persistent_new(r["cells_before"], r["log"][-1]["cells"])) if r["log"] else None,
                                          cells_end=r["log"][-1]["cells"] if r["log"] else None,
                                          n_end=len(r["log"][-1]["cells"]) if r["log"] else None)
    for vs in (None,) + VSTEPS:
        key = "slow" if vs is None else f"{vs:g}"
        rows = []
        for f in FIELDS:
            r = runs.get(f"step_f{f}_{key if vs is None else 'v' + key}")
            if not r:
                continue
            cells = r["first_birth"] if vs is None else r["final"]
            seeds_hit = sorted({k for iv in cells for k in on_seed(iv)})
            off = [iv for iv in cells if not on_seed(iv)]
            st = seed_strengths(f)
            k = len(seeds_hit)
            order = sorted(range(4), key=lambda i: -st[i])
            gap = (st[order[k - 1]] - st[order[k]]) if 0 < k < 4 else None
            rows.append(dict(field=f, n_cells=len(cells), seeds_with_cell=k, cells_off_seed=len(off), winners=seeds_hit,
                             strengths_mV=[round(x, 1) for x in st], topk=sorted(order[:k]),
                             winners_are_topk=sorted(order[:k]) == seeds_hit, gap_kth_mV=gap))
        if rows:
            out["step"][key] = dict(rows=rows, all4_no_off=sum(r_["seeds_with_cell"] == 4 and r_["cells_off_seed"] == 0
                                                              for r_ in rows),
                                    exactly1=sum(r_["seeds_with_cell"] == 1 for r_ in rows),
                                    winners_are_topk=sum(r_["winners_are_topk"] for r_ in rows if 0 < r_["seeds_with_cell"] < 4),
                                    n_partial=sum(0 < r_["seeds_with_cell"] < 4 for r_ in rows),
                                    mean_seeds=float(np.mean([r_["seeds_with_cell"] for r_ in rows])),
                                    mean_off=float(np.mean([r_["cells_off_seed"] for r_ in rows])))
    (OUT / "analysis.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    mode = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    if mode == "predict":
        p = OUT / "predictions.json"
        if p.exists():
            print("zaten var, üzerine yazılmadı:", p)
        else:
            p.write_text(json.dumps(PREDICTIONS, indent=1, ensure_ascii=False), encoding="utf-8")
            print("yazıldı:", p)
    elif mode == "run":
        sub = sys.argv[3]
        js = [j for j in jobs(sub) if not (OUT / "runs" / f"{j['name']}.json").exists()]
        js.sort(key=lambda j: -(j.get("Tr") or 0))                  # longest first
        with Pool(int(sys.argv[2])) as pool:
            for line in pool.imap_unordered(run, js):
                print(line, flush=True)
        print("TAMAM", flush=True)
    else:
        print(json.dumps(analyze(), indent=1)[:6000])
