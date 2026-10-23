"""Exp20: controlling the sudden-current multiplication of exp19, and exploring pulses.
A ring  : ring (enclosed-gate) device, W = 600 um, 4 seeds at 75/225/375/525 um, all equivalent by symmetry, with unequal
          strengths 45-55 mV; slow sweep vs abrupt V_DD steps; do the k winners equal the k strongest seeds?
B cap   : common-node capacitance C_c (R_L C from 5 ns to 5 s): exp19 steps (V_s 1.1 / 1.4 V) and the exp18 abrupt kill.
C pulse : V_DD pulse V_p for tau, then V_hold = 1.2 V; does k need only the height, or a minimum duration?
D gate  : (exploration) a global gate pulse -Delta on every segment of an existing population: synchronous division?
E cycle : (exploration) repeated grow + gate-pulse cycles: 1 -> 2 -> 4 -> 8 cells on command?
usage: python run20.py predict <part> | run <nproc> <part> | analyze   ->  results/exp20/"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "exp18_programmable"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "exp19_pulse"))
import regen18 as rg  # noqa: E402  (imports run18: low priority, paths)
import run18  # noqa: E402
import run19  # noqa: E402
from run18 import law, sp, prog_array, DY, LAM, R_D, R_L  # noqa: E402

OUT = run18.ROOT / "results" / "exp20"
BASE = dict(dy=DY, R_L=R_L, r_D=R_D, contacts="distributed", r_v=R_D * LAM**2, down=False, sigma_scale=1.0, R_sq=1e12,
            **law.BASE)
NPHI = law.BASE["ss"] / np.log(10)
RING_FIELDS = tuple(range(500, 520))
RING_SEED_POS = (75.0, 225.0, 375.0, 525.0)
VS_A = (1.0, 1.05, 1.1, 1.2, 1.4)
CCS = (0.0, 1e-12, 1e-9, 1e-6, 1e-4, 1e-3)
TAUS = (0.02, 0.05, 0.1, 0.3, 1.0)


def ring_seeds(field):
    extra = np.random.default_rng(field + 2000).uniform(0.0, 0.010, 4)
    return [(c, 10.0, round(0.045 + float(e), 5)) for c, e in zip(RING_SEED_POS, extra)]


def mismatch(field, W):
    n = int(round(W / DY))
    return np.random.default_rng(field).normal(0.0, 3.356e-3 / np.sqrt(DY * 0.15), n)


def strengths(field, regions, W):
    """per seed: (largest cell overdrive in mV, log10 of the summed subthreshold weight sum exp(overdrive / n phi))"""
    n = int(round(W / DY))
    y = (np.arange(n) + 0.5) * DY
    ov = mismatch(field, W) + prog_array(regions, W)
    out = []
    for c, w, _ in regions:
        m = np.abs(y - c) < w / 2
        out.append((float(ov[m].max() * 1e3), float(np.log10(np.exp(ov[m] / NPHI).sum()))))
    return out


def evolve(V, vc, prev, vdd_of_t, regions_of_t, t_end, t_fine_end, seed, W, extra, log, dt_fine=0.02):
    t = 0.0
    while t < t_end - 1e-9:
        dt = dt_fine if t < t_fine_end - 1e-9 else 0.5
        t += dt
        v = vdd_of_t(t)
        r = sp.simulate(W=W, vdd_up=(v, v), dvdd=1.0, seed=seed, t_relax=dt, n_sub=2, prog=prog_array(regions_of_t(t), W),
                        V0=V, vc0=vc, **BASE, **extra)
        V, vc = np.array(r["V_final"]), r["vc_final"]
        now = rg.intervals(r["hist"][0])
        log.append(dict(t=round(t, 4), vdd=v, vc=vc, I_uA=r["hist"][0]["I_uA"], cells=now,
                        new=run19.new_cells(prev, now)))
        prev = now
    return V, vc


def on_regions(iv, regions):
    return [k for k, (c, w, _) in enumerate(regions) if min(iv[1], c + w / 2 + 2) - max(iv[0], c - w / 2 - 2) > 0]


# ---------------------------------------------------------------- A ring
def job_ring(j):
    regs = ring_seeds(j["field"])
    if j["vs"] is None:
        r = sp.simulate(W=600.0, vdd_up=(0.6, 1.6), dvdd=0.02, seed=j["field"], t_relax=1.0, prog=prog_array(regs, 600.0),
                        periodic=True, **BASE)
        b = next((h for h in r["hist"] if h["n_fil"] > 0), None)
        return dict(cells=rg.intervals(b) if b else [], regions=regs)
    log = []
    evolve(np.zeros(1200), None, [], lambda t: j["vs"], lambda t: regs, 4.0, 1.0, j["field"], 600.0,
           dict(periodic=True), log)
    return dict(cells=log[-1]["cells"], regions=regs, log=log)


# ---------------------------------------------------------------- B capacitance
def job_cap(j):
    log = []
    if j["kind"] == "step":
        evolve(np.zeros(1200), 0.6 if j["C_c"] > 0 else None, [], lambda t: j["vs"], lambda t: run19.SEEDS4, 12.0, 1.0,
               j["field"], 600.0, dict(C_c=j["C_c"]), log)
        return dict(cells=log[-1]["cells"], regions=run19.SEEDS4, log=log)
    g = json.loads((run18.ROOT / "results" / "exp18" / "regen" / "grow_v3.3.json").read_text(encoding="utf-8"))
    cells = rg.intervals(g["hist"][-1])
    c, w = run19.kill_zone(cells[-1])
    evolve(np.array(g["V_final"]), g["hist"][-1]["vds"] if j["C_c"] > 0 else None, cells, lambda t: 3.3,
           lambda t: [(c, w, -0.30)], 12.0, 1.0, 1, 600.0, dict(C_c=j["C_c"]), log)
    return dict(cells_before=cells, cells=log[-1]["cells"], log=log)


# ---------------------------------------------------------------- C pulse
def job_pulse(j):
    log = []
    evolve(np.zeros(1200), None, [], lambda t: j["vp"] if t <= j["tau"] + 1e-9 else 1.2, lambda t: run19.SEEDS4, 4.0,
           1.5, j["field"], 600.0, {}, log)
    return dict(cells=log[-1]["cells"], regions=run19.SEEDS4, log=log)


# ---------------------------------------------------------------- D dips on an existing population (exploration)
def job_dip(j):
    g = json.loads((run18.ROOT / "results" / "exp18" / "regen" / f"grow_v{j['vh']:g}.json").read_text(encoding="utf-8"))
    cells = rg.intervals(g["hist"][-1])
    vh, t0, t1 = j["vh"], 0.5, 0.5 + j["tau"]
    on = lambda t: t0 < t <= t1 + 1e-9  # noqa: E731
    if j["kind"] == "vdd":
        vdd_of_t, regs = (lambda t: vh - j["dv"] if on(t) else vh), (lambda t: [])
    else:
        vdd_of_t, regs = (lambda t: vh), (lambda t: [(300.0, 600.0, -j["dg"])] if on(t) else [])
    log = []
    evolve(np.array(g["V_final"]), None, cells, vdd_of_t, regs, 4.5, 2.0, 1, 600.0,
           dict(implicit_load=True) if j.get("implicit") else {}, log)
    return dict(cells_before=cells, cells=log[-1]["cells"], log=log)


def events(prev_cells, log):
    """splits / deaths / genuine births along a log, with the time of each event"""
    ev = []
    prev = prev_cells
    for x in log:
        now = x["cells"]
        for p in prev:
            kids = [q for q in now if min(p[1], q[1]) - max(p[0], q[0]) > 0]
            if len(kids) >= 2:
                ev.append(("split", x["t"]))
            elif not kids:
                ev.append(("death", x["t"]))
        ev += [("birth", x["t"]) for _ in x["new"]]
        prev = now
    return ev


# ---------------------------------------------------------------- C' carrying capacity after a 1.4 V pulse
def job_carry(j):
    log = []
    evolve(np.zeros(1200), None, [], lambda t: 1.4 if t <= 0.1 + 1e-9 else j["vhold"], lambda t: run19.SEEDS4, 6.0, 1.5,
           j["field"], 600.0, {}, log)
    return dict(cells=log[-1]["cells"], regions=run19.SEEDS4, log=log)


PRED = {
    "ring2": dict(
        basis="exp20 A on fields 500-519: at abrupt steps the winners matched the k seeds with the largest summed "
              "subthreshold weight (sum exp(overdrive / n phi) over the seed) in 14-15 of 19-20 fields, the single "
              "strongest cell in 5-12; this predictor was chosen after the fact and is tested here on new fields",
        A4="ring, new fields 520-539: at V_s = 1.05, 1.1 and 1.2 V the winners are the k seeds with the largest summed "
           "subthreshold weight in >= 70 % of the fields at each V_s; k = 2, 2, 3 as on fields 500-519"),
    "carry": dict(
        C3="after 0.1 s at 1.4 V (four seeds become cells), the number that survive grows with the holding voltage "
           "V_hold and is at most 4; at V_hold = 1.4 V all four survive (step result); at low V_hold cells die because the "
           "current budget (V_hold - V_c) / R_L cannot feed four cells"),
    "dip_imp": dict(
        basis="the explicit dip runs (exp20 D) ended in a non-equilibrium state after sudden current excesses (rail not "
              "self-consistent); they are repeated with implicit_load=True. D1-D3 below were registered before both batches. "
              "D4 comes from one diagnostic run (5.0 V, -60 mV gate dip for 0.5 s: 3 -> 7 cells at the release) and is "
              "tested on the whole batch",
        D1="V_DD dip: a cell divides when the dip depth exceeds its margin (3.3 V: 7.8/10.2/8.2 mV; 5.0 V: 3.2/2.7/2.8 mV), "
           "if the dip lasts long enough",
        D2="during a global gate dip there is no division (the margin V_M - V_down grows)",
        D3="0.5 V dips kill cells instead of dividing them",
        D4="releasing a global gate dip divides cells: the dip widens them (lower current density), the release forces a "
           "sudden shrink; divisions at the release grow with the dip depth and are more frequent at 5.0 V than at 3.3 V; "
           "at 5.0 V with -60 mV, 3 -> at least 6 cells for both durations"),
    "dip": dict(
        margins_mV={"3.3": [7.8, 10.2, 8.2], "5.0": [3.2, 2.7, 2.8]},
        D1="V_DD dip of depth dV: a cell divides when dV exceeds its quasi-static margin dV_c - du (exact configuration, "
           "values above), if the dip lasts long enough for the centre to quench before the fronts retreat; so at 5.0 V "
           "dips of 5 mV and deeper divide all three cells, at 3.3 V dips of 20 mV and deeper",
        D2="global gate dip (-20/-40/-60 mV on every segment): no division at any depth, because the margin V_M - V_down "
           "grows to 19.2/22.1/25.0 mV; the cell current density falls by 9/18/27 %, the common node rises and the cells "
           "widen; at -60 mV new cells may nucleate",
        D3="very deep V_DD dips (0.5 V) kill cells instead of dividing them"),
    "ring": dict(
        A1="ring, slow sweep: the first cell is born at the seed with the largest single-cell overdrive in >= 16 of 20 fields",
        A2="ring, abrupt steps: whenever 0 < k < 4 seeds become cells, they are the k seeds with the largest single-cell "
           "overdrive in >= 70 % of the cases (exp19 without a ring: 20-30 %)",
        A3="k rises with V_s and all four seeds become cells at V_s = 1.4 V",
        secondary="the summed subthreshold weight per seed is evaluated as a second predictor, without a prediction"),
    "cap": dict(
        B1="steps (V_s 1.1 and 1.4 V): R_L C <= 5 us (on-chip node, 1 pF -> 5 ns) gives exactly the k of C_c = 0 in every "
           "field; R_L C >= 0.5 s makes k fall towards 1 (the node rises like a slow sweep)",
        B2="abrupt kill at 3.3 V: births unchanged for R_L C <= 5 us, none for R_L C >= 0.5 s"),
    "pulse": dict(
        C1="after a pulse of height V_p and duration tau (then 1.2 V), the number of seed cells grows with tau and saturates "
           "at the exp19 step value (V_p 1.1 V: 2, 1.4 V: 4)",
        C2="the shortest tau that reaches the step value is shorter for V_p = 1.4 V than for 1.1 V"),
}


def jobs(part):
    if part == "ring":
        return [dict(part=part, name=f"ring_f{f}_{'slow' if vs is None else f'v{vs:g}'}", field=f, vs=vs)
                for f in RING_FIELDS for vs in (None,) + VS_A]
    if part == "cap":
        out = [dict(part=part, name=f"cap_step_f{f}_v{vs:g}_C{cc:.0e}", kind="step", field=f, vs=vs, C_c=cc)
               for f in (400, 401, 402, 403, 404) for vs in (1.1, 1.4) for cc in CCS]
        return out + [dict(part=part, name=f"cap_kill_C{cc:.0e}", kind="kill", C_c=cc) for cc in CCS]
    if part == "pulse":
        return [dict(part=part, name=f"pulse_f{f}_vp{vp:g}_tau{tau:g}", field=f, vp=vp, tau=tau)
                for f in (400, 401, 402, 403, 404) for vp in (1.1, 1.4) for tau in TAUS]
    if part == "ring2":
        return [dict(part="ring", name=f"ring_f{f}_v{vs:g}", field=f, vs=vs) for f in range(520, 540) for vs in (1.05, 1.1, 1.2)]
    if part == "carry":
        return [dict(part=part, name=f"carry_f{f}_vh{vh:g}", field=f, vhold=vh)
                for f in (400, 401, 402, 403, 404) for vh in (0.95, 1.0, 1.1, 1.2, 1.3, 1.4)]
    if part == "dip_imp":
        return [dict(j, part="dip", implicit=True, name=j["name"].replace("dip_", "dipimp_")) for j in jobs("dip")]
    if part == "dip":
        out = [dict(part=part, name=f"dip_vdd_v{vh:g}_dv{dv:g}_tau{tau:g}", kind="vdd", vh=vh, dv=dv, tau=tau)
               for vh in (3.3, 5.0) for dv in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5) for tau in (0.02, 0.1, 0.5)]
        return out + [dict(part=part, name=f"dip_gate_v{vh:g}_dg{dg:g}_tau{tau:g}", kind="gate", vh=vh, dg=dg, tau=tau)
                      for vh in (3.3, 5.0) for dg in (0.02, 0.04, 0.06) for tau in (0.1, 0.5)]
    raise ValueError(part)


def run(j):
    t0 = time.time()
    res = {"ring": job_ring, "cap": job_cap, "pulse": job_pulse, "dip": job_dip, "carry": job_carry}[j["part"]](j)
    res.update(job=j, wall=time.time() - t0)
    (OUT / "runs").mkdir(parents=True, exist_ok=True)
    (OUT / "runs" / f"{j['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    return f"{j['name']}: {len(res['cells'])} hücre ({res['wall']:.0f} s)"


def analyze():
    runs = {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in (OUT / "runs").glob("*.json")}
    a = {"ring": {}, "cap": {}, "pulse": {}}
    for vs in (None,) + VS_A:
        key = "slow" if vs is None else f"{vs:g}"
        rows = []
        for f in RING_FIELDS:
            r = runs.get(f"ring_f{f}_{key if vs is None else 'v' + key}")
            if not r:
                continue
            regs = [tuple(x) for x in r["regions"]]
            win = sorted({k for iv in r["cells"] for k in on_regions(iv, regs)})
            off = [iv for iv in r["cells"] if not on_regions(iv, regs)]
            st = strengths(f, regs, 600.0)
            k = len(win)
            top_max = sorted(sorted(range(4), key=lambda i: -st[i][0])[:k])
            top_sum = sorted(sorted(range(4), key=lambda i: -st[i][1])[:k])
            rows.append(dict(field=f, k=k, off=len(off), winners=win, top_max=top_max, top_sum=top_sum,
                             hit_max=win == top_max, hit_sum=win == top_sum, st=st))
        part = [x for x in rows if 0 < x["k"] < 4]
        a["ring"][key] = dict(n=len(rows), mean_k=float(np.mean([x["k"] for x in rows])) if rows else None,
                              k_counts={str(k): sum(x["k"] == k for x in rows) for k in range(5)},
                              off=sum(x["off"] for x in rows), n_partial=len(part),
                              hit_max=sum(x["hit_max"] for x in part), hit_sum=sum(x["hit_sum"] for x in part), rows=rows)
    for vs in (1.1, 1.4):
        for cc in CCS:
            ks = []
            for f in (400, 401, 402, 403, 404):
                r = runs.get(f"cap_step_f{f}_v{vs:g}_C{cc:.0e}")
                if r:
                    ks.append(len({k for iv in r["cells"] for k in run19.on_seed(iv)}))
            if ks:
                a["cap"][f"step|{vs:g}|{cc:.0e}"] = dict(k=ks, mean=float(np.mean(ks)))
    for cc in CCS:
        r = runs.get(f"cap_kill_C{cc:.0e}")
        if r:
            a["cap"][f"kill|{cc:.0e}"] = dict(births=sum(len(x["new"]) for x in r["log"]),
                                              persistent=len(run19.persistent_new(r["cells_before"], r["cells"])),
                                              vc_max=max(x["vc"] for x in r["log"]))
    for vp in (1.1, 1.4):
        for tau in TAUS:
            ks = []
            for f in (400, 401, 402, 403, 404):
                r = runs.get(f"pulse_f{f}_vp{vp:g}_tau{tau:g}")
                if r:
                    ks.append(len({k for iv in r["cells"] for k in run19.on_seed(iv)}))
            if ks:
                a["pulse"][f"{vp:g}|{tau:g}"] = dict(k=ks, mean=float(np.mean(ks)))
    for vs in (1.05, 1.1, 1.2):                       # ring2: new fields 520-539 (pre-registered summed predictor)
        rows = []
        for f in range(520, 540):
            r = runs.get(f"ring_f{f}_v{vs:g}")
            if not r:
                continue
            regs = [tuple(x) for x in r["regions"]]
            win = sorted({k for iv in r["cells"] for k in on_regions(iv, regs)})
            st = strengths(f, regs, 600.0)
            k = len(win)
            top_sum = sorted(sorted(range(4), key=lambda i: -st[i][1])[:k])
            top_max = sorted(sorted(range(4), key=lambda i: -st[i][0])[:k])
            rows.append(dict(field=f, k=k, hit_sum=win == top_sum, hit_max=win == top_max,
                             off=sum(1 for iv in r["cells"] if not on_regions(iv, regs))))
        if rows:
            part = [x for x in rows if 0 < x["k"] < 4]
            a.setdefault("ring2", {})[f"{vs:g}"] = dict(n=len(rows), k_counts={str(k): sum(x["k"] == k for x in rows)
                                                                              for k in range(5)},
                                                       off=sum(x["off"] for x in rows), n_partial=len(part),
                                                       hit_sum=sum(x["hit_sum"] for x in part),
                                                       hit_max=sum(x["hit_max"] for x in part))
    for vh in (0.95, 1.0, 1.1, 1.2, 1.3, 1.4):
        ks = []
        for f in (400, 401, 402, 403, 404):
            r = runs.get(f"carry_f{f}_vh{vh:g}")
            if r:
                ks.append(len({k for iv in r["cells"] for k in run19.on_seed(iv)}))
        if ks:
            a.setdefault("carry", {})[f"{vh:g}"] = dict(k=ks, mean=float(np.mean(ks)))
    a["dip"] = {}
    for name, r in sorted(runs.items()):
        if name.startswith("dip_"):
            ev = events(r["cells_before"], r["log"])
            a["dip"][name[4:]] = dict(n_before=len(r["cells_before"]), n_end=len(r["cells"]),
                                      splits=sum(e[0] == "split" for e in ev), deaths=sum(e[0] == "death" for e in ev),
                                      births=sum(e[0] == "birth" for e in ev),
                                      vc_max=max(x["vc"] for x in r["log"]), vc_min=min(x["vc"] for x in r["log"]))
    (OUT / "analysis.json").write_text(json.dumps(a, indent=1), encoding="utf-8")
    return a


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[1] == "predict":
        p = OUT / f"predictions_{sys.argv[2]}.json"
        if p.exists():
            print("zaten var, üzerine yazılmadı:", p)
        else:
            p.write_text(json.dumps(PRED[sys.argv[2]], indent=1, ensure_ascii=False), encoding="utf-8")
            print("yazıldı:", p)
    elif sys.argv[1] == "run":
        js = [j for j in jobs(sys.argv[3]) if not (OUT / "runs" / f"{j['name']}.json").exists()]
        with Pool(int(sys.argv[2])) as pool:
            for line in pool.imap_unordered(run, js):
                print(line, flush=True)
        print("TAMAM", flush=True)
    else:
        a = analyze()
        for k, v in a["ring"].items():
            print(f"ring {k:5s}: k dağılımı {v['k_counts']} tohum dışı {v['off']} | en güçlü-k (tek hücre) {v['hit_max']}/{v['n_partial']}"
                  f" | en güçlü-k (toplam ağırlık) {v['hit_sum']}/{v['n_partial']}")
        for k, v in a["cap"].items():
            print("cap", k, v)
        for k, v in a["pulse"].items():
            print("pulse", k, v)
