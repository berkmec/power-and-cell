"""Exp20: does the implicit load line (implicit_load=True) change earlier results?
Re-runs a subset of exp17 (quasi-static sweeps), exp19 (abrupt steps, abrupt kill) and exp20 A (ring) with the implicit
solver and compares with the stored explicit-solver results.  ->  results/exp20/validate_implicit.json"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run20 as r  # noqa: E402
from run20 import sp, law, rg, run19, run18  # noqa: E402

sys.path.insert(0, str(run18.ROOT / "experiments" / "exp17_silicon_split"))
from run import split_events  # noqa: E402  (exp17)

OUT = r.OUT / "validate"
E17 = ("H_lam60_s1", "E_lam60", "A_lam60", "G_rD20_lam60", "U1_x1.5")


def job(j):
    t0 = time.time()
    if j["kind"] == "e17":
        old = json.loads((run18.ROOT / "results" / "exp17" / "runs" / f"{j['run']}.json").read_text(encoding="utf-8"))
        jb = old["job"]
        p = jb.get("p", law.BASE)
        kw = dict(W=jb["W"], dy=0.5, R_sq=1e12, R_L=jb["R_L"], vdd_up=(0.6, jb["vdd_max"]),
                  dvdd=jb.get("dvdd", 0.01 if jb["R_L"] < 4e3 else 0.02), seed=jb["seed"], r_D=jb["r_D"],
                  contacts="distributed", r_v=jb["r_D"] * jb["lam"] ** 2, down=False, sigma_scale=jb["sigma_scale"],
                  t_relax=jb.get("t_relax", 1.0), **p)
        new = sp.simulate(implicit_load=True, **kw)
        so, sn = old["splits"], split_events(new["hist"])
        res = dict(old=[(round(e["vdd"], 2), round(e["parent_w"], 1)) for e in so],
                   new=[(round(e["vdd"], 2), round(e["parent_w"], 1)) for e in sn],
                   n_end_old=old["hist"][-1]["n_fil"], n_end_new=new["hist"][-1]["n_fil"])
    elif j["kind"] == "step":
        log = []
        r.evolve(np.zeros(1200), None, [], lambda t: j["vs"], lambda t: run19.SEEDS4, 4.0, 1.0, j["field"], 600.0,
                 dict(implicit_load=True), log)
        old = json.loads((run18.ROOT / "results" / "exp19" / "runs" / f"step_f{j['field']}_v{j['vs']:g}.json")
                         .read_text(encoding="utf-8"))
        res = dict(old=sorted({k for iv in old["final"] for k in run19.on_seed(iv)}),
                   new=sorted({k for iv in log[-1]["cells"] for k in run19.on_seed(iv)}),
                   off_new=sum(1 for iv in log[-1]["cells"] if not run19.on_seed(iv)))
    elif j["kind"] == "ring":
        regs = r.ring_seeds(j["field"])
        log = []
        r.evolve(np.zeros(1200), None, [], lambda t: j["vs"], lambda t: regs, 4.0, 1.0, j["field"], 600.0,
                 dict(periodic=True, implicit_load=True), log)
        old = json.loads((r.OUT / "runs" / f"ring_f{j['field']}_v{j['vs']:g}.json").read_text(encoding="utf-8"))
        res = dict(old=sorted({k for iv in old["cells"] for k in r.on_regions(iv, regs)}),
                   new=sorted({k for iv in log[-1]["cells"] for k in r.on_regions(iv, regs)}))
    else:  # abrupt kill, exp19 protocol
        g = json.loads((run18.ROOT / "results" / "exp18" / "regen" / f"grow_v{j['vh']:g}.json").read_text(encoding="utf-8"))
        cells = rg.intervals(g["hist"][-1])
        c, w = run19.kill_zone(cells[-1])
        log = []
        r.evolve(np.array(g["V_final"]), None, cells, lambda t: j["vh"], lambda t: [(c, w, -0.30)], 11.0, 1.0, 1, 600.0,
                 dict(implicit_load=True), log)
        old = json.loads((run18.ROOT / "results" / "exp19" / "runs" / f"kill_v{j['vh']:g}_Tr0.json").read_text(encoding="utf-8"))
        res = dict(old_births=sum(len(x["new"]) for x in old["log"]), new_births=sum(len(x["new"]) for x in log),
                   old_end=len(old["log"][-1]["cells"]), new_end=len(log[-1]["cells"]))
    res.update(job=j, wall=time.time() - t0)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{j['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    return f"{j['name']}: {({k: v for k, v in res.items() if k not in ('job', 'wall')})} ({res['wall']:.0f} s)"


if __name__ == "__main__":
    js = [dict(kind="e17", run=x, name=f"e17_{x}") for x in E17]
    js += [dict(kind="step", field=f, vs=vs, name=f"step_f{f}_v{vs:g}") for f in (400, 401, 402, 403, 404)
           for vs in (1.05, 1.1, 1.2, 1.4)]
    js += [dict(kind="ring", field=f, vs=vs, name=f"ring_f{f}_v{vs:g}") for f in range(500, 510) for vs in (1.05, 1.2)]
    js += [dict(kind="kill", vh=vh, name=f"kill_v{vh:g}") for vh in (3.3, 5.0)]
    js = [j for j in js if not (OUT / f"{j['name']}.json").exists()]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as pool:
        for line in pool.imap_unordered(job, js):
            print(line, flush=True)
    print("TAMAM", flush=True)
