"""Exp20: re-check of the published Exp17 statistics with the implicit load line (spatial.simulate implicit_load=True).
Re-runs all 40 Exp17 jobs (same parameters as experiments/exp17_silicon_split/run.py) and recomputes, exactly as
mechanism_check.py and figures17.py do: the centre-to-edge drop / dV_c at every split (exact configuration one step
before the split) by protocol, and the error of the first split width for the clean runs (R_L >= 4 kOhm, normal
mismatch, at most two filaments present).
usage: python recheck17.py run [nproc] | analyze   ->  results/exp20/recheck17/"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
import run as r17  # noqa: E402  (exp17 jobs, split_events, low priority)
from run import law, sp  # noqa: E402

OUT = ROOT / "results" / "exp20" / "recheck17"


def all_jobs():
    return r17.jobs() + r17.jobs_r5k() + r17.jobs_univ() + r17.jobs_nomax() + r17.jobs_slow()


def run(job):
    t0 = time.time()
    p = job.get("p", law.BASE)
    kw = dict(W=job["W"], dy=0.5, R_sq=1e12, R_L=job["R_L"], vdd_up=(0.6, job["vdd_max"]),
              dvdd=job.get("dvdd", 0.01 if job["R_L"] < 4e3 else 0.02), seed=job["seed"], r_D=job["r_D"],
              contacts="distributed", r_v=job["r_D"] * job["lam"] ** 2, down=False, sigma_scale=job["sigma_scale"],
              t_relax=job.get("t_relax", 1.0), implicit_load=True, **p)
    r = sp.simulate(**kw)
    c = law.constants(p)
    res = dict(job=job, prediction=dict(w_split=law.w_split(job["lam"], job["r_D"], c)), splits=r17.split_events(r["hist"]),
               n_max=max(h["n_fil"] for h in r["hist"]), hist=r["hist"], wall=time.time() - t0)
    (OUT / "runs").mkdir(parents=True, exist_ok=True)
    (OUT / "runs" / f"{job['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    return f"{job['name']}: {len(res['splits'])} bölünme ({res['wall']:.0f} s)"


def protocol(j):
    if j.get("p") and j["p"]["g_sub"] < 1e-14:
        return "low leakage"
    if j.get("p"):
        return "universality"
    return "stiff load (R_L 2k)" if j["R_L"] < 4e3 else "clean (R_L 5k)"


def analyze():
    c0 = law.constants()
    cache, events, first = {}, [], []
    for f in sorted((OUT / "runs").glob("*.json")):
        r = json.loads(f.read_text(encoding="utf-8"))
        j, h = r["job"], r["hist"]
        key = json.dumps(j.get("p"), sort_keys=True)
        if key not in cache:
            cache[key] = law.constants(j["p"]) if j.get("p") else c0
        c = cache[key]
        first_done = False
        for k in range(1, len(h)):
            P = [(cc - w / 2, cc + w / 2) for cc, w in zip(h[k - 1]["centres_um"], h[k - 1]["widths_um"])]
            Cn = [(cc - w / 2, cc + w / 2) for cc, w in zip(h[k]["centres_um"], h[k]["widths_um"])]
            for p in P:
                kids = [q for q in Cn if min(p[1], q[1]) - max(p[0], q[0]) > 0]
                if len(kids) >= 2:
                    du = law.du_config(P, p, j["W"], j["lam"], j["r_D"], c)
                    events.append(dict(run=f.stem, protocol=protocol(j), mm=j["sigma_scale"], ratio=du / c["dV_c"]))
                    if not first_done and c["regime"] == "Maxwell":
                        first_done = True
                        if j["R_L"] >= 4e3 and j["sigma_scale"] == 1.0 and len(P) <= 2:
                            first.append(dict(run=f.stem, err=(p[1] - p[0] - r["prediction"]["w_split"]) /
                                              r["prediction"]["w_split"]))
    groups = {}
    for prot in ("clean (R_L 5k)", "stiff load (R_L 2k)", "universality", "low leakage"):
        for mm in (1.0, 0.1):
            x = [e["ratio"] for e in events if e["protocol"] == prot and e["mm"] == mm]
            if x:
                groups[f"{prot}, mm={mm}"] = dict(n=len(x), mean=float(np.mean(x)), sd=float(np.std(x)),
                                                  min=float(np.min(x)), max=float(np.max(x)))
    err = [x["err"] for x in first]
    out = dict(groups=groups, first_split=dict(n=len(err), mean=float(np.mean(err)) if err else None,
                                               rms=float(np.sqrt(np.mean(np.square(err)))) if err else None, rows=first),
               published=dict(clean=dict(n=22, mean=1.001, sd=0.069), stiff=dict(n=22, mean=1.038, sd=0.062),
                              universality=dict(n=9, mean=0.979, sd=0.133), first_split=dict(n=14, mean=0.016, rms=0.091)))
    (OUT / "analysis.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


if __name__ == "__main__":
    if sys.argv[1] == "run":
        js = [j for j in all_jobs() if not (OUT / "runs" / f"{j['name']}.json").exists()]
        js.sort(key=lambda j: -j["vdd_max"])
        with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 3) as pool:
            for line in pool.imap_unordered(run, js):
                print(line, flush=True)
        print("TAMAM", flush=True)
    else:
        a = analyze()
        for k, v in a["groups"].items():
            print(f"{k:28s} n={v['n']:3d} oran {v['mean']:.3f} ± {v['sd']:.3f} ({v['min']:.2f}-{v['max']:.2f})")
        fs = a["first_split"]
        print(f"ilk bölünme (temiz): n={fs['n']} ortalama {fs['mean']:+.3f} RMS {fs['rms']:.3f} | yayımlanan n=14 +0.016 RMS 0.091")
