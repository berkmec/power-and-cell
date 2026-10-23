"""Exp18 predictions, written BEFORE the runs they predict (write-once files in results/exp18/).
calib: which smoothing length ell of the random field (mismatch d + prog) predicts the first birth site best (seeds 10-29).
place/wta/two: per-run predicted first birth site = argmax of the smoothed field with that ell, and Monte Carlo
P(first birth in the seed region | offset).
usage: python predict18.py calib | place | wta | two"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run18  # noqa: E402

RES = run18.ROOT / "results" / "exp18"
SIG = 3.356e-3 / np.sqrt(run18.DY * 0.15)            # per-cell Pelgrom sigma used by spatial.simulate (sigma_scale 1)
N = int(round(run18.W / run18.DY))
Y = (np.arange(N) + 0.5) * run18.DY
ELLS = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0)


def field(seed, regions):
    """Exactly the per-cell overdrive offsets of spatial.simulate: same RNG stream, plus the programmed regions."""
    return np.random.default_rng(seed).normal(0.0, SIG, N) + run18.prog_array(regions)


def smooth(f, ell):
    return gaussian_filter1d(f, ell / run18.DY, mode="reflect") if ell > 0 else f


def born(run):
    b = run["birth"]
    return [(c - w / 2, c + w / 2) for c, w in zip(b["centres_um"], b["widths_um"])] if b else []


def separated_maxima(f, k, sep=20.0):
    out = []
    for i in np.argsort(f)[::-1]:
        if all(abs(Y[i] - y) > sep for y in out):
            out.append(float(Y[i]))
        if len(out) == k:
            break
    return out


def calibrate():
    runs = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((RES / "runs").glob("calib_s*.json"))]
    table = {}
    for ell in ELLS:
        top1 = sum(any(s - 1 <= separated_maxima(smooth(field(r["job"]["seed"], []), ell), 1)[0] <= e + 1
                       for s, e in born(r)) for r in runs)
        # every born filament contains one of the k highest separated maxima (k = number of births)
        allb = sum(all(any(s - 1 <= y <= e + 1 for y in separated_maxima(smooth(field(r["job"]["seed"], []), ell),
                                                                           len(born(r))))
                       for s, e in born(r)) for r in runs)
        table[str(ell)] = dict(top1_hit=top1, all_births_explained=allb, n_runs=len(runs))
    best = max(ELLS, key=lambda e: (table[str(e)]["top1_hit"], table[str(e)]["all_births_explained"], -e))
    n_births = [len(born(r)) for r in runs]
    return dict(table=table, ell_best=best, births_per_run=n_births,
                birth_vdd=[r["birth"]["vdd"] for r in runs])


def in_region(y, c, w, tol=2.0):
    return abs(y - c) <= w / 2 + tol


def predict(mode, ell):
    jobs = {"place": run18.jobs_place, "wta": run18.jobs_wta, "two": run18.jobs_two}[mode]()
    per_run = {}
    for j in jobs:
        f = smooth(field(j["seed"], j["regions"]), ell)
        y_pred = separated_maxima(f, 1)[0]
        hit = [k for k, (c, w, dv) in enumerate(j["regions"]) if in_region(y_pred, c, w)]
        per_run[j["name"]] = dict(y_pred=y_pred, region=hit[0] if hit else None)
    out = dict(ell=ell, per_run=per_run)
    if mode == "place":
        rng = np.random.default_rng(12345)
        curve = {}
        for dv in sorted({j["regions"][0][2] for j in jobs}):
            reg = run18.prog_array([(150.0, 10.0, dv)])
            ok = 0
            for _ in range(4000):
                f = smooth(rng.normal(0.0, SIG, N) + reg, ell)
                ok += in_region(Y[np.argmax(f)], 150.0, 10.0)
            curve[f"{dv:.3f}"] = ok / 4000
        out["P_birth_in_seed_MC"] = curve
    if mode == "wta":
        out["strongest_input"] = {j["name"]: int(np.argmax([dv for c, w, dv in j["regions"]])) for j in jobs}
    return out


def write_once(path, obj):
    if path.exists():
        print("zaten var, üzerine yazılmadı:", path)
        return
    path.write_text(json.dumps(obj, indent=1), encoding="utf-8")
    print("yazıldı:", path)


if __name__ == "__main__":
    mode = sys.argv[1]
    RES.mkdir(parents=True, exist_ok=True)
    if mode == "calib":
        cal = calibrate()
        print(json.dumps(cal["table"], indent=1), "\nell_best =", cal["ell_best"], "| births per run:", cal["births_per_run"])
        write_once(RES / "calibration.json", cal)
    else:
        cal = json.loads((RES / "calibration.json").read_text(encoding="utf-8"))
        pr = predict(mode, cal["ell_best"])
        if mode == "place":
            print("MC P(birth in seed):", pr["P_birth_in_seed_MC"])
        print("per-run region hits:", sum(v["region"] is not None for v in pr["per_run"].values()), "/", len(pr["per_run"]))
        write_once(RES / f"predictions_{mode}.json", pr)
