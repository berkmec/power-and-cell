"""Exp18: programmable cells in the exp17 filament model (wide SKY130 nFET, reduced 4 K model, distributed drain ballast).
A deterministic per-cell gate overdrive offset ("prog", +: locally lower V_T) is added to the random V_T mismatch.
Step 1 (modes calib, place, wta, two): where is the first filament born, and can a seed region choose that place?
Step 2 (mode cut): can a weak spot (negative offset) inside a filament choose where it divides?
usage: python run18.py [nproc] [mode]  ->  results/exp18/runs/<name>.json"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
if sys.platform == "win32":
    import ctypes
    from ctypes import wintypes

    _k = ctypes.windll.kernel32
    _k.GetCurrentProcess.restype = wintypes.HANDLE
    _k.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    _k.SetPriorityClass(_k.GetCurrentProcess(), 0x00004000)
import json  # noqa: E402
import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))
sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
import law  # noqa: E402
import spatial as sp  # noqa: E402
from run import split_events  # noqa: E402  (exp17)

OUT = ROOT / "results" / "exp18" / "runs"
W, DY, LAM, R_D, R_L = 600.0, 0.5, 60.0, 10.0, 5e3     # exp17 clean protocol, reference device set
INPUTS = (50.0, 150.0, 250.0, 350.0, 450.0, 550.0)       # wta input regions (centres, um)


def prog_array(regions, W=W, dy=DY):
    """regions: [(centre_um, width_um, offset_V)] -> per-cell overdrive offset."""
    n = int(round(W / dy))
    y = (np.arange(n) + 0.5) * dy
    p = np.zeros(n)
    for c, w, dv in regions:
        p[np.abs(y - c) < w / 2] += dv
    return p


def wta_inputs(seed):
    return [float(v) for v in np.round(np.random.default_rng(seed + 1000).uniform(0.0, 0.030, len(INPUTS)), 4)]


def jobs_calib():
    return [dict(name=f"calib_s{s}", seed=s, regions=[], vdd_max=1.6) for s in range(10, 30)]


def jobs_place():
    return [dict(name=f"place_d{int(round(dv * 1e3)):02d}_s{s}", seed=s, regions=[(150.0, 10.0, dv)], vdd_max=1.6)
            for dv in (0.0, 0.005, 0.010, 0.015, 0.020, 0.030) for s in range(100, 110)]


def jobs_wta():
    return [dict(name=f"wta_s{s}", seed=s, regions=[(c, 10.0, dv) for c, dv in zip(INPUTS, wta_inputs(s))],
                 vdd_max=1.6) for s in range(200, 210)]


def jobs_two():
    return [dict(name=f"two_s{s}", seed=s, regions=[(150.0, 10.0, 0.030), (450.0, 10.0, 0.030)], vdd_max=1.6)
            for s in range(300, 310)]


def run(job):
    t0 = time.time()
    r = sp.simulate(W=W, dy=DY, R_sq=1e12, R_L=R_L, vdd_up=(0.6, job["vdd_max"]), dvdd=0.02, seed=job["seed"], r_D=R_D,
                    contacts="distributed", r_v=R_D * LAM**2, down=False, sigma_scale=1.0, t_relax=1.0,
                    prog=prog_array(job["regions"]), **law.BASE)
    b = next((h for h in r["hist"] if h["n_fil"] > 0), None)
    birth = dict(vdd=b["vdd"], I_uA=b["I_uA"], centres_um=b["centres_um"], widths_um=b["widths_um"]) if b else None
    res = dict(job=job, birth=birth, splits=split_events(r["hist"]), n_max=max(h["n_fil"] for h in r["hist"]),
               hist=r["hist"], wall=time.time() - t0)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{job['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    where = [(round(c - w / 2), round(c + w / 2)) for c, w in zip(birth["centres_um"], birth["widths_um"])] if birth else "-"
    return f"{job['name']}: doğum V_DD={birth['vdd'] if birth else '-'} yer={where} ({res['wall']:.0f} s)"


if __name__ == "__main__":
    mode = sys.argv[2] if len(sys.argv) > 2 else "calib"
    js = {"calib": jobs_calib, "place": jobs_place, "wta": jobs_wta, "two": jobs_two}[mode]()
    todo = [j for j in js if not (OUT / f"{j['name']}.json").exists()]      # resumable after an interruption
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    print(f"{mode}: {len(todo)} / {len(js)} koşu", flush=True)
    with Pool(nproc) as p:
        for line in p.imap_unordered(run, todo):
            print(line, flush=True)
    print("TAMAM", flush=True)
