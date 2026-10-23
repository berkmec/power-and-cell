"""Exp32 C3: onset of the second subband at 4 K vs substrate doping, with the project's own Poisson-Schroedinger solver
(src/qcell/mos1d.py, same occupation logic as exp22/q22.py). Measured anchor: Ray & Martin (IBM), (100) p-Si inversion
layers, second subband in Shubnikov-de Haas above ~6e12 cm^-2 (low-doped samples; doping not given in the abstract).
usage: python sub32.py NA_cm3 zmax_nm dz_nm  ->  results/exp32/subband_NA<NA>.json"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "exp22_quantum"))
from qcell.mos1d import Mos  # noqa: E402
from q22 import levels, fermi_at  # noqa: E402

OUT = ROOT / "results" / "exp32"

if __name__ == "__main__":
    na, zmax, dz = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    span = float(sys.argv[4]) if len(sys.argv) > 4 else 0.35
    mos = Mos(na_cm3=na, zmax_nm=zmax, dz_nm=dz, T=300.0)
    rows = []
    phif = float(mos.phif)
    for psi_s in np.linspace(2 * phif - 0.05, 2 * phif + span, 21):
        sol = mos.solve(float(psi_s))
        ns = sol["ns_cm2"]
        if ns < 5e10:
            continue
        lv = levels(sol)
        e4, occ4, _ = fermi_at(lv, ns * 1e4, 4.0)
        n_occ = int(sum(o > 0.01 for o in occ4))
        rows.append(dict(psi_s=float(psi_s), ns_cm2=float(ns), occupied_4K=n_occ, valleys=[l_["valley"] for l_ in lv[:3]],
                         E1_minus_E0_meV=(lv[1]["E"] - lv[0]["E"]) * 1e3, EF_minus_E0_meV=(e4 - lv[0]["E"]) * 1e3))
        print(f"NA {na:.1e}: ns {ns:.2e}  occupied@4K {n_occ}  E1-E0 {rows[-1]['E1_minus_E0_meV']:.1f} meV  "
              f"EF-E0 {rows[-1]['EF_minus_E0_meV']:.1f} meV  next valley {lv[1]['valley']}", flush=True)
    first2 = next((r for r in rows if r["occupied_4K"] >= 2), None)
    last1 = max((r for r in rows if r["occupied_4K"] == 1), key=lambda r: r["ns_cm2"], default=None)
    out = dict(NA_cm3=na, zmax_nm=zmax, dz_nm=dz, rows=rows,
               onset_bracket_cm2=[last1["ns_cm2"] if last1 else None, first2["ns_cm2"] if first2 else None])
    (OUT / f"subband_NA{na:.0e}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("onset bracket", out["onset_bracket_cm2"])
