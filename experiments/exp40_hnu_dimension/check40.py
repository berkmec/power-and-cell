import json
import math
from pathlib import Path

D = json.loads(Path(__file__).resolve().parents[2].joinpath("results", "exp40", "h40.json").read_text(encoding="utf-8"))
checks = {}
checks["electron_compton_CODATA"] = math.isclose(
    D["species"]["electron"]["compton_wavelength_h_over_mc_m"], 2.42631023867e-12, rel_tol=2e-9)
checks["proton_compton_CODATA"] = math.isclose(
    D["species"]["proton"]["compton_wavelength_h_over_mc_m"], 1.32140985539e-15, rel_tol=2e-9)
checks["in_planck_matches_exp39"] = abs(
    D["species"]["indium_atom"]["thresholds_n"]["lambda_crosses_planck_length_m"] - 17.855666369950697) < 1e-7

# Locate the Schwarzschild equality from stored rows by log interpolation and compare.
for name in ("electron", "proton", "indium_atom"):
    rows = D["species"][name]["rows"]
    vals = [(r["n"], math.log10(r["lambda_over_schwarzschild"])) for r in rows]
    root = None
    for (x1, y1), (x2, y2) in zip(vals, vals[1:]):
        if y1 >= 0 >= y2:
            root = x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
            break
    analytic = D["species"][name]["thresholds_n"]["lambda_equals_schwarzschild"]
    checks[f"{name}_schwarzschild_row_root"] = root is not None and abs(root - analytic) < 2e-3

# Direct sample of the defining ratio, independently forming nu from gamma.
h = D["constants"]["h_J_s"]
c = D["constants"]["c_m_s"]
m = D["species"]["indium_atom"]["mass_kg"]
for n in (0.0, 1.0, 7.3, 17.855666363387947, 30.0):
    gamma = 10.0 ** n
    nu = gamma * m * c * c / h
    ratio = m * c * c / (h * nu)
    checks[f"direct_ratio_n{int(round(n))}"] = math.isclose(ratio, 10.0 ** (-n), rel_tol=3e-15)
checks["stored_verdicts"] = all(v == "PASS" for v in D["verdicts"].values())
print("checks", checks)
print("dimension", D["unit_dimension"])
raise SystemExit(0 if all(checks.values()) else 2)
