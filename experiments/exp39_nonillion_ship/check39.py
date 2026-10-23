import json
import math
from pathlib import Path

p = Path(__file__).resolve().parents[2].joinpath("results", "exp39", "ship39.json")
D = json.loads(p.read_text(encoding="utf-8"))
by = {int(r["n"]): r for r in D["rungs"]}
C = D["constants"]["c_m_s"]
checks = {}
# Independent algebra: beta from gamma, photon ratio from gamma, Planck crossing from constants.
checks["beta_n1"] = math.isclose(math.sqrt(1 - 1/by[1]["gamma"]**2), by[1]["beta"], rel_tol=0, abs_tol=1e-15)
checks["beta_n30_deficit"] = math.isclose(by[30]["one_minus_beta"], 0.5e-60, rel_tol=2e-3)
photon_R_n1 = by[1]["gamma"] + math.sqrt(by[1]["gamma"]**2 - 1)
checks["photon_R_n1"] = math.isclose(10 ** by[1]["mass_ratio_log10"]["photon_lower_bound"], photon_R_n1, rel_tol=1e-12)
checks["planck_cross_bracket"] = by[17]["in_wavelength_over_planck"] > 1 > by[18]["in_wavelength_over_planck"]
# Nonrelativistic cross-check of actual burns: beta*c approximately ve*ln(10/8) to O(beta^3).
for name, b in D["actual_burns"].items():
    approx = b["ve_m_s"] * math.log(10/8) / C
    checks[f"burn_{name}"] = math.isclose(b["final_beta"], approx, rel_tol=2e-9)
checks["stored_verdicts"] = all(v == "PASS" for v in D["verdicts"].values())
print("checks", checks)
print("photon_R_n1", photon_R_n1)
print("planck_cross", D["planck_crossing"])
print("literal", D["literal_n_equals_nonillion"])
raise SystemExit(0 if all(checks.values()) else 2)
