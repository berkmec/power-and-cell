import json
import math
from pathlib import Path

D = json.loads(Path(__file__).resolve().parents[2].joinpath("results", "exp41", "lab41.json").read_text(encoding="utf-8"))
C = D["constants"]["c_m_s"]
YEAR_S = 365.25 * 24 * 3600
G0 = 9.80665
checks = {}

for name, runs in D["finite_ship"].items():
    lane = D["lab_config"]["finite_ship"]["lanes"][name]
    phi_exact = lane["ve_m_s"] / C * math.log(10 / 8)
    dv_exact = math.tanh(phi_exact) * C / 1000
    checks[f"finite_{name}"] = math.isclose(runs["FULL"]["final"]["delta_v_km_s"], dv_exact, rel_tol=2e-8)
    checks[f"finite_{name}_profile_spread"] = D["finite_summary"][name]["delta_v_spread_fraction"] < 1e-7

for key, run in D["external_bench"].items():
    if run["profile"] != "FULL":
        continue
    final = run["final"]
    a = run["base_g"] * G0
    factor = C / a / YEAR_S
    phi = a * final["proper_years"] * YEAR_S / C
    checks[f"external_{key}_phi"] = math.isclose(final["rapidity"], phi, rel_tol=1e-9, abs_tol=1e-9)
    checks[f"external_{key}_earth"] = math.isclose(final["earth_years"], factor * math.sinh(phi), rel_tol=2e-9)
    checks[f"external_{key}_distance"] = math.isclose(final["distance_ly"], factor * (math.cosh(phi) - 1), rel_tol=2e-9)

ramp_phi = G0 / C * YEAR_S * (85.0 - 0.5)
pulse_phi = 0.5 * G0 / C * YEAR_S * 85.0
checks["ramp_phi_integral"] = math.isclose(D["external_bench"]["1g_RAMP"]["final"]["rapidity"], ramp_phi, rel_tol=1e-7)
checks["pulse_phi_integral"] = math.isclose(D["external_bench"]["1g_PULSE50"]["final"]["rapidity"], pulse_phi, rel_tol=2e-9)
checks["gamma10_beta"] = math.isclose(D["external_events"]["1g_FULL"]["gamma_10"]["beta"], math.sqrt(0.99), rel_tol=1e-15)
checks["stored_verdicts"] = all(v == "PASS" for v in D["verdicts"].values())
print("checks", checks)
raise SystemExit(0 if all(checks.values()) else 2)
