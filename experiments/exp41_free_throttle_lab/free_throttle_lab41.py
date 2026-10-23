"""Exp41: free-speed virtual propulsion laboratory.

No gamma/n value is a control target. Bench A integrates the finite Exp34 ion
ship in Earth time until its 2 kg propellant is exhausted. Bench B integrates an
externally powered 8 kg test article in proper time under prescribed proper
acceleration; its energy bill is recorded and it is not a closed fuel system.
Registered in results/exp41/predictions_exp41.json before this run.
Output: results/exp41/lab41.json.
"""
import json
import math
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent


def find_root():
    for p in [SCRIPT_DIR, *SCRIPT_DIR.parents]:
        if (p / "results" / "exp34" / "thrust34.json").exists():
            return p
    return Path(__file__).resolve().parents[2]


ROOT = find_root()
C = 299_792_458.0
H = 6.62607015e-34
HBAR = H / (2 * math.pi)
G_NEWTON = 6.67430e-11
G0 = 9.80665
E_CHARGE = 1.602176634e-19
U = 1.66053906660e-27
YEAR_S = 365.25 * 24 * 3600
PLANCK_L = math.sqrt(HBAR * G_NEWTON / C**3)
PLANCK_M = math.sqrt(HBAR * C / G_NEWTON)
PLANCK_E = PLANCK_M * C * C
M_IN = 114.818 * U
M_WET = 10.0
M_DRY = 8.0
M_PROP = M_WET - M_DRY
POWER_W = 100.0


def exhaust_lanes():
    thrust = json.loads((ROOT / "results" / "exp34" / "thrust34.json").read_text(encoding="utf-8"))
    emi = thrust["optima_s1"]["EMI+_precharged"]["constr_beta1"]
    in_knee = [x for x in thrust["pareto_In_s1"] if x["Isp_at_least_s"] == 4570][0]
    ve_5kv = math.sqrt(2 * E_CHARGE * 5000.0 / M_IN)
    return {
        "EMI_R1_feasible": dict(ve_m_s=float(emi["Isp_s"] * G0), isp_s=float(emi["Isp_s"]), fp_uN_W=float(emi["FP_uN_W"])),
        "In_knee_R2": dict(ve_m_s=float(in_knee["Isp_s"] * G0), isp_s=float(in_knee["Isp_s"]), fp_uN_W=float(in_knee["FP_uN_W"])),
        "In_5kV_formula": dict(ve_m_s=float(ve_5kv), isp_s=float(ve_5kv / G0), fp_uN_W=19.61),
    }


def finite_throttle(profile, t_years, tfull_years):
    if profile == "FULL":
        return 1.0
    if profile == "ECO25":
        return 0.25
    if profile == "RAMP":
        return min(1.0, t_years / (0.10 * tfull_years))
    if profile == "PULSE50":
        return 1.0 if (t_years / (0.05 * tfull_years)) % 2.0 < 1.0 else 0.0
    raise ValueError(profile)


def simulate_finite(lane, profile):
    u = lane["ve_m_s"]
    thrust_n = lane["fp_uN_W"] * 1e-6 * POWER_W
    mdot_full = thrust_n / u
    tfull_years = M_PROP / mdot_full / YEAR_S
    dt = tfull_years / 4000.0
    # state: mass, rapidity, distance_ly, proper_years, electrical_J
    y = np.array([M_WET, 0.0, 0.0, 0.0, 0.0], dtype=float)
    t = 0.0

    def deriv(tt, yy):
        mass, phi = yy[0], yy[1]
        thr = finite_throttle(profile, tt, tfull_years)
        mdot_year = mdot_full * YEAR_S * thr
        return np.array([
            -mdot_year,
            (u / C) * mdot_year / max(mass, M_DRY),
            math.tanh(phi),
            1.0 / math.cosh(phi),
            thr * POWER_W * YEAR_S,
        ])

    samples = []
    step = 0
    while y[0] > M_DRY + 1e-12 and step < 100000:
        h = dt
        mdot_now = mdot_full * YEAR_S * finite_throttle(profile, t, tfull_years)
        if mdot_now > 0 and y[0] - mdot_now * h < M_DRY:
            h = (y[0] - M_DRY) / mdot_now
        k1 = deriv(t, y)
        k2 = deriv(t + h / 2, y + h * k1 / 2)
        k3 = deriv(t + h / 2, y + h * k2 / 2)
        k4 = deriv(t + h, y + h * k3)
        y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        t += h
        step += 1
        if step % 20 == 0 or y[0] <= M_DRY + 1e-9:
            phi = float(y[1])
            samples.append(dict(
                earth_years=float(t), mass_kg=float(y[0]), rapidity=phi,
                beta=float(math.tanh(phi)), gamma=float(math.cosh(phi)),
                delta_v_km_s=float(math.tanh(phi) * C / 1000.0),
                distance_ly=float(y[2]), proper_years=float(y[3]),
                electrical_energy_J=float(y[4]),
                throttle=float(finite_throttle(profile, t, tfull_years)),
            ))
    phi = float(y[1])
    gamma = math.cosh(phi)
    return dict(
        profile=profile, ve_m_s=float(u), thrust_full_N=float(thrust_n), mdot_full_kg_s=float(mdot_full),
        estimated_full_burn_years=float(tfull_years), steps=step,
        final=dict(earth_years=float(t), mass_kg=float(y[0]), rapidity=phi, beta=float(math.tanh(phi)),
                   gamma=float(gamma), epsilon=float(1.0 / gamma),
                   delta_v_km_s=float(math.tanh(phi) * C / 1000.0),
                   distance_ly=float(y[2]), proper_years=float(y[3]), electrical_energy_J=float(y[4])),
        samples=samples,
    )


def external_throttle(profile, tau_years):
    if profile == "FULL":
        return 1.0
    if profile == "RAMP":
        return min(1.0, tau_years / 1.0)
    if profile == "PULSE50":
        return 1.0 if (tau_years / 0.10) % 2.0 < 1.0 else 0.0
    raise ValueError(profile)


def log10_cosh(phi):
    return (phi - math.log(2.0) + math.log1p(math.exp(-2.0 * phi))) / math.log(10.0) if phi > 1 else math.log10(math.cosh(phi))


def simulate_external(base_g, profile):
    a0 = base_g * G0
    if profile == "FULL":
        tau_max = min(85.0, 80.0 * C / a0 / YEAR_S)
    else:
        tau_max = 85.0
    nsteps = 6000
    dtau = tau_max / nsteps
    # state: rapidity, earth_years, distance_ly
    y = np.array([0.0, 0.0, 0.0], dtype=float)
    tau = 0.0

    def deriv(tt, yy):
        thr = external_throttle(profile, tt)
        phi = float(yy[0])
        return np.array([a0 * thr / C * YEAR_S, math.cosh(phi), math.sinh(phi)])

    samples = []
    for step in range(nsteps + 1):
        if step % 10 == 0 or step == nsteps:
            phi = float(y[0])
            gamma = math.cosh(phi)
            beta = math.tanh(phi)
            thr = external_throttle(profile, tau)
            samples.append(dict(
                proper_years=float(tau), rapidity=phi, beta=float(beta), gamma=float(gamma),
                log10_gamma=float(log10_cosh(phi)), epsilon=float(1.0 / gamma),
                log10_one_minus_beta=float(math.log10(2.0) - (2.0 * phi / math.log(10.0)) if phi > 20 else math.log10(2.0 / (math.exp(2.0 * phi) + 1.0))),
                earth_years=float(y[1]), distance_ly=float(y[2]), throttle=float(thr),
                kinetic_energy_J=float((gamma - 1.0) * M_DRY * C * C),
                instantaneous_external_power_W=float(M_DRY * a0 * thr * C * beta),
                in_lambda_over_planck=float(H / (gamma * M_IN * C) / PLANCK_L),
            ))
        if step == nsteps:
            break
        k1 = deriv(tau, y)
        k2 = deriv(tau + dtau / 2, y + dtau * k1 / 2)
        k3 = deriv(tau + dtau / 2, y + dtau * k2 / 2)
        k4 = deriv(tau + dtau, y + dtau * k3)
        y = y + dtau * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        tau += dtau
    return dict(base_g=base_g, profile=profile, proper_years_max=float(tau_max), samples=samples, final=samples[-1])


def external_events():
    gammas = {
        "gamma_10": 10.0,
        "in_lambda_1e-18_m": H / (M_IN * C) / 1e-18,
        "in_planck_energy": PLANCK_E / (M_IN * C * C),
        "in_lambda_equals_schwarzschild": math.sqrt(math.pi) * PLANCK_M / M_IN,
        "in_lambda_planck_length": H / (M_IN * C) / PLANCK_L,
        "gamma_1e30_readout": 1e30,
    }
    out = {}
    for base_g in (0.1, 1.0, 10.0):
        a = base_g * G0
        factor_years = C / a / YEAR_S
        out[f"{base_g:g}g_FULL"] = {}
        for label, gamma in gammas.items():
            phi = math.acosh(gamma)
            out[f"{base_g:g}g_FULL"][label] = dict(
                gamma=float(gamma), rapidity=float(phi), beta=float(math.tanh(phi)),
                one_minus_beta=float(1.0 / (gamma * gamma) / (1.0 + math.tanh(phi))),
                proper_years=float(phi * factor_years),
                earth_years=float(math.sinh(phi) * factor_years),
                distance_ly=float((gamma - 1.0) * factor_years),
                ship8_kinetic_energy_J=float((gamma - 1.0) * M_DRY * C * C),
            )
    return out


def main():
    t0 = time.time()
    lanes = exhaust_lanes()
    profiles = ("FULL", "ECO25", "RAMP", "PULSE50")
    finite = {name: {p: simulate_finite(lane, p) for p in profiles} for name, lane in lanes.items()}
    external_runs = {}
    for base_g, profile in ((0.1, "FULL"), (1.0, "FULL"), (10.0, "FULL"), (1.0, "RAMP"), (1.0, "PULSE50")):
        external_runs[f"{base_g:g}g_{profile}"] = simulate_external(base_g, profile)
    events = external_events()

    expected_full = {"EMI_R1_feasible": 0.8498651241, "In_knee_R2": 10.0007611922, "In_5kV_formula": 20.4554952725}
    g1 = True
    g2 = True
    finite_summary = {}
    for name, runs in finite.items():
        finals = {p: r["final"] for p, r in runs.items()}
        dvs = np.array([f["delta_v_km_s"] for f in finals.values()])
        ees = np.array([f["electrical_energy_J"] for f in finals.values()])
        spread_dv = float((dvs.max() - dvs.min()) / dvs.mean())
        spread_ee = float((ees.max() - ees.min()) / ees.mean())
        ratios = {p: float(finals[p]["earth_years"] / finals["FULL"]["earth_years"]) for p in profiles}
        finite_summary[name] = dict(delta_v_spread_fraction=spread_dv, electrical_energy_spread_fraction=spread_ee, duration_ratios=ratios)
        g1 &= spread_dv <= 0.005 and abs(finals["FULL"]["delta_v_km_s"] / expected_full[name] - 1) <= 0.01 and finals["FULL"]["gamma"] - 1 <= 2.5e-9
        g2 &= abs(ratios["ECO25"] / 4.0 - 1) <= 0.02 and abs(ratios["PULSE50"] / 2.0 - 1) <= 0.05 and 1.0 <= ratios["RAMP"] <= 1.10 and spread_ee <= 0.02
    e1 = events["1g_FULL"]
    scaling_ok = True
    for label in ("gamma_10", "in_lambda_planck_length", "gamma_1e30_readout"):
        scaling_ok &= abs(events["0.1g_FULL"][label]["proper_years"] / e1[label]["proper_years"] / 10.0 - 1) < 1e-12
        scaling_ok &= abs(events["10g_FULL"][label]["earth_years"] / e1[label]["earth_years"] / 0.1 - 1) < 1e-12
        scaling_ok &= events["0.1g_FULL"][label]["ship8_kinetic_energy_J"] == e1[label]["ship8_kinetic_energy_J"]
    verdicts = {
        "G1": "PASS" if g1 else "FAIL",
        "G2": "PASS" if g2 else "FAIL",
        "G3": "PASS" if 2.8 < e1["gamma_10"]["proper_years"] < 3.0 and 9.5 < e1["gamma_10"]["earth_years"] < 10.5 and 8.5 < e1["gamma_10"]["distance_ly"] < 9.5 and abs(e1["gamma_10"]["ship8_kinetic_energy_J"] / 6.4710372869e18 - 1) < 0.01 else "FAIL",
        "G4": "PASS" if scaling_ok else "FAIL",
        "G5": "PASS" if 65 < e1["gamma_1e30_readout"]["proper_years"] < 70 and e1["gamma_1e30_readout"]["earth_years"] > 1e29 and abs(e1["gamma_1e30_readout"]["ship8_kinetic_energy_J"] / 7.1900414299e47 - 1) < 0.02 and abs(e1["gamma_1e30_readout"]["one_minus_beta"] / 5e-61 - 1) < 0.01 else "FAIL",
        "G6": "PASS" if 39 < e1["in_lambda_planck_length"]["proper_years"] < 42 and 6e17 < e1["in_lambda_planck_length"]["earth_years"] < 8e17 else "FAIL",
    }
    out = dict(
        exp="41", date="2026-10-19", task="free-speed virtual propulsion laboratory", registered="predictions_exp41.json",
        lab_config=dict(
            finite_ship=dict(wet_mass_kg=M_WET, dry_mass_kg=M_DRY, propellant_kg=M_PROP, power_W=POWER_W, profiles=list(profiles), lanes=lanes),
            external_bench=dict(test_mass_kg=M_DRY, runs=list(external_runs), note="externally powered; dry mass constant; energy recorded"),
            readouts_not_targets=["epsilon=1/gamma", "log10_gamma", "in_lambda_over_planck"],
        ),
        constants=dict(c_m_s=C, h_J_s=H, planck_length_m=PLANCK_L, in_mass_kg=M_IN),
        finite_ship=finite, finite_summary=finite_summary,
        external_bench=external_runs, external_events=events,
        verdicts=verdicts, runtime_s=round(time.time() - t0, 2),
    )
    outdir = ROOT / "results" / "exp41"
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        outdir = Path.home() / "_exp41_out"
        outdir.mkdir(parents=True, exist_ok=True)
    dst = outdir / "lab41.json"
    dst.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("WROTE", dst, "runtime", out["runtime_s"], "verdicts", verdicts, flush=True)
    for name, runs in finite.items():
        print("FINITE", name, {p: (round(r["final"]["delta_v_km_s"], 6), round(r["final"]["earth_years"], 6)) for p, r in runs.items()}, finite_summary[name], flush=True)
    for key, ev in e1.items():
        print("EVENT_1G", key, "gamma", ev["gamma"], "proper_y", ev["proper_years"], "earth_y", ev["earth_years"], "ly", ev["distance_ly"], flush=True)
    for key, run in external_runs.items():
        print("EXTERNAL_FINAL", key, run["final"], flush=True)
    return 0 if all(v == "PASS" for v in verdicts.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())


