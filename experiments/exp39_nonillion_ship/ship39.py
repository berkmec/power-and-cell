"""Exp39: Denklem-1 nonillion ladder for a theoretical ion-propelled spaceship.

The three Denklemler images are treated as one Equation-1 group:
sqrt(1-beta^2) = m0*c^2/E = m0*c^2/(h*nu) = 0.1, plus the notes E=-mc^2
and f(r_s). The simulation generalizes 0.1 to epsilon_n=0.1^n=10^-n, so
gamma_n=10^n and the coefficient reaches one short-scale nonillion at n=30.
It does not iterate to nonillion; it works in rapidity/log space.

Registered in results/exp39/predictions_exp39.json before this simulation was run.
Output: results/exp39/ship39.json.
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
E_CHARGE = 1.602176634e-19
U = 1.66053906660e-27
G0 = 9.80665
YEAR_S = 365.25 * 24 * 3600
PLANCK_L = 1.616255e-35
M_IN = 114.818 * U
M_WET = 10.0
M_DRY = 8.0
M_PROP = M_WET - M_DRY
POWER_W = 100.0
ALPHA_CEN_LY = 4.37
NONILLION_LOG10 = 30.0


def exhaust_lanes():
    thrust = json.loads((ROOT / "results" / "exp34" / "thrust34.json").read_text(encoding="utf-8"))
    emi = thrust["optima_s1"]["EMI+_precharged"]["constr_beta1"]
    in_knee = [x for x in thrust["pareto_In_s1"] if x["Isp_at_least_s"] == 4570][0]
    return {
        "EMI_R1_feasible": dict(ve_m_s=float(emi["Isp_s"] * G0), isp_s=float(emi["Isp_s"]), fp_uN_W=float(emi["FP_uN_W"]), source="Exp34 EMI feasible optimum"),
        "In_knee_R2": dict(ve_m_s=float(in_knee["Isp_s"] * G0), isp_s=float(in_knee["Isp_s"]), fp_uN_W=float(in_knee["FP_uN_W"]), source="Exp34 In Pareto knee"),
        "In_5kV_formula": dict(ve_m_s=float(math.sqrt(2 * E_CHARGE * 5000.0 / M_IN)), isp_s=float(math.sqrt(2 * E_CHARGE * 5000.0 / M_IN) / G0), fp_uN_W=19.61, source="Exp34 In 5 kV exhaust-speed formula"),
        "photon_lower_bound": dict(ve_m_s=C, isp_s=C / G0, fp_uN_W=3.3356409519815204e-3, source="photon rocket F/P=1/c in µN/W"),
    }


def ladder_rows(ns, lanes):
    rows = []
    e_rest_dry = M_DRY * C * C
    nu0_in = M_IN * C * C / H
    for n in ns:
        n = float(n)
        gamma = 10.0 ** n
        eps = 10.0 ** (-n)
        beta = math.sqrt(max(0.0, 1.0 - eps * eps)) if n < 8 else 1.0
        one_minus_beta = (eps * eps) / (1.0 + beta) if n > 0 else 1.0
        beta_stable = 1.0 - one_minus_beta
        rapidity = math.acosh(gamma) if gamma > 1 else 0.0
        ekin = (gamma - 1.0) * e_rest_dry
        row = dict(
            n=n, epsilon=eps, gamma=gamma, log10_gamma=n, beta=beta_stable,
            one_minus_beta=one_minus_beta, log10_one_minus_beta=(math.log10(one_minus_beta) if one_minus_beta > 0 else None),
            c_minus_v_m_s=C * one_minus_beta, rapidity=rapidity,
            ship_rest_energy_J=e_rest_dry, ship_total_energy_J=gamma * e_rest_dry,
            ship_kinetic_energy_J=ekin, log10_kinetic_energy_J=(math.log10(ekin) if ekin > 0 else None),
            ideal_100W_supply_years=(ekin / POWER_W / YEAR_S),
            in_ion_frequency_Hz=gamma * nu0_in,
            in_ion_wavelength_m=(H / (gamma * M_IN * C)),
            in_wavelength_over_planck=(H / (gamma * M_IN * C) / PLANCK_L),
            negative_branch_ratio_m0c2_over_Eminus=-eps,
            negative_kinetic_equivalent_kg=-((gamma - 1.0) * M_DRY),
            gravity_analogue_r_over_rs=(1.0 / (1.0 - eps * eps) if eps < 1 else None),
            alpha_centauri_earth_years=(ALPHA_CEN_LY / beta_stable if beta_stable > 0 else None),
            alpha_centauri_ship_proper_years=(ALPHA_CEN_LY / (gamma * beta_stable) if beta_stable > 0 else None),
            mass_ratio_log10={},
        )
        for name, lane in lanes.items():
            row["mass_ratio_log10"][name] = float((C / lane["ve_m_s"]) * rapidity / math.log(10.0))
        rows.append(row)
    return rows


def actual_burn(lanes):
    out = {}
    mass = np.linspace(M_WET, M_DRY, 2001)
    for name, lane in lanes.items():
        if name == "photon_lower_bound":
            continue
        u = lane["ve_m_s"]
        thrust_n = lane["fp_uN_W"] * 1e-6 * POWER_W
        mdot = thrust_n / u
        rapidity = (u / C) * np.log(M_WET / mass)
        beta = np.tanh(rapidity)
        t_s = (M_WET - mass) / mdot
        distance_m = np.concatenate([[0.0], np.cumsum(0.5 * (beta[1:] + beta[:-1]) * C * np.diff(t_s))])
        gamma = np.cosh(rapidity)
        out[name] = dict(
            ve_m_s=float(u), thrust_at_100W_N=float(thrust_n), mdot_kg_s=float(mdot),
            burn_years=float(t_s[-1] / YEAR_S), final_beta=float(beta[-1]), final_gamma=float(gamma[-1]),
            final_epsilon=float(1.0 / gamma[-1]), final_delta_v_km_s=float(beta[-1] * C / 1000.0),
            distance_at_burnout_ly=float(distance_m[-1] / (C * YEAR_S)),
            samples=[dict(mass_kg=float(mass[i]), time_years=float(t_s[i] / YEAR_S), beta=float(beta[i]),
                          gamma=float(gamma[i]), distance_ly=float(distance_m[i] / (C * YEAR_S)))
                     for i in range(0, len(mass), 200)],
        )
    return out


def main():
    t0 = time.time()
    lanes = exhaust_lanes()
    rungs = ladder_rows(range(31), lanes)
    continuous = ladder_rows(np.linspace(0.0, 30.0, 301), lanes)
    burns = actual_burn(lanes)
    by_n = {int(r["n"]): r for r in rungs}
    gamma_planck_cross = H / (M_IN * C * PLANCK_L)
    verdicts = {
        "D1": "PASS" if abs(by_n[1]["beta"] - 0.99498743710662) < 1e-12 and by_n[30]["gamma"] == 1e30 and abs(by_n[30]["one_minus_beta"] / 5e-61 - 1) < 0.01 else "FAIL",
        "D2": "PASS" if all(abs(r["negative_branch_ratio_m0c2_over_Eminus"] + r["epsilon"]) < 1e-18 for r in rungs[1:]) else "FAIL",
        "D3": "PASS" if burns["In_knee_R2"]["final_gamma"] - 1 < 1e-9 and burns["In_knee_R2"]["final_epsilon"] > 0.999999999 else "FAIL",
        "D4": "PASS" if by_n[1]["mass_ratio_log10"]["In_knee_R2"] > 8000 and by_n[1]["mass_ratio_log10"]["In_5kV_formula"] > 4000 and 19 < 10 ** by_n[1]["mass_ratio_log10"]["photon_lower_bound"] < 21 else "FAIL",
        "D5": "PASS" if 6e18 < by_n[1]["ship_kinetic_energy_J"] < 7e18 and by_n[1]["ideal_100W_supply_years"] > 2e9 and 47 < by_n[30]["log10_kinetic_energy_J"] < 48 else "FAIL",
        "D6": "PASS" if by_n[17]["in_wavelength_over_planck"] > 1 and by_n[18]["in_wavelength_over_planck"] < 1 else "FAIL",
        "D7": "PASS" if 0.40 < by_n[1]["alpha_centauri_ship_proper_years"] < 0.48 and 4e-30 < by_n[30]["alpha_centauri_ship_proper_years"] < 5e-30 and abs(by_n[30]["alpha_centauri_earth_years"] - ALPHA_CEN_LY) < 1e-12 else "FAIL",
    }
    literal_n_nonillion = dict(
        interpretation="If n itself is set to one nonillion rather than gamma reaching one nonillion",
        n=1e30, log10_gamma=1e30, log10_epsilon=-1e30,
        rapidity_approx=math.log(2.0) + math.log(10.0) * 1e30,
        log10_one_minus_beta_approx=math.log10(0.5) - 2e30,
        note="Only logarithmic descriptors are finite/representable; beta itself is indistinguishable from 1 in floating point.",
    )
    selected = [by_n[n] for n in (0, 1, 2, 3, 5, 10, 15, 17, 18, 20, 25, 30)]
    out = dict(
        exp="39", date="2026-10-19", task="Denklem-1 nonillion ladder theoretical ion spaceship simulation",
        registered="predictions_exp39.json",
        equation_group=dict(
            main="epsilon_n=sqrt(1-beta^2)=m0*c^2/E=m0*c^2/(h*nu)=0.1^n",
            gamma="gamma_n=1/epsilon_n=10^n",
            nonillion_short_scale=1e30, nonillion_rung_n=30,
            negative_branch_note="E_-=-gamma*m0*c^2 gives m0*c^2/E_-=-10^-n; simultaneous +0.1^n equality requires redefining E as a magnitude or allowing negative frequency/mass.",
            f_rs_note="f(r_s) is undefined in the images; gravity rows are only the Schwarzschild redshift analogue epsilon_g=sqrt(1-r_s/r).",
        ),
        constants=dict(c_m_s=C, h_J_s=H, in_atom_mass_kg=M_IN, planck_length_m=PLANCK_L),
        reference_ship=dict(wet_mass_kg=M_WET, dry_mass_kg=M_DRY, propellant_kg=M_PROP, power_W=POWER_W),
        exhaust_lanes=lanes, actual_burns=burns, rungs=rungs, continuous_ladder=continuous,
        planck_crossing=dict(gamma=gamma_planck_cross, n=math.log10(gamma_planck_cross)),
        literal_n_equals_nonillion=literal_n_nonillion, verdicts=verdicts,
        runtime_s=round(time.time() - t0, 1),
    )
    outdir = ROOT / "results" / "exp39"
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        outdir = Path.home() / "_exp39_out"
        outdir.mkdir(parents=True, exist_ok=True)
    dst = outdir / "ship39.json"
    dst.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("WROTE", dst, "runtime", out["runtime_s"], "verdicts", verdicts, flush=True)
    print("n gamma beta 1-beta Ekin_J log10R_InKnee proper_AC_years lambda_over_Planck")
    for r in selected:
        print(f"{r['n']:.0f} {r['gamma']:.3e} {r['beta']:.15f} {r['one_minus_beta']:.3e} "
              f"{r['ship_kinetic_energy_J']:.3e} {r['mass_ratio_log10']['In_knee_R2']:.1f} "
              f"{r['alpha_centauri_ship_proper_years']} {r['in_wavelength_over_planck']:.3e}", flush=True)
    for name, b in burns.items():
        print("BURN", name, "final_beta", b["final_beta"], "gamma", b["final_gamma"], "burn_years", b["burn_years"], flush=True)
    return 0 if all(v == "PASS" for v in verdicts.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())

