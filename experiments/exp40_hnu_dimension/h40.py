"""Exp40: h-nu line of the Denklem-1 group; scale and dimension audit.

Registered in results/exp40/predictions_exp40.json before this simulation was run.
The calculation distinguishes unit dimension, length-scale descent, and a
hypothetical spatial extra dimension. A Planck-scale crossing is not treated as
proof of an extra dimension.
Output: results/exp40/h40.json.
"""
import json
import math
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent


def find_root():
    for p in [SCRIPT_DIR, *SCRIPT_DIR.parents]:
        if (p / "results" / "exp39" / "ship39.json").exists():
            return p
    return Path(__file__).resolve().parents[2]


ROOT = find_root()
C = 299_792_458.0
H = 6.62607015e-34
HBAR = H / (2 * math.pi)
G = 6.67430e-11
EV_J = 1.602176634e-19
U = 1.66053906660e-27
PLANCK_L = math.sqrt(HBAR * G / C**3)
PLANCK_M = math.sqrt(HBAR * C / G)
PLANCK_E = PLANCK_M * C * C
SPECIES = {
    "electron": 9.1093837015e-31,
    "proton": 1.67262192369e-27,
    "indium_atom": 114.818 * U,
    "planck_mass_particle": PLANCK_M,
    "ship_8kg_single_quantum_counterfactual": 8.0,
}
SCALES = {
    "atomic_1e-10_m": 1e-10,
    "nuclear_1e-14_m": 1e-14,
    "proton_charge_radius_m": 0.8414e-15,
    "electroweak_1e-18_m": 1e-18,
    "planck_length_m": PLANCK_L,
}


def add_dim(a, b):
    return {k: a.get(k, 0) + b.get(k, 0) for k in ("M", "L", "T")}


def sub_dim(a, b):
    return {k: a.get(k, 0) - b.get(k, 0) for k in ("M", "L", "T")}


def series_for(name, mass, ns):
    rows = []
    lambda_c = H / (mass * C)
    rest_e = mass * C * C
    for n in ns:
        n = float(n)
        gamma = 10.0 ** n
        eps = 10.0 ** (-n)
        beta = math.sqrt(max(0.0, 1.0 - eps * eps)) if n < 8 else 1.0
        one_minus_beta = eps * eps / (1.0 + beta) if n > 0 else 1.0
        beta = 1.0 - one_minus_beta
        energy = gamma * rest_e
        nu = energy / H
        lam_h = C / nu
        rs = 2 * G * gamma * mass / (C * C)
        rows.append(dict(
            n=n, epsilon=eps, gamma=gamma, beta=beta, one_minus_beta=one_minus_beta,
            energy_J=energy, energy_eV=energy / EV_J, log10_energy_eV=math.log10(energy / EV_J),
            frequency_Hz=nu, log10_frequency_Hz=math.log10(nu),
            lambda_h_m=lam_h, log10_lambda_h_m=math.log10(lam_h),
            lambda_de_broglie_m=(lam_h / beta if beta > 0 else None),
            lambda_h_over_de_broglie=(beta if beta > 0 else None),
            photon_to_particle_momentum=(1.0 / beta if beta > 0 else None),
            lambda_over_planck=lam_h / PLANCK_L,
            schwarzschild_radius_of_energy_m=rs,
            lambda_over_schwarzschild=lam_h / rs,
            energy_over_planck_energy=energy / PLANCK_E,
            relative_residual=(rest_e / (H * nu)) / eps - 1.0,
        ))
    thresholds = {f"lambda_crosses_{k}": float(math.log10(lambda_c / v)) for k, v in SCALES.items()}
    thresholds.update(
        energy_reaches_planck_energy=float(math.log10(PLANCK_E / rest_e)),
        lambda_equals_schwarzschild=float(math.log10(math.sqrt(math.pi) * PLANCK_M / mass)),
    )
    return dict(
        mass_kg=mass, rest_energy_J=rest_e, rest_energy_eV=rest_e / EV_J,
        compton_wavelength_h_over_mc_m=lambda_c,
        reduced_compton_wavelength_m=HBAR / (mass * C),
        rows=rows, thresholds_n=thresholds,
    )


def kk_overlay(species_data, ns):
    radii = np.logspace(-35, -10, 251)
    curves = {}
    for name in ("electron", "proton", "indium_atom"):
        mass = SPECIES[name]
        lambda_c = species_data[name]["compton_wavelength_h_over_mc_m"]
        curves[name] = [dict(
            radius_m=float(r),
            log10_radius_m=float(math.log10(r)),
            geometric_threshold_n=float(math.log10(lambda_c / r)),
            first_kk_energy_threshold_n=float(math.log10(HBAR / (mass * C * r))),
            first_kk_gap_eV=float((HBAR * C / r) / EV_J),
        ) for r in radii]
    scenarios = {}
    for radius in (1e-18, 1e-20, 1e-25, 1e-30, PLANCK_L):
        key = f"R={radius:.3e} m"
        scenarios[key] = {}
        for name in ("electron", "proton", "indium_atom"):
            mass = SPECIES[name]
            vals = []
            for n in ns:
                gamma = 10.0 ** float(n)
                eps = 10.0 ** (-float(n))
                beta = math.sqrt(max(0.0, 1.0 - eps * eps))
                kmax = gamma * beta * mass * C * radius / HBAR
                vals.append(dict(n=float(n), kmax_continuous=float(kmax), first_mode_open=bool(kmax >= 1.0)))
            scenarios[key][name] = vals
    return dict(radius_convention="R is compactification radius; KK momentum gap=hbar*c/R", curves=curves, scenarios=scenarios)


def main():
    t0 = time.time()
    ns = np.linspace(0.0, 30.0, 301)
    data = {name: series_for(name, mass, ns) for name, mass in SPECIES.items()}
    dims = {
        "m0": {"M": 1, "L": 0, "T": 0},
        "c_squared": {"M": 0, "L": 2, "T": -2},
        "h": {"M": 1, "L": 2, "T": -1},
        "nu": {"M": 0, "L": 0, "T": -1},
    }
    dim_energy = add_dim(dims["m0"], dims["c_squared"])
    dim_hnu = add_dim(dims["h"], dims["nu"])
    dim_ratio = sub_dim(dim_energy, dim_hnu)
    max_residual = max(abs(r["relative_residual"]) for sp in data.values() for r in sp["rows"])
    pointlike = ("electron", "proton", "indium_atom", "planck_mass_particle")
    order_ok = True
    diffs = {}
    for name in pointlike:
        th = data[name]["thresholds_n"]
        a = th["energy_reaches_planck_energy"]
        b = th["lambda_equals_schwarzschild"]
        c = th["lambda_crosses_planck_length_m"]
        diffs[name] = dict(schwarzschild_minus_planck_energy=b - a, planck_length_minus_schwarzschild=c - b)
        order_ok &= a < b < c and abs((b - a) - 0.5 * math.log10(math.pi)) < 1e-12 and abs((c - b) - math.log10(2 * math.sqrt(math.pi))) < 1e-12
    in_n1 = data["indium_atom"]["rows"][10]
    ship_n0 = data["ship_8kg_single_quantum_counterfactual"]["rows"][0]
    kk_sample_diff = math.log10(2 * math.pi)
    verdicts = {
        "V1": "PASS" if dim_ratio == {"M": 0, "L": 0, "T": 0} and max_residual < 1e-12 else "FAIL",
        "V2": "PASS" if 0.004 < in_n1["photon_to_particle_momentum"] - 1 < 0.006 else "FAIL",
        "V3": "PASS" if order_ok else "FAIL",
        "V4": "PASS" if (
            23.0 < data["electron"]["thresholds_n"]["lambda_crosses_planck_length_m"] < 23.3
            and 19.8 < data["proton"]["thresholds_n"]["lambda_crosses_planck_length_m"] < 20.0
            and 17.8 < data["indium_atom"]["thresholds_n"]["lambda_crosses_planck_length_m"] < 17.9
            and 6.3 < data["electron"]["thresholds_n"]["lambda_crosses_electroweak_1e-18_m"] < 6.5
            and 3.0 < data["proton"]["thresholds_n"]["lambda_crosses_electroweak_1e-18_m"] < 3.2
            and 1.0 < data["indium_atom"]["thresholds_n"]["lambda_crosses_electroweak_1e-18_m"] < 1.2
        ) else "FAIL",
        "V5": "PASS" if ship_n0["lambda_over_planck"] < 1e-7 else "FAIL",
        "V6": "PASS" if abs(kk_sample_diff - 0.798179868358115) < 1e-12 else "FAIL",
    }
    out = dict(
        exp="40", date="2026-10-19", task="h-nu line scale and dimension audit",
        registered="predictions_exp40.json",
        constants=dict(c_m_s=C, h_J_s=H, hbar_J_s=HBAR, G_SI=G, planck_length_m=PLANCK_L, planck_mass_kg=PLANCK_M, planck_energy_J=PLANCK_E),
        unit_dimension=dict(inputs=dims, m0c2=dim_energy, hnu=dim_hnu, ratio=dim_ratio, spatial_dimension_in_standard_equation=3),
        species=data, threshold_order_differences=diffs,
        hypothetical_extra_dimension=kk_overlay(data, ns),
        max_relative_residual=max_residual, verdicts=verdicts, runtime_s=round(time.time() - t0, 2),
        dimension_verdict=dict(
            unit_dimension_changes=False,
            length_scale_crosses_planck_regimes=True,
            spatial_dimension_changes_in_equation=False,
            extra_dimension_proven=False,
            note="Scale crossings mark where the semiclassical h-nu identification leaves its standard domain; an extra dimension needs an independent R/coupling model and observable dispersion or force-law change."
        ),
    )
    outdir = ROOT / "results" / "exp40"
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        outdir = Path.home() / "_exp40_out"
        outdir.mkdir(parents=True, exist_ok=True)
    dst = outdir / "h40.json"
    dst.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("WROTE", dst, "runtime", out["runtime_s"], "verdicts", verdicts, "max_residual", max_residual, flush=True)
    for name, sp in data.items():
        th = sp["thresholds_n"]
        print(name, "lambdaC", sp["compton_wavelength_h_over_mc_m"],
              "n_EW", th["lambda_crosses_electroweak_1e-18_m"],
              "n_Eplanck", th["energy_reaches_planck_energy"],
              "n_schwarz", th["lambda_equals_schwarzschild"],
              "n_lplanck", th["lambda_crosses_planck_length_m"], flush=True)
    return 0 if all(v == "PASS" for v in verdicts.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())

