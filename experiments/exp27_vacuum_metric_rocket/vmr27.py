"""Exp27 (Vakum-Metrik Roket): an energy-conserving rocket equation for a craft that takes its energy from the recyclable
dividing cells, its thrust from vacuum control (static / dynamic Casimir) and from 'metric engineering' by programmable atoms.
Everything follows from energy-momentum conservation of the closed system craft + emitted fields; every number used is either
measured (references in summary_exp27.md) or a stated model input.
 derive : relativistic rocket that radiates a directed fraction eta of the energy it spends -> beta = tanh(eta ln(M0/M1));
          numerical integration with full lab-frame bookkeeping of craft + radiation energy and momentum
 casimir: static Casimir energy E = -A pi^2 hbar c / (720 a^3) is a state function -> any closed cycle in (a, A) gives zero
          work, and the plate forces are internal -> zero impulse
 dce    : dynamical Casimir radiation. Perfect mirror (EM): f = -hbar A q^(5) / (30 pi^2 c^4) (TE + 11 TM, Dalvit, Maia Neto,
          Mazzitelli review) -> P = hbar A dq^2 W^6 / (60 pi^2 c^4) for q = dq sin(W t); superconducting circuit at the
          measured level of Wilson et al. 2011 (few K per unit bandwidth); the dividing cells as switched plasma mirrors
 metric : gravitational-wave thrust bound (quadrupole formula) and the warp-metric negative-energy requirement versus the
          Casimir negative energy density -pi^2 hbar c / (720 a^4)
 energy : the cell data (Exp23 recycling, Exp25/26 energy per atom), onboard energy sources and the beamed-power door
usage: python vmr27.py derive|casimir|dce|metric|energy|fig  -> results/exp27/*.json"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[2] / "results" / "exp27"
HBAR, C, KB, G = 1.054571817e-34, 2.99792458e8, 1.380649e-23, 6.67430e-11
MSUN = 1.98892e30


def save(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(obj, indent=1, default=float), encoding="utf-8")
    print(json.dumps(obj, indent=1, default=float)[:3000])


def derive():
    """rest-frame steps: the craft turns dE = -c^2 dM into radiation, eta dE collimated backwards, (1-eta) dE isotropic"""
    rows = []
    for eta in (1.0, 0.5, 0.1):
        M0, M, beta = 1.0, 1.0, 0.0
        E_rad, P_rad = 0.0, 0.0                     # lab-frame energy and momentum of all emitted radiation (units c = 1)
        n = 200000
        Mend = 0.1
        dM = (M0 - Mend) / n
        worst_E, worst_P = 0.0, 0.0
        for _ in range(n):
            g = 1 / np.sqrt(1 - beta**2)
            dE = dM                                  # rest-frame radiated energy (c = 1)
            Ed, Ei = eta * dE, (1 - eta) * dE
            E_rad += g * (1 - beta) * Ed + g * Ei    # backward photons: E_lab = g E'(1 - beta); isotropic: g E'
            P_rad += -g * (1 - beta) * Ed + g * beta * Ei
            M -= dM
            dbp = Ed / M                             # rest-frame recoil: M dv' = eta dE
            beta = (beta + dbp) / (1 + beta * dbp)
            g = 1 / np.sqrt(1 - beta**2)
            worst_E = max(worst_E, abs(g * M + E_rad - M0))
            worst_P = max(worst_P, abs(g * beta * M + P_rad))
        exact = np.tanh(eta * np.log(M0 / M))
        rows.append(dict(eta=eta, M0_over_M1=M0 / M, beta_numeric=beta, beta_closed_form=exact, rel_err=abs(beta - exact) / exact,
                         max_energy_violation=worst_E, max_momentum_violation=worst_P))
    out = dict(equation="beta = tanh(eta * ln(M0/M1)); nonrelativistic: dv = eta c ln(M0/M1)", rows=rows,
               thrust_per_power_bound_N_per_W=1 / C)
    save("derive", out)


def derive_exact():
    """second check after the registered first-order run (derive): every step conserves the 4-momentum exactly in the
    craft's rest frame (remaining craft: energy M - dE, momentum eta dE -> new rest mass and recoil), then a Lorentz
    transformation to the lab; energy and momentum are then conserved to rounding, and beta converges to the closed form"""
    rows = []
    for eta in (1.0, 0.5, 0.1):
        for n in (20000, 200000):
            M, beta = 1.0, 0.0
            E_rad, P_rad = 0.0, 0.0
            dE = 0.9 / n
            worst_E = worst_P = 0.0
            for _ in range(n):
                g = 1 / np.sqrt(1 - beta**2)
                Ed, Ei = eta * dE, (1 - eta) * dE
                E_rad += g * (1 - beta) * Ed + g * Ei
                P_rad += -g * (1 - beta) * Ed + g * beta * Ei
                e_rem, p_rem = M - dE, Ed                     # rest-frame energy and momentum of what remains
                M = np.sqrt(e_rem**2 - p_rem**2)
                bp = p_rem / e_rem
                beta = (beta + bp) / (1 + beta * bp)
                g = 1 / np.sqrt(1 - beta**2)
                worst_E = max(worst_E, abs(g * M + E_rad - 1.0))
                worst_P = max(worst_P, abs(g * beta * M + P_rad))
            exact = np.tanh(eta * np.log(1.0 / M))
            rows.append(dict(eta=eta, steps=n, M0_over_M1=1 / M, beta_numeric=beta, beta_closed_form=exact,
                             rel_err=abs(beta - exact) / exact, max_energy_violation=worst_E, max_momentum_violation=worst_P))
    save("derive_exact", dict(rows=rows))


def casimir():
    """closed cycles in (gap a, overlap area A): work = -oint dE = 0; internal forces -> no impulse"""
    e = lambda a: -np.pi**2 * HBAR * C / (720 * a**3)             # energy per area (J/m^2)
    rows = []
    for (a1, a2, A1, A2) in ((10e-9, 100e-9, 1e-12, 4e-12), (1e-9, 5e-9, 1e-6, 2e-6), (20e-9, 30e-9, 1e-4, 1e-4 * 1.5)):
        path = []
        for t in np.linspace(0, 1, 4001):
            path.append((a1 + (a2 - a1) * t, A1))
        for t in np.linspace(0, 1, 4001):
            path.append((a2, A1 + (A2 - A1) * t))
        for t in np.linspace(0, 1, 4001):
            path.append((a2 + (a1 - a2) * t, A2))
        for t in np.linspace(0, 1, 4001):
            path.append((a1, A2 + (A1 - A2) * t))
        path = np.array(path)
        a, A = path[:, 0], path[:, 1]
        Fa = -A * (3 * np.pi**2 * HBAR * C / (720 * a**4))            # -dE/da (attractive, toward smaller a)
        FA = -e(a)                                                    # -dE/dA
        W = np.sum(0.5 * (Fa[1:] + Fa[:-1]) * np.diff(a) + 0.5 * (FA[1:] + FA[:-1]) * np.diff(A))
        scale = np.max(np.abs(A * e(a)))
        # second check (after the registered trapezoid run): exact work along each straight leg = E(start) - E(end)
        corners = [(a1, A1), (a2, A1), (a2, A2), (a1, A2), (a1, A1)]
        W_exact = sum(A_s * e(a_s) - A_e * e(a_e) for (a_s, A_s), (a_e, A_e) in zip(corners[:-1], corners[1:]))
        rows.append(dict(a1_nm=a1 * 1e9, a2_nm=a2 * 1e9, A1_m2=A1, A2_m2=A2, cycle_work_J=W, work_over_energy_scale=abs(W) / scale,
                         exact_leg_work_sum_J=W_exact, exact_over_energy_scale=abs(W_exact) / scale,
                         pressure_at_a1_Pa=float(np.pi**2 * HBAR * C / (240 * a1**4))))
    out = dict(rows=rows, impulse="plate forces are action-reaction pairs inside the craft: net impulse 0 for any cycle",
               conclusion="static Casimir energy is a state function: no net work, no thrust")
    save("casimir", out)


def p_mirror(A, dq, f):
    """EM perfect plane mirror, both sides: P = hbar A dq^2 W^6 / (60 pi^2 c^4)"""
    W = 2 * np.pi * f
    return HBAR * A * dq**2 * W**6 / (60 * np.pi**2 * C**4)


def dce():
    out = {}
    grid = []
    for f in (1e6, 1e8, 1e9, 1e10):
        for dq in (1e-12, 1e-9, 1e-6):
            v = 2 * np.pi * f * dq
            P = p_mirror(1.0, dq, f)
            grid.append(dict(f_Hz=f, dq_m=dq, peak_velocity_m_s=v, v_over_c=v / C, power_W_per_m2=P,
                             one_sided_thrust_N_per_m2=P / C, photons_per_s_per_m2=P / (HBAR * 2 * np.pi * f)))
    out["mechanical_mirror"] = grid
    # superconducting circuit, measured level of Wilson et al. 2011: few K spectral power (photons double the noise)
    for T, B in ((3.0, 50e6), (3.0, 1e9), (10.0, 5e9)):
        P = KB * T * B
        out.setdefault("superconducting_circuit", []).append(dict(noise_temperature_K=T, bandwidth_Hz=B, power_W=P, thrust_N=P / C))
    # dividing cells as switched plasma mirrors: filament width change dq over the switching time (Exp18-20 time scales)
    cells = []
    for dq, t_sw, A in ((10e-6, 1e-3, 10e-6 * 100e-6), (1e-6, 1e-6, 1e-6 * 10e-6), (10e-6, 1e-9, 10e-6 * 100e-6)):
        f = 1 / (2 * t_sw)
        P = p_mirror(A, dq, f)
        cells.append(dict(width_change_m=dq, switching_time_s=t_sw, area_m2=A, power_W=P, thrust_N=P / C))
    out["cells_as_plasma_mirror"] = cells
    out["bound"] = "any DCE thrust <= P_DCE / c, and P_DCE comes from the drive: the vacuum supplies neither energy nor momentum"
    save("dce", out)


def metric():
    out = {}
    # rotating dumbbell (two masses m at distance L/2 from the axis): P = (32/5) G mu^2 r^4 w^6 / c^5 with mu = m/2, r = L
    gw = []
    for m, L, f in ((1.0, 1.0, 1e3), (1e3, 10.0, 1e2), (9e-24, 1e-5, 1e3)):
        w = 2 * np.pi * f
        mu = m / 2
        P = 32 / 5 * G * mu**2 * L**4 * w**6 / C**5
        gw.append(dict(mass_kg=m, size_m=L, rotation_Hz=f, P_GW_W=P, thrust_bound_N=P / C))
    out["gravitational_waves"] = gw
    rho = lambda a: -np.pi**2 * HBAR * C / (720 * a**4)            # Casimir energy density (J/m^3)
    cas = [dict(gap_nm=a * 1e9, energy_density_J_m3=rho(a), mass_density_kg_m3=rho(a) / C**2) for a in (1e-9, 10e-9, 100e-9)]
    out["casimir_negative_energy"] = cas
    need = dict(alcubierre_100m_bubble_kg=6.2e62, van_den_broeck_kg=3 * MSUN)
    have = abs(cas[0]["mass_density_kg_m3"]) * 1.0                 # 1 m^3 filled with 1 nm gaps (plates ignored)
    out["warp_gap_orders"] = {k: float(np.log10(v / have)) for k, v in need.items()}
    out["warp_requirements_kg"] = need
    out["casimir_1m3_1nm_kg"] = -have
    out["momentum"] = ("Bobrick & Martire 2021: every warp drive is a shell moving inertially and needs propulsion; it obeys a "
                       "rocket equation. Asymptotic (ADM) momentum of a closed system changes only by emitted momentum")
    save("metric", out)


def energy():
    out = {}
    out["cells"] = dict(recycling_transfer_eta=[0.985, 1.023], architecture_loss_fraction=[0.75, 0.83],
                        charge_flip_energy_meV=2.8, programmed_29Si_energy_J=[63e-6, 1.5e-3],
                        note="Exp23/25/26: recycling moves energy between cells; it does not create energy")
    sources = dict(li_ion=0.9e6, fission_U235_full_burn=8.2e13, fusion_DT=3.4e14, antimatter_ideal=9.0e16)
    rows = []
    for name, eps in sources.items():
        for f in (0.5, 0.9):
            xi = f * eps / C**2                                     # converted rest-mass fraction
            for eta in (1.0, 0.3):
                beta = np.tanh(eta * np.log(1 / (1 - xi)))
                rows.append(dict(source=name, specific_energy_J_kg=eps, energy_mass_fraction=f, eta=eta, delta_v_m_s=beta * C))
    out["onboard"] = rows
    out["thrust_per_power"] = dict(photon_collimated_N_per_W=1 / C, light_sail_R1_N_per_W=2 / C,
                                   chip_1W_photon_thrust_N=1 / C, cell_array_1kW_N=1e3 / C)
    # thermal photon thrust from the cells' waste heat (architecture loss) with a directional radiator
    out["waste_heat_radiator"] = dict(P_in_W=1.0, heat_fraction=0.8, directionality=[0.0, 0.5, 1.0],
                                      thrust_N=[0.0, 0.5 * 0.8 / C, 0.8 / C])
    out["beamed"] = dict(laser_power_W=100e9, sail_reflectivity=1.0, force_N=2 * 100e9 / C,
                         note="momentum and energy arrive from outside; Delta v is not limited by onboard energy mass")
    save("energy", out)


def fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    D = {k: json.loads((OUT / f"{k}.json").read_text(encoding="utf-8")) for k in ("derive", "dce", "metric", "energy")}
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    a = ax[0]
    R = np.linspace(1, 10, 200)
    for eta in (1.0, 0.5, 0.1):
        a.plot(R, np.tanh(eta * np.log(R)), label=f"η = {eta}")
    for r in D["derive"]["rows"]:
        a.plot(r["M0_over_M1"], r["beta_numeric"], "ko", ms=4)
    a.set(xlabel="M₀ / M₁ (dinlenme kütlesi oranı)", ylabel="β = v / c", title="(a) β = tanh(η ln(M₀/M₁)); noktalar: sayısal")
    a.legend(fontsize=8)
    a = ax[1]
    labels, vals = [], []
    for r in D["dce"]["mechanical_mirror"]:
        if r["dq_m"] == 1e-9:
            labels.append(f"ayna {r['f_Hz']:.0e} Hz, 1 nm")
            vals.append(r["one_sided_thrust_N_per_m2"])
    for r in D["dce"]["superconducting_circuit"]:
        labels.append(f"süperiletken {r['noise_temperature_K']:g} K × {r['bandwidth_Hz']:.0e} Hz")
        vals.append(r["thrust_N"])
    for r in D["dce"]["cells_as_plasma_mirror"]:
        labels.append(f"hücre anahtarı {r['switching_time_s']:.0e} s")
        vals.append(r["thrust_N"])
    labels.append("1 W foton roketi")
    vals.append(1 / C)
    a.barh(range(len(vals)), np.log10(np.maximum(vals, 1e-300)))
    a.set_yticks(range(len(vals)), labels, fontsize=7)
    a.set(xlabel="log₁₀ itki (N; ayna için N/m²)", title="(b) vakum ayağı: dinamik Casimir")
    a = ax[2]
    names = ["Casimir 1 m³, 1 nm (|negatif kütle|)", "Van Den Broeck (~3 M☉)", "Alcubierre 100 m balon"]
    vals = [abs(D["metric"]["casimir_1m3_1nm_kg"]), D["metric"]["warp_requirements_kg"]["van_den_broeck_kg"],
            D["metric"]["warp_requirements_kg"]["alcubierre_100m_bubble_kg"]]
    a.barh(range(3), np.log10(vals), color=["tab:green", "tab:orange", "tab:red"])
    a.set_yticks(range(3), names, fontsize=8)
    a.set(xlabel="log₁₀ negatif kütle-enerji (kg)", title="(c) metrik ayağı: elde edilebilen ve gereken")
    fig.tight_layout()
    fig.savefig(OUT / "fig_exp27.png", dpi=130)
    print("saved", OUT / "fig_exp27.png")


if __name__ == "__main__":
    {"derive": derive, "derive_exact": derive_exact, "casimir": casimir, "dce": dce, "metric": metric, "energy": energy,
     "fig": fig}[sys.argv[1]]()
