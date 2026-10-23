"""Exp22: which quantum states exist in SKY130 silicon, and what would show them? (from the project's own solvers/data)
1. Inversion layer: self-consistent Poisson-Schroedinger (src/qcell/mos1d.py, SKY130 stack: tox 4.148 nm, N_A 1.7e17)
   at 300 K; subband energies vs n_s. 4 K occupation on the same subbands (the confining field is set by depletion and
   inversion charge; dopants in the depletion layer stay ionised): Fermi level for each n_s, number of occupied subbands,
   n_s where the second subband starts to fill, and the resulting step of the gate capacitance.
2. Landau quantisation (Shubnikov-de Haas / quantum Hall): field needed for omega_c tau > 1 and hbar omega_c > k T, and
   the fields of the integer filling factors at typical n_s.
3. Single-electron charging: E_C of gate-defined islands (exp1 FEM constant) vs temperature.
4. A filament "cell": number of electrons and of body holes -> classical object.
usage: python q22.py  ->  results/exp22/quantum_signatures.json, fig_exp22.png"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.mos1d import Mos, Q, KB, HBAR, M0, VALLEYS  # noqa: E402

OUT = ROOT / "results" / "exp22"
EPS0, EPS_SI = 8.8541878128e-12, 11.7


def levels(sol):
    """flat list of subbands: energy rel. E_F (eV), valley index, 2D DOS (1/(eV m^2)), centroid (nm)"""
    out = []
    for vi, (w, chi2, v) in enumerate(sol["bands"]):
        dos = v["g"] * v["md"] * M0 / (np.pi * HBAR**2) * Q        # spin included in g? g = valley x spin factor of mos1d
        for E, c2 in zip(w, chi2):
            z = np.arange(len(c2)) * 0.1
            out.append(dict(E=float(E), valley=vi, dos=float(dos), zc=float(np.sum(z * c2) * 0.1)))
    return sorted(out, key=lambda x: x["E"])


def fermi_at(levels_, ns_m2, T):
    """E_F (shift relative to the 300 K reference, eV) giving ns at temperature T on fixed subbands"""
    kT = KB * T / Q

    def n_of(ef):
        return sum(l_["dos"] * kT * np.log1p(np.exp((ef - l_["E"]) / kT)) for l_ in levels_)
    lo, hi = -0.5, 0.5
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if n_of(m) < ns_m2 else (lo, m)
    ef = 0.5 * (lo + hi)
    occ = [l_["dos"] * kT * np.log1p(np.exp((ef - l_["E"]) / kT)) / ns_m2 for l_ in levels_]
    dos_ef = sum(l_["dos"] / (1 + np.exp((l_["E"] - ef) / kT)) for l_ in levels_)
    return ef, occ, dos_ef


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mos = Mos(T=300.0)
    rows = []
    for psi_s in np.linspace(0.80, 1.25, 19):
        sol = mos.solve(float(psi_s))
        ns = sol["ns_cm2"]
        if ns < 5e10:
            continue
        lv = levels(sol)
        ns_m2 = ns * 1e4
        e300, occ300, _ = fermi_at(lv, ns_m2, 300.0)
        e4, occ4, dos4 = fermi_at(lv, ns_m2, 4.0)
        n_occ4 = sum(o > 0.01 for o in occ4)
        zc = sum(o * l_["zc"] for o, l_ in zip(occ4, lv))
        cq = Q * dos4                                                  # F/m^2 (dos in 1/(eV m^2) -> * e)
        cox = mos.cox
        c_tot = 1.0 / (1.0 / cox + 1.0 / max(cq, 1e-9) + zc * 1e-9 / (EPS_SI * EPS0))
        rows.append(dict(psi_s=float(psi_s), Vg=float(sol["Vg"]), ns_cm2=float(ns),
                         E0=lv[0]["E"], E1=lv[1]["E"], E2=lv[2]["E"], valleys=[l_["valley"] for l_ in lv[:3]],
                         dE10_meV=(lv[1]["E"] - lv[0]["E"]) * 1e3, ground_fraction_300K=float(occ300[0]),
                         ground_fraction_4K=float(occ4[0]), occupied_subbands_4K=int(n_occ4),
                         EF_minus_E0_4K_meV=(e4 - lv[0]["E"]) * 1e3, zc_nm=float(zc),
                         C_total_uF_cm2=float(c_tot * 100), C_ox_uF_cm2=float(cox * 100), C_Q_uF_cm2=float(cq * 100)))
        print(f"ns={ns:.2e} cm^-2  E1-E0={rows[-1]['dE10_meV']:.1f} meV  ground 300K {occ300[0]:.2f} / 4K {occ4[0]:.2f}  "
              f"subbands@4K {n_occ4}  C {c_tot*100:.3f} uF/cm2", flush=True)
    n2 = next((r_["ns_cm2"] for r_ in rows if r_["occupied_subbands_4K"] >= 2), None)
    # Landau quantisation (unprimed valleys, m* = 0.19 m0)
    m_eff = 0.19 * M0
    hw_per_T_meV = HBAR * Q / m_eff / Q * 1e3
    landau = dict(hbar_omega_c_meV_per_T=hw_per_T_meV, kT_4K_meV=KB * 4 / Q * 1e3,
                  B_for_hw_gt_kT_4K=KB * 4 / Q * 1e3 / hw_per_T_meV,
                  B_for_mu_B_gt_1={f"{mu}": 1.0 / (mu * 1e-4) for mu in (500, 1000, 2000, 5000)},
                  B_filling_nu={f"ns={ns:.0e}": {f"nu={nu}": ns * 1e4 * 6.62607015e-34 / (Q * nu) for nu in (4, 8, 12)}
                                for ns in (5e11, 1e12, 2e12)})
    # single-electron charging, exp1 FEM constant: E_C/k_B = 1171 K nm / a (eps 3.9), ~585 K nm / a (eps_eff 7.8)
    sizes = (10, 30, 150)
    charging = {f"a={a} nm": dict(EC_over_kB_K_eps3p9=1171.0 / a, EC_over_kB_K_eps7p8=585.0 / a,
                                  T_for_EC_10kT_K=dict(eps3p9=117.1 / a, eps7p8=58.5 / a)) for a in sizes}
    # a filament cell: 100 um wide, 0.15 um long channel at n_s ~ 1e12 cm^-2; body: c_b = 1 fF/um, swing ~0.5 V
    cell = dict(electrons_in_channel=1e12 * (100e-4 * 0.15e-4), body_holes_per_um=1e-15 * 0.5 / Q,
                body_holes_100um=1e-15 * 0.5 / Q * 100)
    out = dict(rows=rows, ns_second_subband_4K_cm2=n2, landau=landau, charging=charging, filament_cell=cell)
    (OUT / "quantum_signatures.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("second subband starts filling at 4 K near n_s =", n2)
    print("Landau:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in landau.items() if k != "B_filling_nu"})
    print("cell:", cell)


if __name__ == "__main__":
    main()
