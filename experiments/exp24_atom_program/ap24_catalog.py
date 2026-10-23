"""Exp24 (catalogue): techniques for programming natural atoms (0/1) with the 4 K silicon cells, each with its governing
equation evaluated with the project's numbers. Criteria: (S) two stable states, (C) external control by the cell,
(A) single-atom addressability, (R) read-out. Literature values are marked 'lit' and cited in the summary.
T2 field ionisation of a bound hole (B acceptor) in the cell's front field:
     saddle field F_c = pi eps E_b^2 / e^3 ; tunnel rate G(F) = (4/tau_a)(F_a/F) exp(-2 F_a/(3F)), F_a = 2E_b/(e a*),
     tau_a = hbar/(2E_b), a* = e^2/(8 pi eps E_b)
T3 nuclear spin of a single donor (31P) selected by the cell's charge-state control + one global RF tone:
     nu(D+) = gamma_n B ; nu(D0, m_s) = |A/2 -/+ gamma_n B| ; selective if |nu(D0) - nu(D+)| >> Rabi frequency
T4 single-atom contact (Schirm 2013) in a branch of the cell current: switches if |I_branch| = w j_on >= I_s (10-50 uA, lit)
T5 O-vacancy chain in SKY130 ReRAM written through a body-gated 1T1R: n = G/G0 = I_C/(V_C G0) (compliance law, lit),
     I_C = EKV saturation current of the write FET at the cell's body voltage; kick dV = C_gdo V_BL / (c_b dy + C_g,w)
T6 Si-H bond at the interface by multiple vibrational excitation: rate ~ (Phi p_in tau_v)^N / tau_v, N = E_bond / hbar w
T7 tunnel load/unload of one atom by the cell-edge reservoir: G(d) = G0 exp(-2d/a*) -> addressing shell thickness
T8 hole spin of one acceptor driven by EDSR (RF line), addressed by the cell's static field through the g-factor shift
usage: python ap24_catalog.py  ->  results/exp24/catalog.json"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp18_programmable"))
import regen18 as rg  # noqa: E402
from run18 import law, sp, DY  # noqa: E402

OUT = ROOT / "results" / "exp24"
Q, HBAR, H, KB, EPS0, M0 = 1.602176634e-19, 1.054571817e-34, 6.62607015e-34, 1.380649e-23, 8.8541878128e-12, 9.1093837015e-31
EPS = 11.7 * EPS0
G0 = 7.748091729e-5
MUB = 9.2740100783e-24


def t2():
    Eb = 0.045 * Q
    a_star = Q**2 / (8 * np.pi * EPS * Eb)
    Fc = np.pi * EPS * Eb**2 / Q**3
    Fa = 2 * Eb / (Q * a_star)
    tau = HBAR / (2 * Eb)
    rate = lambda F: 4 / tau * (Fa / F) * np.exp(-2 * Fa / (3 * F))
    F1 = next(F for F in np.logspace(5, 8, 3000) if rate(F) >= 1.0)       # field for 1 ionisation per second
    front = {f"front {w} um": 0.9 / (w * 1e-6) for w in (0.5, 2.2)}
    return dict(a_star_nm=a_star * 1e9, F_saddle_kV_cm=Fc / 1e5, F_atomic_kV_cm=Fa / 1e5,
                F_for_1_per_s_kV_cm=F1 / 1e5, rate_at={f"{k}: {v/1e5:.0f} kV/cm": float(rate(v)) for k, v in front.items()},
                depletion_field_kV_cm=0.9 / 83e-9 / 1e5)


def t3():
    A, gn = 117.53e6, 17.235e6                     # Hz, Hz/T (31P, lit)
    out = {}
    for B in (0.5, 1.0, 1.5):
        nu_p = gn * B
        nu_dn, nu_up = abs(A / 2 + gn * B), abs(A / 2 - gn * B)
        out[f"B={B} T"] = dict(nu_ionized_MHz=nu_p / 1e6, nu_neutral_MHz=[nu_up / 1e6, nu_dn / 1e6],
                              min_separation_MHz=min(abs(nu_p - nu_up), abs(nu_p - nu_dn)) / 1e6,
                              electron_ESR_GHz=27.97 * B)
    return dict(levels=out, T2n_ionized_28Si_s=35.6, note="T2n from Muhonen 2014 (lit); Rabi ~ kHz << MHz separation")


def t4():
    j_on = law.constants()["j_on"]
    return dict(I_s_uA=[10, 50], V_s_mV=[50, 200], dG_G0=[0.2, 1.5], j_on_uA_um=j_on * 1e6,
                branch_width_um=[10e-6 / j_on, 50e-6 / j_on], bipolar="two cells (forward / reverse branch)")


def ekv(Vg, VT, W_um, Vds=1.0):
    n, U, mu, cox, L = 1.2, 7.25e-3, 0.06, 8.325e-3, 0.15e-6
    ispec = 2 * n * mu * cox * (W_um * 1e-6 / L) * U**2
    vp = (Vg - VT) / n
    f = lambda v: np.logaddexp(0.0, v / 2.0) ** 2
    return ispec * (f(vp / U) - f((vp - Vds) / U))


def t5():
    g = json.loads((ROOT / "results" / "exp18" / "regen" / "grow_v3.3.json").read_text(encoding="utf-8"))
    Vb = np.array(g["V_final"])
    cells = Vb > 0.5
    out = {}
    for W_w, VT in ((1.0, 0.60), (0.42, 0.70)):
        rng = np.random.default_rng(7)
        vt = VT + rng.normal(0, 3.356e-3 / np.sqrt(W_w * 0.15), len(Vb))
        Ic = ekv(Vb, vt, W_w)
        R_hrs, V_set, V_C = 1e5, 1.46, 0.4
        written = Ic >= V_set / R_hrs                  # the FET can drive the HRS ReRAM to its set voltage
        n = np.where(written, Ic / (V_C * G0), (1 / R_hrs) / G0)
        c_g = 8.325e-3 * W_w * 1e-6 * 0.15e-6 + 2 * 2.449e-10 * W_w * 1e-6
        kick = 2.449e-10 * W_w * 1e-6 * 2.4 / (1e-15 * DY + c_g)
        out[f"W={W_w} um, V_T={VT} V"] = dict(
            fidelity=float(np.mean(written == cells)), bits_1=int(written.sum()), cell_segments=int(cells.sum()),
            n_channels_under_cells=dict(mean=float(n[cells].mean()), p10=float(np.percentile(n[cells], 10)),
                                        p90=float(np.percentile(n[cells], 90))),
            kick_V_at_2p4V=float(kick), low_state_margin_V=0.42)
    return out


def t6():
    phi = 7.2e9                                     # electrons per second past one interface site (exp23 numbers)
    N = int(np.ceil(3.6 / 0.26))
    out = {}
    for tau_v in (1e-9, 1e-8):
        p_crit = (1e-5 * tau_v) ** (1 / N) / (phi * tau_v)          # p_in giving one desorption per ~day
        out[f"tau_v={tau_v:g} s"] = dict(N_quanta=N, p_in_needed_for_1_per_day=float(p_crit))
    return out


def t7():
    a_star = Q**2 / (8 * np.pi * EPS * 0.045 * Q)
    G0r = 1e12
    d = {f"rate {r:g}/s": a_star / 2 * np.log(G0r / r) * 1e9 for r in (1e3, 1.0, 1e-3)}
    return dict(a_star_nm=a_star * 1e9, shell_depth_nm=d, acceptor_spacing_nm=18.05,
                shell_thickness_nm_rate_1e3_to_1em3=(d["rate 0.001/s"] - d["rate 1000/s"]))


def t8():
    out = {}
    for B in (0.3, 1.0):
        fL = 2.0 * MUB * B / H
        out[f"B={B} T"] = {"larmor_GHz": fL / 1e9, "shift_MHz_dg_1e-3": 1e-3 * MUB * B / H / 1e6,
                           "shift_MHz_dg_1e-2": 1e-2 * MUB * B / H / 1e6}
    return dict(levels=out, rabi_lit_MHz=[3, 40], note="Rabi from Crippa 2018 hole qubit (lit); dg range assumed")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = dict(T2_field_ionisation=t2(), T3_nuclear_spin_31P=t3(), T4_single_atom_contact=t4(), T5_reram_chain=t5(),
               T6_SiH_MVE=t6(), T7_tunnel_shell=t7(), T8_acceptor_EDSR=t8())
    (OUT / "catalog.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
