"""Exp24 (spin): program the nuclear spin of ONE natural atom (31P) with a cell + one global RF pulse.
The cell (n-body mirror of the 4 K cell, assumed) sets the charge state of the donors by the voltage address (as for the
acceptors in ap24_atoms.py): target donor ionised (D+), the others neutral (D0). Nuclear frequencies (31P, I = 1/2):
  D+ : f = gamma_n B ;   D0 (electron frozen in m_s = -1/2 or +1/2, gamma_e B >> A): f = A/2 + gamma_n B  or  |A/2 - gamma_n B|
A global RF field at f_RF = gamma_n B with nuclear Rabi frequency Omega rotates the target by pi (bit 0 -> 1). Each
transition is a two-level system in the rotating-wave approximation: P(t) = Omega^2/(Omega^2 + d^2) sin^2(sqrt(Omega^2 + d^2) t/2),
d = 2 pi (f - f_RF); checked against a direct numerical integration of the rotating-frame Schroedinger equation.
Combined single-atom programming fidelity = P(charge configuration) x P(target flip) x (1 - crosstalk).
usage: python ap24_spin.py  ->  results/exp24/spin.json"""
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm

OUT = Path(__file__).resolve().parents[2] / "results" / "exp24"
A, GN = 117.53e6, 17.235e6                       # Hz, Hz/T (31P)


def rabi(omega, det, t):
    w = np.sqrt(omega**2 + det**2)
    return omega**2 / w**2 * np.sin(w * t / 2) ** 2


def numeric(omega, det, t, steps=2000):
    """rotating frame: H = (det/2) sz + (omega/2) sx, piecewise-constant integration (check of the formula)"""
    sx = np.array([[0, 1], [1, 0]], complex)
    sz = np.array([[1, 0], [0, -1]], complex)
    U = expm(-1j * (det / 2 * sz + omega / 2 * sx) * t / steps)
    psi = np.array([1, 0], complex)
    for _ in range(steps):
        psi = U @ psi
    return float(abs(psi[1]) ** 2)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pf = OUT / "predictions_spin.json"
    if not pf.exists():
        pf.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), predictions=[
            "SP1: at B = 1 T and a nuclear Rabi frequency of 10 kHz, the target (ionised) 31P nucleus flips with P >= 0.999 "
            "after a pi pulse while every neutral donor nucleus flips with P <= 1e-6",
            "SP2: the analytic Rabi formula and the numerical integration agree to 1e-6"]), indent=1), encoding="utf-8")
    res = {}
    for B in (0.5, 1.0, 1.5):
        f_rf = GN * B
        others = {"D0, m_s=-1/2": A / 2 + GN * B, "D0, m_s=+1/2": abs(A / 2 - GN * B)}
        for om_khz in (1.0, 10.0, 100.0):
            om = 2 * np.pi * om_khz * 1e3
            t_pi = np.pi / om
            p_t = rabi(om, 0.0, t_pi)
            tt = np.linspace(0, t_pi, 4001)
            cross = {k: float(np.max(rabi(om, 2 * np.pi * (f - f_rf), tt))) for k, f in others.items()}
            row = dict(t_pi_us=t_pi * 1e6, P_target=float(p_t), max_crosstalk=cross,
                       detunings_MHz={k: (f - f_rf) / 1e6 for k, f in others.items()})
            if om_khz == 10.0:
                d = 2 * np.pi * (others["D0, m_s=+1/2"] - f_rf)
                row["numeric_check"] = dict(target=numeric(om, 0.0, t_pi), neutral=numeric(om, d, t_pi, 20000),
                                            formula_target=float(p_t), formula_neutral=float(rabi(om, d, t_pi)))
            res[f"B={B} T, Omega={om_khz:g} kHz"] = row
            print(f"B {B} T, Rabi {om_khz:g} kHz: t_pi {t_pi*1e6:.1f} us, P_target {p_t:.6f}, crosstalk "
                  + ", ".join(f"{k}: {v:.2e}" for k, v in cross.items()))
    (OUT / "spin.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res["B=1.0 T, Omega=10 kHz"]["numeric_check"], indent=1))


if __name__ == "__main__":
    main()
