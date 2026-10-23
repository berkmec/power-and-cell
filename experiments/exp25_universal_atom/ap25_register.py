"""Exp25 (universal atom): a programmable atom that interacts with ALL atoms around it and controls them universally.
Central atom: a 31P donor (charge state set by the cell: D0 = coupled hub, D+ = decoupled; n-body mirror cell assumed).
Neighbours: 29Si nuclei of natural silicon (abundance 4.685 %) on the diamond lattice within 6 nm.
Hyperfine (Fermi contact) from the Kohn-Luttinger state (Koiller et al. 2004: a = 2.5 nm, b = 1.4 nm, k0 = 0.844 (2 pi/a0)):
  A(r) = C [sum_{i=x,y,z} F_i(r) cos(k0 r_i)]^2,  F_i = exp(-sqrt(r_perp,i^2/a^2 + r_i^2/b^2)),
  C fixed by the measured [100] first shell A = 2.981 MHz (Hale & Mieher 1969); KL is ~50 % accurate (Feher 1959).
Nuclear frequencies at B: f_k(m_s = +-1/2) = |nu_n +- A_k/2|, nu_n = 8.465 MHz/T x B. Nucleus k is individually addressable
if f_k(+) is >= Delta_min from every other transition (all nuclei, both electron states, including its own f_k(-)).
Universal control: ESR (electron), selective NMR (each nucleus), electron-conditioned NMR = CNOT (electron -> nucleus k);
nucleus-nucleus gates through the electron. Test: GHZ of the electron and the 4 most strongly coupled addressable nuclei with one
global RF coil (all nuclei see every pulse; exact two-level propagators in each nucleus frame, RWA), fidelity up to local Z
phases, electron dephasing exp(-T/T2e).
usage: python ap25_register.py  ->  results/exp25/register.json (predictions_register.json written first)"""
import json
import time
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[2] / "results" / "exp25"
A0_NM, A_T, B_L = 0.5431, 2.5, 1.4
K0 = 0.844 * 2 * np.pi / A0_NM
GN, ABUND, B_FIELD = 8.465e6, 0.04685, 1.0
OMEGA = 2 * np.pi * 10e3                       # nuclear Rabi (rad/s)
DMIN = 100e3                                   # Hz


def S(r):
    out = 0.0
    for i in range(3):
        perp = np.sqrt(np.sum(np.delete(r, i, axis=-1) ** 2, axis=-1))
        F = np.exp(-np.sqrt(perp**2 / A_T**2 + r[..., i] ** 2 / B_L**2))
        out = out + F * np.cos(K0 * r[..., i])
    return out


C_CAL = 2.981e6 / S(np.array([0.0, 0.0, A0_NM])) ** 2


def lattice(R=6.0):
    basis = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0], [.25, .25, .25], [.25, .75, .75], [.75, .25, .75],
                      [.75, .75, .25]])
    n = int(np.ceil(R / A0_NM)) + 1
    g = np.arange(-n, n + 1)
    cells = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    pts = (cells[:, None, :] + basis[None, :, :]).reshape(-1, 3) * A0_NM
    d = np.linalg.norm(pts, axis=1)
    return pts[(d > 1e-6) & (d <= R)]


SITES = lattice()
A_SITES = C_CAL * S(SITES) ** 2


def addressable(A, nu):
    fp, fm = np.abs(nu + A / 2), np.abs(nu - A / 2)
    allf = np.concatenate([fp, fm])
    ok = []
    for k in range(len(A)):
        d = np.abs(allf - fp[k])
        d[k] = np.inf                                      # itself
        ok.append(bool(np.min(d) >= DMIN))
    return np.array(ok), fp, fm


def ghz(A, sel, nu):
    fp, fm = np.abs(nu + A / 2), np.abs(nu - A / 2)
    tau = np.pi / OMEGA
    amp_p, amp_m = [], []
    for j in sel:                                          # only the register nuclei matter for the GHZ projection
        res = {}
        for sec, f in (("+", fp[j]), ("-", fm[j])):
            psi = np.array([1.0 + 0j, 0.0])
            for p, k in enumerate(sel):
                fpulse = fp[k]
                dl = 2 * np.pi * (f - fpulse)
                w = np.sqrt(OMEGA**2 + dl**2)
                Ud = np.cos(w * tau / 2) * np.eye(2) - 1j * np.sin(w * tau / 2) * np.array([[dl, OMEGA], [OMEGA, -dl]]) / w
                ts, te = p * tau, (p + 1) * tau
                rz = lambda th: np.diag([np.exp(1j * th / 2), np.exp(-1j * th / 2)])
                psi = rz(dl * te) @ Ud @ rz(-dl * ts) @ psi
            res[sec] = psi
        amp_p.append(abs(res["+"][1]) ** 2)
        amp_m.append(abs(res["-"][0]) ** 2)
    Pp, Pm = float(np.prod(amp_p)), float(np.prod(amp_m))
    T = len(sel) * tau
    out = {}
    for t2 in (None, 0.3e-3, 10e-3):
        D = 1.0 if t2 is None else np.exp(-T / t2)
        out["no_dephasing" if t2 is None else f"T2e={t2*1e3:g}ms"] = 0.25 * (Pp + Pm + 2 * D * np.sqrt(Pp * Pm))
    return out, T


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pf = OUT / "predictions_register.json"
    if not pf.exists():
        pf.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), predictions=[
            "RG1: at B = 1 T a P donor in natural Si has on average >= 3 individually addressable 29Si nuclei (Delta_min = "
            "100 kHz) and >= 50 % of the configurations have >= 4",
            "RG2: GHZ of the electron + 4 addressable nuclei with one global RF coil (Rabi 10 kHz) has fidelity >= 0.95 "
            "without electron dephasing in >= 80 % of the configurations that have 4",
            "RG3: with T2e = 0.3 ms (natural Si, Hahn) the GHZ fidelity is < 0.9; with T2e = 10 ms it is >= 0.9",
            "RG4: every 29Si within 2 nm is coupled above 1 kHz (the hub reaches all its near neighbours)"]), indent=1),
            encoding="utf-8")
    nu = GN * B_FIELD
    rng = np.random.default_rng(25)
    rows, fids = [], []
    for cfg in range(50):
        occ = rng.random(len(SITES)) < ABUND
        A = A_SITES[occ]
        r = np.linalg.norm(SITES[occ], axis=1)
        ok, fp, fm = addressable(A, nu)
        near = r <= 2.0
        row = dict(n_29Si=int(occ.sum()), n_addressable=int(ok.sum()), A_addr_kHz=sorted((A[ok] / 1e3).round(1).tolist())[::-1],
                   r_addr_nm=r[ok].round(2).tolist(), near_2nm=int(near.sum()),
                   near_coupled_1kHz=int(np.sum(near & (A > 1e3))), coupled_10kHz_within_6nm=int(np.sum(A > 1e4)))
        if ok.sum() >= 4:
            sel = list(np.argsort(-np.where(ok, A, -1))[:4])
            f, T = ghz(A, sel, nu)
            row.update(ghz=f, ghz_time_us=T * 1e6)
            fids.append(f)
        rows.append(row)
    n_addr = np.array([x["n_addressable"] for x in rows])
    summ = dict(mean_addressable=float(n_addr.mean()), frac_ge4=float(np.mean(n_addr >= 4)),
                ghz_configs=len(fids),
                ghz_no_deph_ge_0p95=float(np.mean([f["no_dephasing"] >= 0.95 for f in fids])) if fids else None,
                ghz_mean={k: float(np.mean([f[k] for f in fids])) for k in fids[0]} if fids else None,
                near_all_coupled=float(np.mean([x["near_coupled_1kHz"] == x["near_2nm"] for x in rows])),
                calibration_C_Hz=float(C_CAL), A_first_100_shell_MHz=2.981,
                stark_note="E3 in energy.json: the cell field shifts the donor hyperfine by MHz (A-gate addressing)")
    print(json.dumps(summ, indent=1))
    for x in rows[:5]:
        print(x["n_addressable"], x["A_addr_kHz"][:6], x.get("ghz"))
    (OUT / "register.json").write_text(json.dumps(dict(summary=summ, configs=rows), indent=1), encoding="utf-8")


def spectator_overlap(A, sel, nu, skip):
    """product over all non-register nuclei of |<phi+|phi->|: their small, electron-state-dependent rotations under the
    register pulses reduce the electron coherence of the GHZ state (not in the registered RG6 metric)"""
    fp, fm = np.abs(nu + A / 2), np.abs(nu - A / 2)
    tau = np.pi / OMEGA
    ov = 1.0
    rz = lambda th: np.diag([np.exp(1j * th / 2), np.exp(-1j * th / 2)])
    for j in range(len(A)):
        if j in skip:
            continue
        states = []
        for f in (fp[j], fm[j]):
            psi = np.array([1.0 + 0j, 0.0])
            for p, k in enumerate(sel):
                dl = 2 * np.pi * (f - fp[k])
                w = np.sqrt(OMEGA**2 + dl**2)
                Ud = np.cos(w * tau / 2) * np.eye(2) - 1j * np.sin(w * tau / 2) * np.array([[dl, OMEGA], [OMEGA, -dl]]) / w
                psi = rz(dl * (p + 1) * tau) @ Ud @ rz(-dl * p * tau) @ psi
            states.append(psi)
        ov *= abs(np.vdot(states[0], states[1]))
    return ov


def spectators():
    """GHZ fidelity with and without the spectator factor on fresh enriched configurations (0.1 % and 0.05 % 29Si)"""
    global DMIN, OMEGA
    nu = GN * B_FIELD
    rng = np.random.default_rng(7525)
    out = {}
    for ab in (0.001, 0.0005):
        for dmin in (30e3, 10e3):
            DMIN, OMEGA = dmin, 2 * np.pi * dmin / 10
            rows = []
            for cfg in range(60):
                occ = rng.random(len(SITES)) < ab
                A = A_SITES[occ]
                ok, _, _ = addressable(A, nu)
                if ok.sum() < 4:
                    continue
                sel = list(np.argsort(-np.where(ok, A, -1))[:4])
                f, T = ghz(A, sel, nu)
                sp_ov = spectator_overlap(A, sel, nu, set(sel))
                Pp_Pm = None
                rows.append(dict(registered=f["no_dephasing"], with_spectators_no_deph=None, spect_overlap=sp_ov, T_us=T * 1e6,
                                 f10=f["T2e=10ms"]))
                # fidelity with spectators: coherence term scaled by the spectator overlap
                base = f["no_dephasing"]
                # recover Pp, Pm from the dephasing-free value is not unique; recompute directly
                fpv, fmv = np.abs(nu + A / 2), np.abs(nu - A / 2)
                rows[-1]["with_spectators_no_deph"] = _ghz_with_factor(A, sel, nu, sp_ov)
            key = f"29Si={ab*100:g}%, Dmin={dmin/1e3:g}kHz"
            out[key] = dict(n=len(rows), registered_mean=float(np.mean([x["registered"] for x in rows])) if rows else None,
                            spectator_mean=float(np.mean([x["with_spectators_no_deph"] for x in rows])) if rows else None,
                            spect_overlap_mean=float(np.mean([x["spect_overlap"] for x in rows])) if rows else None)
            print(key, out[key], flush=True)
    DMIN, OMEGA = 100e3, 2 * np.pi * 10e3
    (OUT / "register_spectators.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def _ghz_with_factor(A, sel, nu, factor):
    fp, fm = np.abs(nu + A / 2), np.abs(nu - A / 2)
    tau = np.pi / OMEGA
    rz = lambda th: np.diag([np.exp(1j * th / 2), np.exp(-1j * th / 2)])
    Pp, Pm = 1.0, 1.0
    for j in sel:
        for sec, f in (("+", fp[j]), ("-", fm[j])):
            psi = np.array([1.0 + 0j, 0.0])
            for p, k in enumerate(sel):
                dl = 2 * np.pi * (f - fp[k])
                w = np.sqrt(OMEGA**2 + dl**2)
                Ud = np.cos(w * tau / 2) * np.eye(2) - 1j * np.sin(w * tau / 2) * np.array([[dl, OMEGA], [OMEGA, -dl]]) / w
                psi = rz(dl * (p + 1) * tau) @ Ud @ rz(-dl * p * tau) @ psi
            if sec == "+":
                Pp *= abs(psi[1]) ** 2
            else:
                Pm *= abs(psi[0]) ** 2
    return float(0.25 * (Pp + Pm + 2 * factor * np.sqrt(Pp * Pm)))


def isotopes():
    """design scan after the natural-Si result (RG1 failed by spectral crowding): 29Si abundance and selectivity"""
    global DMIN, OMEGA
    pf = OUT / "predictions_register2.json"
    if not pf.exists():
        pf.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), basis=(
            "natural Si: 0.82 addressable nuclei per donor (RG1 failed): ~150 29Si within 2.5 nm crowd the spectrum"),
            predictions=[
                "RG5: at 0.1 % 29Si and Delta_min = 30 kHz (Rabi 3 kHz) the mean number of addressable nuclei is >= 3 and "
                ">= 80 % of the 29Si within 3 nm that have A >= 2 Delta_min are individually addressable",
                "RG6: for the 0.1 % configurations with >= 4 addressable nuclei the GHZ (electron + 4 nuclei, one global RF "
                "coil) fidelity is >= 0.95 without electron dephasing in >= 80 % of them, >= 0.9 with T2e = 10 ms and < 0.9 "
                "with T2e = 0.3 ms",
                "RG7 (exploration): addressable fraction versus abundance and Delta_min"]), indent=1), encoding="utf-8")
    nu = GN * B_FIELD
    rng = np.random.default_rng(2525)
    out = {}
    for ab in (0.04685, 0.01, 0.001, 0.0005):
        for dmin in (100e3, 30e3, 10e3):
            DMIN, OMEGA = dmin, 2 * np.pi * dmin / 10
            n_addr, frac_strong, fids = [], [], []
            for cfg in range(50):
                occ = rng.random(len(SITES)) < ab
                A = A_SITES[occ]
                r = np.linalg.norm(SITES[occ], axis=1)
                ok, _, _ = addressable(A, nu)
                strong = (r <= 3.0) & (A >= 2 * dmin)
                n_addr.append(int(ok.sum()))
                if strong.any():
                    frac_strong.append(float(np.mean(ok[strong])))
                if ok.sum() >= 4:
                    sel = list(np.argsort(-np.where(ok, A, -1))[:4])
                    f, T = ghz(A, sel, nu)
                    fids.append(dict(f, T_us=T * 1e6))
            key = f"29Si={ab*100:g}%, Dmin={dmin/1e3:g}kHz"
            out[key] = dict(mean_addressable=float(np.mean(n_addr)), frac_ge4=float(np.mean(np.array(n_addr) >= 4)),
                            frac_strong_within_3nm_addressable=float(np.mean(frac_strong)) if frac_strong else None,
                            ghz_n=len(fids), ghz_mean={k: float(np.mean([f[k] for f in fids])) for k in fids[0]} if fids else None,
                            ghz_no_deph_ge_0p95=float(np.mean([f["no_dephasing"] >= 0.95 for f in fids])) if fids else None)
            print(key, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out[key].items() if k != "ghz_mean"},
                  {k: round(v, 3) for k, v in (out[key]["ghz_mean"] or {}).items()}, flush=True)
    DMIN, OMEGA = 100e3, 2 * np.pi * 10e3
    (OUT / "register_isotopes.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "isotopes":
        isotopes()
    elif len(sys.argv) > 1 and sys.argv[1] == "spectators":
        spectators()
    else:
        main()
