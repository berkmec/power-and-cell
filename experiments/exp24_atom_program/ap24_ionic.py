"""Exp24 (ionic cell): a programmable current-filament cell that CAN be one atom wide, written by the silicon cell.
The SKY130 ReRAM (HfO2, metal1-metal2) holds a conductive filament made of oxygen vacancies. Model: 2D lattice cross-section
(Nx x Nz sites, a = 0.3 nm, 4.8 nm thick), site = O or vacancy V_O; bottom electrode below z = 0, top electrode above z = Nz-1.
  conduction: spanning vacancy cluster -> G = n G0, n = narrowest cross-section of the cluster (one vacancy column = one
              channel, Long et al. APL 2013); otherwise tunnelling through a gap of d layers: G = w_tip G0 exp(-kappa d)
  events    : an O site next to the cluster becomes V_O (O2- drifts to the top reservoir) with rate
              nu exp(-max(0, E_a - alpha F_loc) / (k T_loc)); F_loc = gap field V_R/((d + 1/2) a) at the tip, constriction
              field V_R/(L_c a) at the narrowest layers of a connected filament, body field V_R/(Nz a) elsewhere;
              T_loc = T_amb + R_th G V_R^2 (Joule heating)
  circuit   : bit line ramp V_W -> ReRAM -> write nFET (EKV, 4 K) whose gate is the silicon cell's body voltage V_b ->
              the cell sets the compliance current I_C and so the final width.
Calibration: alpha from the verified 4 K SET voltage 1.46 V (Lu et al. 2026) for a 2-layer gap. Assumptions: kappa = 3.15
per layer (one-layer gap with a 3-wide tip = 100 kOhm, a typical HfO2 HRS, not a SkyWater number), nu = 1e13 /s,
E_a = 1.0 eV, R_th = 1e6 K/W (sensitivity runs vary E_a and R_th with alpha re-calibrated).
Expected from the model (derived before the runs): after the gap closes the constriction widens while alpha V_R/(L_c a)
exceeds E_a - k T ln(nu t); with the current held at I_C this stops at V_C = (E_a - k T ln(nu t)) a / alpha, so
n = max(1, I_C / (G0 V_C)).
usage: python ap24_ionic.py  ->  results/exp24/ionic.json (predictions_ionic.json written first)"""
import json
import time
from collections import deque
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[2] / "results" / "exp24"
G0, KB = 7.748091729e-5, 8.617333e-5           # S, eV/K
NX, NZ, A_NM = 9, 16, 0.3
KAPPA, NU, T_AMB = 3.15, 1e13, 4.0
V_SET, GAP0 = 1.46, 2


def ekv(vg, vds, vt=0.70, w_um=1.0):
    n, U = 1.2, 7.25e-3
    ispec = 2 * n * 0.06 * 8.325e-3 * (w_um / 0.15) * U**2
    vp = (vg - vt) / n
    f = lambda v: np.logaddexp(0.0, v / 2.0) ** 2
    return ispec * (f(vp / U) - f((vp - vds) / U))


def calibrate_alpha(Ea, Rth):
    """alpha such that the first gap-closing event on a 1 s ramp happens at V_R = 1.46 V at 4 K (2-layer gap, HRS)"""
    G_hrs = 3 * G0 * np.exp(-KAPPA * GAP0)
    T = T_AMB + Rth * G_hrs * V_SET**2
    barrier = KB * T * np.log(NU * 1.0)
    F = V_SET / ((GAP0 + 0.5) * A_NM)
    return (Ea - barrier) / F


def initial():
    s = np.zeros((NZ, NX), int)
    s[: NZ - GAP0, 3:6] = 1                       # formed filament, 3 columns wide, reset gap of 2 layers at the top
    return s


def cluster(s):
    seen = np.zeros_like(s, bool)
    q = deque((0, x) for x in range(NX) if s[0, x])
    for z, x in q:
        seen[z, x] = True
    while q:
        z, x = q.popleft()
        for dz, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            zz, xx = z + dz, x + dx
            if 0 <= zz < NZ and 0 <= xx < NX and s[zz, xx] and not seen[zz, xx]:
                seen[zz, xx] = True
                q.append((zz, xx))
    return seen


def conduct(c):
    widths = c.sum(axis=1)
    if widths[NZ - 1] > 0 and np.all(widths > 0):
        n = int(widths.min())
        return n * G0, True, n, int(np.sum(widths == n)), None
    zs = np.nonzero(widths > 0)[0]
    if len(zs) == 0:
        return 1e-9, False, 0, 0, NZ
    z_tip = int(zs.max())
    d = NZ - 1 - z_tip
    return widths[z_tip] * G0 * np.exp(-KAPPA * d) + 1e-9, False, 0, 0, d


def solve_vr(G, vw, vg):
    lo, hi = 0.0, vw
    for _ in range(60):
        m = 0.5 * (lo + hi)
        if G * m > ekv(vg, vw - m):
            hi = m
        else:
            lo = m
    return 0.5 * (lo + hi)


def trial(vg, Ea, Rth, rng, t_ramp=1.0, vmax=2.4, dt=1e-3):
    alpha = calibrate_alpha(Ea, Rth)
    s = initial()
    changed = 0
    t_set = None
    times = np.arange(0.0, 2 * t_ramp + 1e-12, dt)
    for t in times:
        vw = vmax * (t / t_ramp if t <= t_ramp else max(0.0, 2 - t / t_ramp))
        t_end, tt = t + dt, t
        while True:
            c = cluster(s)
            G, conn, n, lc, d = conduct(c)
            vr = solve_vr(G, vw, vg) if vw > 0 else 0.0
            T = T_AMB + Rth * G * vr**2
            cand, rates = [], []
            widths = c.sum(axis=1)
            nmin = n if conn else None
            for z in range(NZ):
                for x in range(NX):
                    if s[z, x]:
                        continue
                    nb = any(0 <= z + dz < NZ and 0 <= x + dx < NX and c[z + dz, x + dx]
                             for dz, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)))
                    if not nb:
                        continue
                    if not conn:
                        z_tip = NZ - 1 - d
                        F = vr / ((d + 0.5) * A_NM) if z > z_tip else vr / (NZ * A_NM)
                    else:
                        F = vr / (lc * A_NM) if widths[z] == nmin else vr / (NZ * A_NM)
                    cand.append((z, x))
                    rates.append(NU * np.exp(-max(0.0, Ea - alpha * F) / (KB * T)))
            R = float(np.sum(rates)) if rates else 0.0
            if R <= 0:
                break
            step = -np.log(rng.random()) / R
            if tt + step > t_end:
                break
            tt += step
            k = rng.choice(len(cand), p=np.array(rates) / R)
            s[cand[k]] = 1
            changed += 1
            if t_set is None and conduct(cluster(s))[1]:
                t_set = tt
    c = cluster(s)
    G, conn, n, lc, d = conduct(c)
    return dict(connected=bool(conn), n=n, G_over_G0=G / G0, sites_changed=changed,
                V_set=None if t_set is None else float(vmax * t_set / t_ramp), I_C_uA=float(ekv(vg, 1.0) * 1e6))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pf = OUT / "predictions_ionic.json"
    Ea, Rth = 1.0, 1e6
    alpha = calibrate_alpha(Ea, Rth)
    Vc = (Ea - KB * (T_AMB + Rth * 1e-4 * 0.5) * np.log(NU * 0.1)) * A_NM / alpha
    if not pf.exists():
        pf.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), alpha_eV_per_V_per_nm=alpha,
                                      V_C_estimate_V=Vc, predictions=[
            f"IC1: median n over 20 trials = max(1, round(I_C / (G0 V_C))) within +-1 for every body voltage, V_C = {Vc:.2f} V",
            "IC2: for I_C < G0 V_C (single-channel regime) the final constriction is exactly one site (one atom) wide in "
            ">= 80 % of the trials",
            "IC3: atoms moved per write (O sites converted) <= GAP0 + n + 2",
            "IC4: across the cell's own body-voltage window 0.805-1.006 V the median n spans from 1 to >= 3"]), indent=1),
            encoding="utf-8")
    rng = np.random.default_rng(24)
    res = dict(alpha=alpha, V_C_estimate=Vc, runs={})
    for vg in (0.81, 0.85, 0.90, 0.95, 1.00):
        rows = [trial(vg, Ea, Rth, rng) for _ in range(20)]
        ns = [r["n"] for r in rows]
        ic = rows[0]["I_C_uA"]
        res["runs"][f"V_b={vg:.2f}"] = dict(I_C_uA=ic, n=ns, median_n=float(np.median(ns)),
                                            predicted_n=int(max(1, round(ic * 1e-6 / (G0 * Vc)))),
                                            frac_single=float(np.mean(np.array(ns) == 1)),
                                            sites_changed=[r["sites_changed"] for r in rows],
                                            V_set=[r["V_set"] for r in rows], connected=sum(r["connected"] for r in rows))
        x = res["runs"][f"V_b={vg:.2f}"]
        print(f"V_b {vg:.2f} V: I_C {ic:6.1f} uA, n {ns} -> median {x['median_n']}, predicted {x['predicted_n']}, "
              f"single-atom {x['frac_single']:.2f}, set at {np.median([v for v in x['V_set'] if v]) if any(x['V_set']) else None} V, "
              f"sites changed median {np.median(x['sites_changed'])}", flush=True)
    sens = {}
    for Ea_s, Rth_s in ((0.8, 1e6), (1.2, 1e6), (1.0, 1e5), (1.0, 1e7)):
        rows = [trial(0.95, Ea_s, Rth_s, rng) for _ in range(10)]
        sens[f"Ea={Ea_s}, Rth={Rth_s:g}"] = dict(median_n=float(np.median([r["n"] for r in rows])),
                                                 V_set=float(np.median([r["V_set"] for r in rows if r["V_set"]] or [np.nan])))
        print("sensitivity", f"Ea={Ea_s}, Rth={Rth_s:g}", sens[f"Ea={Ea_s}, Rth={Rth_s:g}"], flush=True)
    res["sensitivity_Vb_0.95"] = sens
    (OUT / "ionic.json").write_text(json.dumps(res, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
