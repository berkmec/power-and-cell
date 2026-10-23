"""Exp29: closure test of the Exp28 atom + positive-energy equation against independently MEASURED published data (no fit).
 MC1 thrust   F = gamma J_b sqrt(2 m V_b / q) + F_neutral   vs NEXT thrust-stand data, 40 throttle levels
              (Diamant 2010 Table 2: J_b, V_bps, F_meas; Pollard 2010: measured gamma ranges; Diamant: neutralizer 0.41 mN)
 MC2 velocity v = sqrt(2 q V / m)  ->  same field, time of flight t ~ sqrt(m): Ca+ vs CaO+ (Schnitzler 2009)
 MC3 on demand  predetermined number of ions extracted; success fraction vs detector efficiency (Schnitzler 2009)
 MC4 I = e f    silicon quantum-dot pump, 1 GHz (Zhao 2017)
 MC5 programmable atoms gate vacuum emission: electroformed SiOx nanodiode (Wu 2018)
 MC6 cells      no SKY130 measurement -> model level (pre-registered)
Fixed inputs and criteria: results/exp29/predictions_exp29.json (written 14:44 before this script).
usage: python mc29.py [fig]   -> results/exp29/closure.json, fig_exp29.png"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[2] / "results" / "exp29"
AMU, E, ME_U = 1.66053906660e-27, 1.602176634e-19, 5.48579909e-4
M_XE, M_CA40, M_O16 = 131.293, 39.962591, 15.994915
G_LOW, G_HIGH, SIG_G = 0.945, 0.955, 0.01          # Pollard 2010 measured ranges (midpoints), assumed 1% systematic
F_NEUT, SIG_N, SIG_V = 0.41e-3, 0.3e-3, 20.0        # Diamant 2010 neutralizer thrust; discharge-only thrust and V_b offset (assumed)
POLLARD_UTIL = {1.00: 0.87, 1.20: 0.93, 1.60: 0.93, 2.00: 0.93, 2.35: 0.89, 2.70: 0.89, 3.10: 0.89, 3.52: 0.89}


def ion_thrust(J, V, m_u=M_XE):
    return J * np.sqrt(2 * m_u * AMU * V / E)


def mc1():
    T = json.loads((OUT / "next_diamant2010_table2.json").read_text(encoding="utf-8"))["rows"]
    rows = []
    for r in T:
        g = G_LOW if r["TL"] <= 22 else G_HIGH
        Fi = ion_thrust(r["Jb_A"], r["Vbps_V"])
        Fp = g * Fi + F_NEUT
        Fm = r["F_meas_mN"] * 1e-3
        d = Fp / Fm - 1
        s = np.sqrt((r["unc_pct"] / 100) ** 2 + SIG_G**2 + (SIG_V / (2 * r["Vbps_V"])) ** 2 + (SIG_N / Fm) ** 2)
        g_eff = (Fm - F_NEUT) / Fi
        util = r["Jb_A"] * M_XE * AMU / E / (r["md_mg_s"] * 1e-6)
        rows.append(dict(TL=r["TL"], Jb_A=r["Jb_A"], Vbps_V=r["Vbps_V"], gamma_used=g, F_meas_mN=Fm * 1e3, F_pred_mN=Fp * 1e3,
                         dev_pct=d * 100, sigma_pct=s * 100, within_2sigma=bool(abs(d) <= 2 * s), gamma_eff=g_eff,
                         d_corrTT10_pct=r["d_corr_pct"], util_discharge=util))
    dev = np.array([x["dev_pct"] for x in rows])
    frac = np.mean([x["within_2sigma"] for x in rows])
    hiV = np.array([x["dev_pct"] for x in rows if x["Vbps_V"] >= 1000])
    util = {}
    for x in rows:
        util.setdefault(x["Jb_A"], x["util_discharge"])
    crit = dict(a_frac_within_2sigma=float(frac), a_pass=bool(frac >= 0.9), b_mean_dev_pct=float(dev.mean()), b_pass=bool(abs(dev.mean()) <= 1.5),
                c_rms_dev_pct=float(np.sqrt(np.mean(dev**2))), c_pass=bool(np.sqrt(np.mean(dev**2)) <= 2.0))
    diag = dict(max_abs_dev_pct=float(np.abs(dev).max()), max_TL=int(rows[int(np.abs(dev).argmax())]["TL"]),
                rms_dev_Vbps_ge_1000_pct=float(np.sqrt(np.mean(hiV**2))), n_Vbps_ge_1000=int(hiV.size),
                gamma_eff_range=[float(min(x["gamma_eff"] for x in rows)), float(max(x["gamma_eff"] for x in rows))],
                corr_dev_vs_Vbps=float(np.corrcoef([x["Vbps_V"] for x in rows], dev)[0, 1]),
                discharge_utilisation_from_Jb_and_md={str(k): [round(v, 4), POLLARD_UTIL[k]] for k, v in util.items()})
    return dict(rows=rows, criteria=crit, diagnostics=diag, passed=bool(crit["a_pass"] and crit["b_pass"] and crit["c_pass"]))


def mc2():
    m1, m2 = M_CA40 - ME_U, M_CA40 + M_O16 - ME_U
    pred = np.sqrt(m2 / m1)
    t1, t2, st = 12.8, 15.3, 0.05
    meas = t2 / t1
    sig = meas * np.sqrt((st / t1) ** 2 + (st / t2) ** 2)
    return dict(pred_ratio=float(pred), meas_ratio=float(meas), sigma=float(sig), dev_sigma=float((meas - pred) / sig),
                dev_pct=float((meas / pred - 1) * 100), passed=bool(abs(meas - pred) <= 2 * sig),
                note="t ∝ √(m/q) at fixed electrostatic geometry; different start positions in the mixed crystal are an extra systematic")


def mc3():
    succ = 123 / 139
    det, lo, hi = 0.87, 0.11, 0.07
    m = (M_CA40 - ME_U) * AMU
    v, dv = 19.47e3, 6.3
    return dict(requested_vs_detected="1→1, 2→2, 4→4 events (Fig. 2e–g)", success_fraction=succ, detector_efficiency=det,
                passed=bool(det - lo <= succ <= det + hi), impulse_bit_Ns=m * v, impulse_bit_rel_spread=dv / v,
                kinetic_energy_eV=0.5 * m * v**2 / E)


def mc4():
    f, I_plateau_pA = 1e9, 160.0
    return dict(f_Hz=f, e_f_pA=E * f * 1e12, reported_plateau_pA=I_plateau_pA, accuracy_ppm=0.27, passed=bool(0.27 <= 1.0),
                source="Zhao et al. (2017), arXiv 1703.04795: plateau = ef within 0.27 ppm, 1 GHz, B = 0")


def mc5():
    d = dict(switch_on_V=6.9, I_emit_uA=3.2, I_cond_uA=103.8, V_at=21.3, eff_max=0.166, V_eff_max=40.0, t_on_ns=50, t_off_ns=130,
             controls="unformed graphene nanogap: no emission; metal nanogap on SiNx: none up to 210 V",
             array="100 diodes: 5 A/cm^2 global, <3% decay and <4% rms noise in 24 h at 5e-6 Pa", mechanism="SiOx -> Si reduction (oxygen moved by field), ruptured filament = MIM tunnel diode")
    d["efficiency_at_21V"] = d["I_emit_uA"] / d["I_cond_uA"]
    d["passed"] = bool(max(d["t_on_ns"], d["t_off_ns"]) <= 1000)
    d["honest_limit"] = "programming is binary (formed+switched = emitter); emission did not depend on the initial resistance state (Fig. 1d)"
    return d


def mc6():
    return dict(passed=None, status="model level: no measurement of SKY130 cells; measured analogue only (n-GaAs impurity-breakdown filaments imaged at 4.2 K, Brandl 1989)")


def neutralizer():
    """Path note from MC5: a programmable-atom (electroformed SiOx) emitter as the on-chip neutralizer of each cell.
    Measured inputs (Wu 2018): array 73.4 uA over 82 um x 18 um at 24.8 V; single diode 3.1 % at 21.3 V, 16.6 % at 40 V."""
    j_arr = 73.4e-6 / (82e-4 * 18e-4)
    cases = {"hücre (Exp23, 0,268 mA, 5 kV In⁺)": (0.268e-3, 5e3), "NEXT TL40 (3,52 A, 1800 V)": (3.52, 1800.0)}
    out = dict(array_density_A_cm2=j_arr)
    for name, (I, V) in cases.items():
        P_beam = I * V
        out[name] = dict(area_cm2=I / j_arr, side_um=np.sqrt(I / j_arr) * 1e4,
                         P_conservative_W=I / 0.031 * 21.3, P_best_W=I / 0.166 * 40.0,
                         overhead_conservative=I / 0.031 * 21.3 / P_beam, overhead_best=I / 0.166 * 40.0 / P_beam)
    out["NEXT_neutralizer_flow_TT10_sccm"] = [2.5, 4.0]
    out["sccm_to_mg_s"] = 2.6867811e19 / 60 * M_XE * AMU * 1e6
    return out


def fig(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = res["MC1"]["rows"]
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.4))
    a = ax[0]
    Fm = np.array([r["F_meas_mN"] for r in R]); Fp = np.array([r["F_pred_mN"] for r in R])
    a.errorbar(Fm, Fp, xerr=Fm * np.array([r["sigma_pct"] for r in R]) / 100 * 0, yerr=Fp * np.array([r["sigma_pct"] for r in R]) / 100,
               fmt="o", ms=4, capsize=2, label="öngörü (serbest parametre yok)")
    lim = [20, 245]
    a.plot(lim, lim, "k--", lw=0.8, label="ölçüm = öngörü")
    a.set(xlabel="ölçülen itki, tezgâh (mN) — Diamant 2010", ylabel="F = γ J_b √(2mV/q) + F_nötr (mN)",
          title="(a) NEXT, 40 işletme noktası: ölçüm vs denklem")
    a.legend(fontsize=8)
    a = ax[1]
    tl = [r["TL"] for r in R]
    a.errorbar(tl, [r["dev_pct"] for r in R], yerr=[2 * r["sigma_pct"] for r in R], fmt="o", ms=4, capsize=2, label="Exp29 denklem (±2σ)")
    a.plot(tl, [-r["d_corrTT10_pct"] for r in R], "s", ms=3, color="gray", label="Diamant: (düzeltilmiş TT10 − ölçülen)/ölçülen")
    a.axhline(0, color="k", lw=0.8)
    a.set(xlabel="işletme noktası (TL)", ylabel="(öngörü − ölçüm)/ölçüm (%)", title="(b) sapma; noktalar 2σ bantla")
    a.legend(fontsize=7)
    a = ax[2]
    a.plot(tl, [r["gamma_eff"] for r in R], "o", ms=4, label="γ_etkin = (F_ölç − F_nötr)/(J_b √(2mV_bps/q))")
    a.fill_between([0.5, 22.5], 0.94, 0.95, color="C1", alpha=0.3, label="Pollard 2010 ölçülen γ (%93 kullanım)")
    a.fill_between([22.5, 40.5], 0.95, 0.96, color="C2", alpha=0.3, label="Pollard 2010 ölçülen γ (%89 kullanım)")
    a.set(xlabel="işletme noktası (TL)", ylabel="itki düzeltme çarpanı", title="(c) itkiden geri çıkan γ ile bağımsız plüm ölçümü")
    a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_exp29.png", dpi=125)
    print("saved", OUT / "fig_exp29.png")


if __name__ == "__main__":
    res = dict(MC1=mc1(), MC2=mc2(), MC3=mc3(), MC4=mc4(), MC5=mc5(), MC6=mc6(), path_neutralizer=neutralizer())
    (OUT / "closure.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float), encoding="utf-8")
    c = res["MC1"]
    print("MC1 criteria", c["criteria"], "passed", c["passed"])
    print("MC1 diagnostics", c["diagnostics"])
    for r in c["rows"]:
        print(f"TL{r['TL']:02d} Jb {r['Jb_A']:.2f} Vbps {r['Vbps_V']:6.0f}  meas {r['F_meas_mN']:6.1f}  pred {r['F_pred_mN']:6.1f}  "
              f"dev {r['dev_pct']:+5.2f}%  2σ {2 * r['sigma_pct']:4.2f}%  {'ok' if r['within_2sigma'] else 'OUT'}  γeff {r['gamma_eff']:.4f}")
    for k in ("MC2", "MC3", "MC4", "MC5", "MC6", "path_neutralizer"):
        print(k, res[k])
    if "fig" in sys.argv:
        fig(res)
