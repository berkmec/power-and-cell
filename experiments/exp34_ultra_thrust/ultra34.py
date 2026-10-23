"""Exp34: ultra high-efficiency ion thrust equation from ALL project data (Exp0-Exp33).
Criteria fixed in results/exp34/predictions_exp34.json BEFORE this grid was computed.
Core (Exp28 S6 / Exp32): F = gamma * I_tot * sqrt(2 m V/q), I_tot = N * I_c, I_c = 0.75 j_on w_s(lambda,r_D);
P = I_tot (V + V_M + E_entry/e) / eta_rec;  F/P = 2 eta / v_e, eta = gamma * V/(V + V_M + E_entry/e) (eta_rec=1 nominal).
Constraints: Child-Langmuir F/A = (8/9) eps0 (V/d)^2; emission supply field F_th per species (Exp24 barrier /
AR4 tunnelling, stored aer28 outputs); geometry floor d >= 4.148 nm (SKY130 tox); cell constants at body-effect
scales s = 1 (Exp32) / 0.50 / 0.33 (Exp33 exploratory sensitivity, stored cell33_sensitivity.json values).
Validations: Exp29 NEXT closure readback + aer28 measured-thruster eta recomputation (UT6).
Output -> results/exp34/thrust34.json (fallback ~/_exp34_out). Stored E: files are read, never overwritten."""
import json
import sys
import time
from pathlib import Path

import numpy as np

t0 = time.time()
_here = Path(__file__).resolve()
ROOT = _here.parents[2] if (_here.parents[2] / "results" / "exp32").exists() else Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
import law  # noqa: E402

E, AMU, EPS0, CLIGHT, G0 = 1.602176634e-19, 1.66053907e-27, 8.8541878128e-12, 2.99792458e8, 9.80665
GAMMA_B = 0.9                      # beam factor (Exp28 cells / Exp32 cells, gamma=0.9)
D_MIN_NM = 4.148                   # SKY130 gate oxide thickness (project datum, Exp28 explore basis)
ALPHA, Q0 = 0.507, 1.0             # Exp24 barrier law
FTH_PRE = Q0 / ALPHA               # 1.9724 V/nm, pre-charged (ionic) species

# species: mass amu, entry energy eV/charge, supply field threshold V/nm + provenance
SPECIES = {
    "EMI+_precharged": dict(m=111.17, Eentry=1.0, Fth=FTH_PRE, Fth_src="Exp24 barrier Q0/alpha"),
    "O-_precharged": dict(m=15.999, Eentry=1.0, Fth=FTH_PRE, Fth_src="Exp24 barrier Q0/alpha"),
    "In+_neutral": dict(m=114.818, Eentry=5.786, Fth=4.186 / ALPHA, Fth_src="In evaporation barrier (L+I-phi)/alpha, aer28 emission"),
    "Cs+_neutral": dict(m=132.905, Eentry=3.894, Fth=3.799305794261914, Fth_src="AR4 WKB field ionisation 1e9/s, nu=1e15 (stored tunnel.json)"),
    "Ga+_neutral": dict(m=69.723, Eentry=5.999, Fth=4.489 / ALPHA, Fth_src="Ga evaporation barrier/alpha, aer28 emission"),
    "Xe+_neutral": dict(m=131.293, Eentry=12.130, Fth=20.888224692063034, Fth_src="AR4 WKB field ionisation 1e9/s, nu=1e15 (stored tunnel.json)"),
}
CONSTSETS = {
    "s1_Exp32": dict(V_M=0.70101199867378, V_down=0.639048, V_up=2.322625, dV_c=0.061964, j_on=6.405893e-05),
    "s0p50_sens": dict(V_M=0.934385, V_down=0.857405, V_up=2.383833, dV_c=0.076980, j_on=6.528232e-06),
    "s0p33_sens": dict(V_M=1.149236, V_down=1.053119, V_up=2.419073, dV_c=0.096117, j_on=8.526427e-07),
}
VGRID = np.logspace(np.log10(0.5), np.log10(5000.0), 400)


def fp_curve(sp, VM):
    C = VM + sp["Eentry"]
    ve = np.sqrt(2 * E * VGRID / (sp["m"] * AMU))
    eta = GAMMA_B * VGRID / (VGRID + C)
    return C, ve, eta, 2 * eta / ve


def opt_block(name, sp, VM):
    C, ve, eta, FP = fp_curve(sp, VM)
    i = int(np.argmax(FP))
    Vb = sp["Fth"] * D_MIN_NM                       # emission-feasible boundary at beta=1, d=d_min
    feas = VGRID >= Vb
    j = int(np.argmax(np.where(feas, FP, -1)))
    ve_star = float(np.sqrt(2 * E * C / (sp["m"] * AMU)))
    return dict(species=name, C_V=float(C),
                unconstr=dict(V_grid=float(VGRID[i]), V_analytic_star=float(C), FP_uN_W=float(FP[i] * 1e6),
                              FP_analytic_uN_W=float(GAMMA_B / ve_star * 1e6), Isp_s=float(ve[i] / G0)),
                emission_boundary_V=float(Vb), beta_needed_at_Vstar=float(Vb / C),
                constr_beta1=dict(V=float(VGRID[j]), FP_uN_W=float(FP[j] * 1e6), Isp_s=float(ve[j] / G0),
                                  CL_F_per_A_N_m2=float(8 / 9 * EPS0 * (sp["Fth"] * 1e9) ** 2),
                                  P_per_A_W_m2=float((8 / 9 * EPS0 * (sp["Fth"] * 1e9) ** 2) / FP[j])))


def pareto(sp, VM):
    C, ve, eta, FP = fp_curve(sp, VM)
    feas = VGRID >= sp["Fth"] * D_MIN_NM
    pts = []
    for isp_min in (300, 500, 1000, 2000, 3000, 4570, 6000, 9000):
        sel = feas & (ve / G0 >= isp_min)
        if sel.any():
            k = int(np.argmax(np.where(sel, FP, -1)))
            pts.append(dict(Isp_at_least_s=isp_min, V=float(VGRID[k]), Isp_s=float(ve[k] / G0), FP_uN_W=float(FP[k] * 1e6)))
    return pts


def cell_grid(c):
    rows = {}
    for r_D in (5.0, 10.0, 20.0):
        for lam in (40.0, 60.0, 120.0):
            ws = law.w_split(lam, r_D, c)
            rows[f"rD{r_D:g}_lam{lam:g}"] = dict(lambda_c_um=law.lambda_c(r_D, c), w_split_um=ws,
                                                 I_c_mA=0.75 * c["j_on"] * ws * 1e3 if np.isfinite(ws) else None,
                                                 splits=bool(np.isfinite(ws)))
    return rows


def main():
    out = dict(equation="F = gamma*I_tot*sqrt(2mV/q); I_tot=N*I_c; I_c=0.75*j_on*w_s(lambda,r_D); P=I_tot*(V+V_M+E_entry/e)/eta_rec; "
                         "F/P=2*eta/v_e, eta=gamma*V/(V+V_M+E_entry/e); caps: F/A<=(8/9)eps0(V/d)^2, beta*V/d>=F_th(species), "
                         "d>=4.148nm, N integer steps, F/P<= (photon door only via beamed sail 2/c - Exp27)",
               gamma_b=GAMMA_B, eta_rec_nominal=1.0, eta_rec_measured_range_Exp23=[0.985, 1.023])
    # ---- species optima at s=1 (VM enters eta) + VM sensitivity of EMI optimum ----
    vm1 = CONSTSETS["s1_Exp32"]["V_M"]
    opts = {n: opt_block(n, sp, vm1) for n, sp in SPECIES.items()}
    out["optima_s1"] = opts
    out["pareto_EMI_s1"] = pareto(SPECIES["EMI+_precharged"], vm1)
    out["pareto_In_s1"] = pareto(SPECIES["In+_neutral"], vm1)
    out["EMI_optimum_vs_body_scale"] = {k: opt_block("EMI+_precharged", SPECIES["EMI+_precharged"], c["V_M"]) for k, c in CONSTSETS.items()}
    emi, inn = opts["EMI+_precharged"], opts["In+_neutral"]
    ut1 = (abs(emi["unconstr"]["V_grid"] / emi["C_V"] - 1) <= 0.10 and 450 <= emi["unconstr"]["FP_uN_W"] <= 600
           and 150 <= emi["unconstr"]["Isp_s"] <= 200 and 220 <= inn["unconstr"]["FP_uN_W"] <= 320)
    print(f"UT1 EMI V*grid={emi['unconstr']['V_grid']:.3f} C={emi['C_V']:.4f} FPmax={emi['unconstr']['FP_uN_W']:.1f}uN/W Isp={emi['unconstr']['Isp_s']:.0f}s | "
          f"In FPmax={inn['unconstr']['FP_uN_W']:.1f} verdict={'pass' if ut1 else 'FAIL'}")
    cb = emi["constr_beta1"]
    ut2_main = (abs(cb["V"] / 8.1815 - 1) <= 0.15 and 350 <= cb["FP_uN_W"] <= 450 and 350 <= cb["Isp_s"] <= 420
                and inn["emission_boundary_V"] >= 34 and inn["constr_beta1"]["FP_uN_W"] < 200)
    ut2_beta = inn["beta_needed_at_Vstar"] >= 50
    print(f"UT2 EMI constr V={cb['V']:.2f} FP={cb['FP_uN_W']:.1f} Isp={cb['Isp_s']:.0f} | In boundary={inn['emission_boundary_V']:.2f}V "
          f"FP@boundary={inn['constr_beta1']['FP_uN_W']:.1f} beta_needed={inn['beta_needed_at_Vstar']:.2f} "
          f"verdict={'pass' if (ut2_main and ut2_beta) else 'partial (beta clause falsified)' if ut2_main else 'FAIL'}")
    ut3 = (2.5e4 <= cb["CL_F_per_A_N_m2"] <= 3.5e4 and cb["P_per_A_W_m2"] > 1e7
           and 100 / cb["P_per_A_W_m2"] <= 1e-5 and 1000 / cb["P_per_A_W_m2"] <= 1e-4)
    print(f"UT3 CL={cb['CL_F_per_A_N_m2']:.3e} N/m2 P/A={cb['P_per_A_W_m2']:.3e} W/m2 A@100W={100/cb['P_per_A_W_m2']:.2e} m2 verdict={'pass' if ut3 else 'FAIL'}")
    verdicts = dict(UT1="pass" if ut1 else "fail",
                    UT2=("pass" if ut2_beta else "partial: beta>=50 sub-clause falsified") if ut2_main else "fail",
                    UT3="pass" if ut3 else "fail")

    # ---- UT4 missions: reference craft 10 kg, 100 W, 2 kg propellant ----
    def mission_row(sp_name, V, VM):
        sp = SPECIES[sp_name]
        ve = float(np.sqrt(2 * E * V / (sp["m"] * AMU)))
        eta_m = float(GAMMA_B * V / (V + VM + sp["Eentry"]))
        dv = ve * float(np.log(10 / 8))
        burn = lambda eta: 0.5 * 2 * ve**2 / (eta * 100) / 3.156e7
        return dict(species=sp_name, V=V, ve_km_s=ve / 1e3, Isp_s=ve / G0, FP_uN_W=2 * eta_m / ve * 1e6,
                    eta_model=eta_m, dv_km_s=dv / 1e3, burn_y_eta06=float(burn(0.6)), burn_y_etamodel=float(burn(eta_m)))
    ve_need = 10000 / float(np.log(10 / 8))
    V_knee_in = float(SPECIES["In+_neutral"]["m"] * AMU * ve_need**2 / (2 * E))
    missions = [mission_row("EMI+_precharged", cb["V"], vm1), mission_row("In+_neutral", V_knee_in, vm1),
                mission_row("In+_neutral", 2000.0, vm1), mission_row("In+_neutral", 5000.0, vm1),
                mission_row("Xe+_neutral", 1800.0, vm1)]
    knee = next(r for r in out["pareto_In_s1"] if r["Isp_at_least_s"] == 4570)
    ut4_dv = missions[3]["dv_km_s"] >= 15 and missions[3]["burn_y_eta06"] <= 5 and missions[0]["dv_km_s"] <= 2.5
    ut4_knee = 2000 <= knee["Isp_s"] <= 4500
    verdicts["UT4"] = "pass" if (ut4_dv and ut4_knee) else ("partial: knee-Isp sub-clause falsified" if ut4_dv else "fail")
    out["missions_10kg_100W_2kg"] = missions
    out["knee_In"] = dict(ve_need_km_s=ve_need / 1e3, V_knee_V=V_knee_in, pareto_point=knee)
    print(f"UT4 dv5kV={missions[3]['dv_km_s']:.2f}km/s burn06={missions[3]['burn_y_eta06']:.2f}y nano dv={missions[0]['dv_km_s']:.3f}km/s "
          f"knee Isp={knee['Isp_s']:.0f}s V={knee['V']:.0f} FP={knee['FP_uN_W']:.1f} verdict={verdicts['UT4']}")

    # ---- UT5 cell grids at the three body-effect constant sets ----
    grids = {k: cell_grid(c) for k, c in CONSTSETS.items()}
    out["cell_grids"] = grids
    g33 = grids["s0p33_sens"]
    split33 = [k for k, v in g33.items() if v["splits"]]
    r = g33.get("rD20_lam120", {})
    ut5 = (len(split33) == 1 and split33[0] == "rD20_lam120" and 450 <= r["w_split_um"] <= 600 and 0.25 <= r["I_c_mA"] <= 0.45
           and all(v["splits"] for v in grids["s0p50_sens"].values()))
    verdicts["UT5"] = "pass" if ut5 else "fail"
    print(f"UT5 s0.33 splits={split33} ws={r.get('w_split_um')} Ic={r.get('I_c_mA')} s0.50 all_split={all(v['splits'] for v in grids['s0p50_sens'].values())} verdict={verdicts['UT5']}")
    f_per_A = 0.9 * float(np.sqrt(2 * 114.818 * AMU * 5e3 / E))
    out["F_per_cell_5kV_In_uN_lam60_rD10"] = {k: (v["rD10_lam60"]["I_c_mA"] * 1e-3 * f_per_A * 1e6 if v["rD10_lam60"]["I_c_mA"] else None)
                                              for k, v in grids.items()}

    # ---- breathing + beamed envelope rows (Exp28 EX3/EX4 formulas) ----
    v_orb = float(np.sqrt(3.986004418e14 / (6.371e6 + 200e3)))
    rho_O = 2.5e15 * 15.999 * AMU
    breath = {}
    for label, fm in (("all_species", 1.0), ("oxygen_only_rho2e-10", rho_O / 2e-10)):
        ve_r = v_orb * (1 + 2.2 / (2 * 0.4 * fm))
        breath[label] = dict(ve_km_s=ve_r / 1e3, P_per_m2_W=float(0.4 * fm * 2e-10 * v_orb * ve_r**2 / (2 * 0.6)),
                             FP_ref_eta06_uN_W=float(2 * 0.6 / ve_r * 1e6))
    beamed = dict(nano_EMI_FP_per_beamed_W_uN=0.5 * cb["FP_uN_W"], nano_ratio_to_sail=float(0.5 * cb["FP_uN_W"] * 1e-6 / (2 / CLIGHT)),
                  knee_In_FP_per_beamed_W_uN=0.5 * knee["FP_uN_W"], knee_ratio_to_sail=float(0.5 * knee["FP_uN_W"] * 1e-6 / (2 / CLIGHT)))
    out["breathing_200km"] = breath
    out["beamed_eta_pv0p5"] = beamed
    print(f"ENVELOPE breath all ve={breath['all_species']['ve_km_s']:.2f}km/s P={breath['all_species']['P_per_m2_W']:.0f}W/m2 | "
          f"beamed nano x{beamed['nano_ratio_to_sail']:.0f} vs sail, knee x{beamed['knee_ratio_to_sail']:.0f}")

    # ---- UT6 validation readbacks (no regression vs measured anchors) ----
    clos = json.loads((ROOT / "results" / "exp29" / "closure.json").read_text(encoding="utf-8"))
    crit = clos["MC1"]["criteria"]
    rock = json.loads((ROOT / "results" / "exp28" / "rocket.json").read_text(encoding="utf-8"))
    dev_eta = max(abs((mm["F_N"] * mm["v_e_m_s"] / (2 * mm["P_W"])) / mm["eta"] - 1) for mm in rock["measured"])
    ut6 = bool(clos["MC1"]["passed"] and crit["a_frac_within_2sigma"] == 1.0 and abs(crit["c_rms_dev_pct"] - 1.44767) < 0.01 and dev_eta < 1e-9)
    verdicts["UT6"] = "pass" if ut6 else "fail"
    out["validation_UT6"] = dict(NEXT_MC1_stored=dict(a_frac=crit["a_frac_within_2sigma"], mean_pct=crit["b_mean_dev_pct"],
                                                       rms_pct=crit["c_rms_dev_pct"], passed=clos["MC1"]["passed"]),
                                  measured_eta_recompute_max_rel_dev=float(dev_eta),
                                  FP_In_5kV_crosscheck_uN_W=missions[3]["FP_uN_W"], note="aer28 cells F/P at 5 kV = 19.61 uN/W")
    print(f"UT6 NEXT a={crit['a_frac_within_2sigma']} rms%={crit['c_rms_dev_pct']:.4f} eta_dev={dev_eta:.2e} verdict={verdicts['UT6']}")

    out["verdicts"] = verdicts
    out["envelope"] = [
        dict(regime="R1 nano-gap pre-charged (max F/P, max F/A)", V_V=cb["V"], species="EMI+", FP_uN_W=cb["FP_uN_W"], Isp_s=cb["Isp_s"],
             F_per_A_N_m2=cb["CL_F_per_A_N_m2"], P_per_A_W_m2=cb["P_per_A_W_m2"], limit="power/thermal area, emission boundary"),
        dict(regime="R2 cell kV (mission)", V_V=5000.0, species="In+", FP_uN_W=missions[3]["FP_uN_W"], Isp_s=missions[3]["Isp_s"],
             dv_km_s=missions[3]["dv_km_s"], limit="I_c steps (body-effect band), propellant"),
        dict(regime="R3 breathing 200 km (no onboard propellant limit)", ve_km_s=breath["all_species"]["ve_km_s"],
             P_per_m2_W=breath["all_species"]["P_per_m2_W"], limit="intake/drag balance, ambient density assumption"),
        dict(regime="R4 beamed + atoms", FP_per_beamed_W_uN=beamed["nano_EMI_FP_per_beamed_W_uN"],
             ratio_to_sail=beamed["nano_ratio_to_sail"], limit="beam delivery, eta_pv=0.5 assumption"),
    ]
    txt = json.dumps(out, indent=1, ensure_ascii=False, default=float)
    dst = ROOT / "results" / "exp34" / "thrust34.json"
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(txt, encoding="utf-8")
        print("WROTE", dst)
    except Exception as ex:
        alt = Path.home() / "_exp34_out"
        alt.mkdir(parents=True, exist_ok=True)
        (alt / "thrust34.json").write_text(txt, encoding="utf-8")
        print(f"WRITE_FAIL {dst} :: {type(ex).__name__}: {ex} -> WROTE {alt / 'thrust34.json'}")
    print(f"ELAPSED_S={time.time()-t0:.1f}")
    return out


if __name__ == "__main__":
    main()
