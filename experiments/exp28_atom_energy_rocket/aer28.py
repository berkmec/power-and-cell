"""Exp28 (atom + positive energy rocket from ALL project data): the dividing cells' positive energy accelerates atoms, and the
atoms carry the momentum. Every project number comes from exp0-exp27 (data_registry.json: value, file, role); literature
numbers are listed in summary_exp28.md.
 cells (Exp12-23)          quantised, self-regulating current (I_c = P_c / V_M) and lossless power sharing (eta_rec)
 programmable atoms (Exp24) emission switched by the field-lowered barrier law k = nu exp(-(Q0 - alpha F) / kT)
 field (Exp5/23/25)        the cells' own field; emitters need an enhancement beta
 artificial atoms (Exp0-4) one charge per cycle (Coulomb blockade): atom-by-atom thrust F = m v_e f
 tunnelling (Exp1/4 WKB)   electron leaves the propellant atom into vacuum: exponent b I^1.5 / F
 rocket                    F = mdot v_e, P = mdot v_e^2 / (2 eta)  ->  F / P = 2 eta / v_e  (photons, Exp27: 1 / c)
usage: python aer28.py registry|rocket|cells|emission|tunnel|turnstile|charge|mission|fig"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[2] / "results" / "exp28"
E, AMU, ME, HBAR, KB, EPS0, C, G0 = 1.602176634e-19, 1.66053907e-27, 9.1093837e-31, 1.054571817e-34, 1.380649e-23, \
    8.8541878128e-12, 2.99792458e8, 9.80665
PROP = {"In": dict(m=114.818, I=5.786, phi=4.09, L=2.49), "Cs": dict(m=132.905, I=3.894, phi=2.14, L=0.79),
        "Xe": dict(m=131.293, I=12.130), "Ga": dict(m=69.723, I=5.999, phi=4.32, L=2.81), "O": dict(m=15.999, I=13.618)}
# project data (see registry for files)
V_M, J_ON = 0.80179, 2.569e-6 / 1e-6                     # V; A/m  (Exp17)
P_CELL = (135.9e-6, 220.0e-6)                             # W per cell at 3.3 / 5.0 V (Exp23)
V_D = 0.82                                                # V carrier energy clamp (Exp23)
ETA_REC = (0.985, 1.009, 1.010, 1.023)                    # dying-cell power transfer (Exp23)
ALPHA, Q_OX, NU = 0.507, 1.0, 1e13                        # eV/(V/nm), eV, 1/s (Exp24 ionic KMC)
F_CELL = (0.4e-3, 1.8e-3, 0.02, 0.07, 0.12, 0.038, 0.051)  # V/nm: Exp25 field, Exp23 lateral/vertical, Exp5 surface
EC_KNM = 1171.0                                           # K nm (eps 3.9), Exp1 FEM
B_WKB = 4 / 3 * np.sqrt(2 * ME) * (E ** 1.5) / (HBAR * E) * 1e-9  # V nm^-1 eV^-3/2  (6.83)


def save(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(obj, indent=1, default=float, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(obj, indent=1, default=float, ensure_ascii=False)[:2500])


def ve(sp, V):
    return np.sqrt(2 * E * V / (PROP[sp]["m"] * AMU))


def registry():
    R = [
        ("exp0", "results/exp0/summary.md", "E_C/k_B, dE/k_B vs size: 100 nm 9/5 K, 12.5 nm 75/322 K", "kullanıldı: tek yük adası boyutu"),
        ("exp1", "results/exp1/summary_exp1.md", "c_Sigma = 22.989 eps0 a (FEM, dolfinx 8e-16) -> E_C/k_B = 1171 K nm / a", "kullanıldı: turnike E_C yasası"),
        ("exp2", "results/exp2/chain.json", "12/6/3 nm adalar E_C/k_B 120/232/404 K; 4 K'de yük kuantumu 0,93/0,98/1,00; 77/300 K'de yok",
         "kullanıldı + kısıt: turnike 4 K'de"),
        ("exp3", "results/exp3/dqd.json", "300 K çift nokta: P <= 0,74, kapı <= 0,78", "olumsuz -> kısıt: kuantum denetim kriyojenik"),
        ("exp4", "results/exp4/stage1.json", "WKB kappa 4,4-5,7 /nm, LZ hatası <= 1,2e-3, yükleme ~59 ps (t = 24,8 meV)",
         "kullanıldı: tünelleme yöntemi, turnike frekansı (GHz)"),
        ("exp5", "results/exp5/mos1d.json", "yüzey alanı 0,38-0,51 MV/cm, alt bant 44-50 meV, W_dep 81/55 nm", "kullanıldı: kanal alan ölçeği"),
        ("exp6", "results/exp6", "RDF sigma 20,6 meV; kazara nokta oranı geri çekildi", "olumsuz -> kısıt: tasarlanmış ada gerekir"),
        ("exp7", "experiments/exp7_testplan", ">= 100 aygıt (4 +- 2 % oranı için)", "kullanıldı: dizi istatistiği / ölçüm planı"),
        ("exp8", "results/exp8/summary_exp8.md", "g_c ~ N^-1,08; hazır tek domain bölünmedi", "olumsuz -> ilke: sayı çatallanmayla artar"),
        ("exp9", "results/exp9/COMPARISON.md", "62 satır, 17 kaynak: araçlar uyumlu", "kullanıldı: araç çapraz denetimi"),
        ("exp10", "results/exp10/summary_exp10.md", "2B gaz: K < 0 ama kapı ekranlı Hartree 11-21x|K|: bölünme yok", "olumsuz -> kısıt: filaman hücre"),
        ("exp11", "results/exp11/review_exp6_exp10_exp11.md", "süperkapasitör 0/28, Gray-Scott 5/42", "olumsuz -> ilke: otokataliz (darbe iyonlaşması)"),
        ("exp12", "results/exp12/summary_exp12.md", "V_M = 0,80 V kenetlenmesi; R_box >= 1e10-1e12", "kullanıldı: hücre gerilim kenedi"),
        ("exp13", "results/exp13/summary_exp13.md", "dağıtık balast lambda = 40 um -> ~4 lambda'da bölünme", "kullanıldı: hücre genişliğini balast seçer"),
        ("exp14", "results/exp14_dimension/summary_exp14.md", "bölünme penceresi 2 < d <~ 3,5-4; 2B parmak dizisinde bölünme yok",
         "kısıt: yayıcılar 1B şeritler"),
        ("exp15", "results/exp15_phase_space/summary_exp15.md", "PIC doğrulandı; 2B/3B yapılar 20-60 w_p^-1'de dağılıyor, B (W_ce = 2 w_p) koruyor",
         "kullanıldı: bulutta kalıcı uzay yükü yapısı yok; manyetik alan düzenler"),
        ("exp16", "results/exp16_equations/summary_exp16.md", "PySINDy bilinen yasaları buldu, gürültüye dayanıksız", "yöntem: denklem fizikten türetildi"),
        ("exp17", "results/exp17/summary_exp17.md", "j_on = 2,569 uA/um; V_M 0,80179, V_down 0,78552; w_s yasası; N ~ I/(0,75 j_on w_s)",
         "kullanıldı: hücre başına akım ve hücre sayısı"),
        ("exp18", "results/exp18/summary_exp18.md", "+30 mV tohum doğum yerini 10/10 seçiyor", "kullanıldı: hangi yayıcının ateşleneceği"),
        ("exp19", "results/exp19/summary_exp19.md", "ani V_DD adımı k = 1, 2, 4 hücre", "kullanıldı: itki basamakları"),
        ("exp20", "results/exp20/summary_exp20.md", "zorla akım düşüşüyle eşzamanlı bölünme 3 -> 7 (5,0 V)", "kullanıldı: hızlı itki artışı"),
        ("exp21", "results/exp21/summary_exp21.md", "300 K: 21 koşuda 0 bölünme", "olumsuz -> kısıt: hücre elektroniği kriyojenik (<= 4 K)"),
        ("exp22", "results/exp22/summary_exp22.md", "oksit tünellemesi 1 mm^2'de 1,8 V'ta fA-pA", "kullanıldı: kaçak kayıp ihmal edilebilir"),
        ("exp23", "results/exp23/summary_exp23.md", "hücre gücü 135,9/220 uW; aktarım 0,985-1,023; yük payı %74,6/83; gövde enerjisi 3-5e-14 J; "
         "taşıyıcı 0,82 eV; yanal 0,02-0,07, dikey 0,12 V/nm", "kullanıldı: enerji, geri dönüşüm, alan"),
        ("exp24", "results/exp24/ionic.json", "alpha = 0,507 eV/(V/nm), E_a = 1 eV, nu = 1e13/s, SET 1,466 V (4 K); 1-2 atom; yazma başına 2-4 O; "
         "pay >= 8 kT/e", "kullanıldı: programlanabilir yayma yasası, yük adresleme"),
        ("exp25", "results/exp25/summary_exp25.md", "hücre alanı 0,4-1,8 V/um; atom başına 2,8 meV (11,5 x Landauer)", "kullanıldı: denetim enerjisi, alan"),
        ("exp26", "results/exp26/summary_exp26.md", "doğal Si spin depolaması ölçülmüş döngüde çöküyor (%31)", "olumsuz -> kısıt: denetim yük durumuyla, spinle değil"),
        ("exp27", "results/exp27/summary_exp27.md", "foton sınırı F/P <= 1/c = 3,34 nN/W; vakum ve metrik ayakları kapalı", "karşılaştırma tabanı"),
    ]
    save("data_registry", [dict(exp=a, file=b, data=c, role=d) for a, b, c, d in R])


def rocket():
    """constant power P, efficiency eta: dM = -mdot dt, M dv = mdot v_e dt, P = mdot v_e^2 / (2 eta); bookkeeping in the lab"""
    out = {}
    M0, Mp, P, eta, v = 10.0, 2.0, 100.0, 0.6, 60e3
    mdot = 2 * eta * P / v**2
    n = 200000
    dt = Mp / mdot / n
    M, u, KE_ex, E_in = M0, 0.0, 0.0, 0.0
    for _ in range(n):
        dm = mdot * dt
        u_ex = u - v                                   # exhaust velocity in the lab (Galilean)
        u_new = (M * u - dm * u_ex) / (M - dm)          # momentum conservation
        KE_ex += 0.5 * dm * u_ex**2
        M, u = M - dm, u_new
        E_in += P * dt
    KE_tot = 0.5 * M * u**2 + KE_ex
    out["ode"] = dict(dv_numeric=u, dv_closed=v * np.log(M0 / M), rel_err=abs(u - v * np.log(M0 / M)) / (v * np.log(M0 / M)),
                      energy_balance=(KE_tot - eta * E_in) / (eta * E_in), burn_time_s=Mp / mdot)
    out["thrust_per_power"] = {f"v_e={x/1e3:g} km/s, eta={e:g}": 2 * e / x for x in (10e3, 30e3, 60e3, 91e3) for e in (0.26, 0.6, 1.0)}
    out["photon_bound_N_per_W"] = 1 / C
    out["atom_over_photon(eta=1)"] = {f"{x/1e3:g} km/s": 2 * C / x for x in (10e3, 30e3, 60e3, 91e3)}
    # measured thrusters: eta = F v_e / (2 P)
    meas = []
    for name, F, Pw, isp in (("NSTAR (DS1)", 92e-3, 2300.0, 3100.0), ("NSTAR (DS1)", 92e-3, 2300.0, 3300.0),
                             ("Enpulsion IFM Nano", 0.35e-3, 40.0, 2000.0), ("Enpulsion IFM Nano", 0.35e-3, 40.0, 6000.0),
                             ("MIT electrospray 2007", 13e-6, 0.275, 3000.0)):
        vex = isp * G0
        meas.append(dict(name=name, F_N=F, P_W=Pw, Isp_s=isp, v_e_m_s=vex, F_over_P_uN_W=F / Pw * 1e6, eta=F * vex / (2 * Pw)))
    out["measured"] = meas
    save("rocket", out)


def cells():
    out = {}
    I_c = [p / V_D for p in P_CELL]
    out["I_cell_A"] = I_c
    out["width_from_j_on_um"] = [i / J_ON * 1e6 for i in I_c]
    rows = []
    for sp in ("In", "Xe", "O"):
        for V in (1e3, 2e3, 5e3):
            for i in I_c:
                mq = PROP[sp]["m"] * AMU / E
                F = 0.9 * i * mq * ve(sp, V)                          # beam factor 0.9 (divergence, mass use)
                Pw = i * (V + V_M + PROP[sp]["I"])                     # acceleration + cell + ionisation (per unit charge)
                rows.append(dict(species=sp, V_acc=V, I_cell_A=i, v_e_km_s=ve(sp, V) / 1e3, F_cell_uN=F * 1e6, P_cell_W=Pw,
                                 F_over_P_uN_W=F / Pw * 1e6, eta_V=V / (V + V_M + PROP[sp]["I"])))
    out["per_cell"] = rows
    i = I_c[1]
    F1 = 0.9 * i * PROP["In"]["m"] * AMU / E * ve("In", 5e3)
    out["throttle_steps_uN"] = {f"N={k}": k * F1 * 1e6 for k in (1, 2, 3, 4, 7)}
    # cell death: total delivered power changes by the measured transfer (Exp23: 392.9 -> 390.2 uW in the worst case)
    out["fault_tolerance"] = dict(total_power_change_worst=(390.2 - 392.9) / 392.9, transfer_eta=ETA_REC,
                                  thrust_change_bound=max(abs(1 - e) for e in ETA_REC) / 3)
    out["architecture"] = ("Exp23: 74.6-83 % of the source power heats the load resistor; with the ion beam as the load this share "
                           "becomes beam power, the cell itself keeps only V_M/V_acc")
    save("cells", out)


def emission():
    out = {}
    rows = []
    for name, Q0, T in (("Exp24 O/HfO2 law, 4 K", Q_OX, 4.0), ("Exp24 O/HfO2 law, 430 K", Q_OX, 430.0),
                        ("In field evaporation (image hump Q0 = L + I - phi), 430 K", PROP["In"]["L"] + PROP["In"]["I"] - PROP["In"]["phi"], 430.0),
                        ("Ga field evaporation, 300 K", PROP["Ga"]["L"] + PROP["Ga"]["I"] - PROP["Ga"]["phi"], 300.0)):
        kT = KB * T / E
        Fth = Q0 / ALPHA                                         # linear-law threshold (barrier -> 0)
        ratio_10 = np.exp(0.1 * Q0 / kT)                         # 10 % field change below threshold
        rows.append(dict(case=name, Q0_eV=Q0, T_K=T, F_threshold_V_nm=Fth, on_off_for_10pct_field=ratio_10 if ratio_10 < 1e300 else 1e300,
                         log10_on_off=0.1 * Q0 / kT / np.log(10)))
    out["switching"] = rows
    out["addressing"] = dict(birth_seed_mV=30, birth_success="10/10 (Exp18)", charge_margin="8 kT/e (Exp24)",
                             control_energy_meV=2.8, landauer_4K_meV=0.24)
    save("emission", out)


def tunnel():
    """WKB through the triangular barrier U(x) = I - e F x (Exp1/4 method with the free-electron mass in vacuum)"""
    out = {}
    out["b_constant"] = B_WKB
    # numeric check of the WKB integral for one case
    I, F = 5.786, 6.0
    x = np.linspace(0, I / F, 200001)
    kappa = np.sqrt(2 * ME * np.maximum(I - F * x, 0) * E) / HBAR * 1e-9          # 1/nm
    out["numeric_vs_closed"] = dict(numeric=2 * np.trapezoid(kappa, x), closed=B_WKB * I**1.5 / F)
    rows = []
    for sp in ("Cs", "In", "Ga", "Xe", "O"):
        I = PROP[sp]["I"]
        for nu in (1e15, 1e16):
            F = B_WKB * I**1.5 / np.log(nu / 1e9)
            rows.append(dict(species=sp, I_eV=I, attempt_hz=nu, F_for_1e9_per_s_V_nm=F,
                             beta_needed_vs_cells=[F / f for f in (min(F_CELL), max(F_CELL))],
                             dlnk_dlnF=np.log(nu / 1e9)))
    out["field_ionisation"] = rows
    out["cell_fields_V_nm"] = F_CELL
    save("tunnel", out)


def turnstile():
    out = {}
    out["island_max_nm_for_EC_10kT_at_4K"] = {"eps3.9": EC_KNM / 40.0, "eps11.7": EC_KNM / 3 / 40.0}
    out["exp2_chain"] = dict(sizes_nm=[12, 6, 3], EC_K=[119.8, 231.5, 404.0], quantisation_4K=[0.926, 0.982, 1.000])
    rows = []
    for sp in ("In", "Xe", "O"):
        for f in (1e8, 1e9, 1e10):
            F = PROP[sp]["m"] * AMU * ve(sp, 5e3) * f
            rows.append(dict(species=sp, f_Hz=f, current_A=E * f, thrust_N=F, turnstiles_for_1mN=1e-3 / F))
    out["atom_by_atom"] = rows
    out["exp4_loading"] = "t = 24.8 meV, rise ~59 ps at gamma = 1 (L 1 nm, d 0.5 nm): GHz cycles possible"
    save("turnstile", out)


def charge():
    out = {}
    rows = []
    for name, F in (("Exp25 cell field 0.4 V/um", 0.4e6), ("Exp25 cell field 1.8 V/um", 1.8e6), ("Exp5 surface 0.38 MV/cm", 0.38e8),
                    ("gridded ion thruster ~1.5 kV/mm", 1.5e6)):
        rows.append(dict(case=name, E_V_m=F, thrust_density_N_m2=8 / 9 * EPS0 * F**2))
    out["child_langmuir"] = rows
    out["nstar_check_N_m2"] = 92e-3 / (np.pi * 0.15**2)
    out["plume"] = "Exp15: 2D/3D phase-space structures break up in 20-60 w_p^-1; a magnetic field (W_ce = 2 w_p) keeps them"
    save("charge", out)


def mission():
    out = {}
    rows = []
    for vex in (30e3, 60e3, 91.7e3):
        for eta in (0.26, 0.6):
            for P in (100.0,):
                M0, Mp = 10.0, 2.0
                t = 0.5 * Mp * vex**2 / (eta * P)
                rows.append(dict(v_e_km_s=vex / 1e3, eta=eta, P_W=P, dv_km_s=vex * np.log(M0 / (M0 - Mp)) / 1e3, burn_years=t / 3.156e7,
                                 thrust_mN=2 * eta * P / vex * 1e3))
    out["craft_10kg_100W_2kg"] = rows
    out["photon_li_ion_dv_m_s"] = 2.7e-3
    best = max(r["dv_km_s"] for r in rows if r["burn_years"] <= 5)
    out["ratio_atom_over_photon"] = best * 1e3 / 2.7e-3
    out["stuhlinger"] = "with power-source mass P/a_s, the optimum v_e ~ sqrt(2 eta a_s t) (characteristic velocity)"
    save("mission", out)


def fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    D = {k: json.loads((OUT / f"{k}.json").read_text(encoding="utf-8")) for k in ("rocket", "cells", "emission", "tunnel", "mission")}
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    a = ax[0, 0]
    v = np.logspace(3.5, 5.3, 100)
    for eta in (0.26, 0.6, 1.0):
        a.loglog(v / 1e3, 2 * eta / v * 1e6, label=f"atom: 2η/v_e, η = {eta}")
    a.axhline(1e6 / C, color="k", ls="--", label="foton sınırı 1/c (Exp27)")
    for m in D["rocket"]["measured"]:
        a.plot(m["v_e_m_s"] / 1e3, m["F_over_P_uN_W"], "o", label=f"{m['name']} (ölçülen)")
    a.set(xlabel="çıkış hızı v_e (km/s)", ylabel="itki/güç (µN/W)", title="(a) atom pozitif enerjiyle: F/P = 2η/v_e")
    a.legend(fontsize=6)
    a = ax[0, 1]
    for sp in ("In", "Xe", "O"):
        r = [x for x in D["cells"]["per_cell"] if x["species"] == sp and abs(x["I_cell_A"] - D["cells"]["I_cell_A"][1]) < 1e-9]
        a.plot([x["V_acc"] / 1e3 for x in r], [x["F_cell_uN"] for x in r], "o-", label=f"{sp}⁺, I_hücre = {D['cells']['I_cell_A'][1]*1e3:.2f} mA")
    a.set(xlabel="hızlandırma gerilimi (kV)", ylabel="hücre başına itki (µN)", title="(b) bir hücre = bir kuantumlu akım kanalı")
    a.legend(fontsize=7)
    a = ax[0, 2]
    st = D["cells"]["throttle_steps_uN"]
    a.bar(list(st.keys()), list(st.values()))
    a.set(ylabel="toplam itki (µN), In⁺ 5 kV", title="(c) hücre sayısıyla itki basamakları (Exp19: 1,2,4; Exp20: 3→7)")
    a = ax[1, 0]
    F = np.linspace(0.5, 2.2, 300)
    for T, ls in ((4.0, "-"), (430.0, "--")):
        kT = KB * T / E
        a.semilogy(F, NU * np.exp(-np.maximum(0, Q_OX - ALPHA * F) / kT), ls, label=f"Exp24 yasası, {T:g} K")
    a.set(xlabel="yerel alan (V/nm)", ylabel="yayma/sıçrama hızı (1/s)", title="(d) programlanabilir atom: k = ν exp(−(Q₀−αF)/kT)", ylim=(1e-30, 1e14))
    a.legend(fontsize=7)
    a = ax[1, 1]
    Fv = np.logspace(-0.5, 1.6, 300)
    for sp in ("Cs", "In", "Xe"):
        I = PROP[sp]["I"]
        a.loglog(Fv, 1e15 * np.exp(-B_WKB * I**1.5 / Fv), label=f"{sp} (I = {I} eV)")
    a.axvspan(min(F_CELL), max(F_CELL), color="grey", alpha=0.3, label="hücrelerin kendi alanları")
    a.set(xlabel="alan (V/nm)", ylabel="alan iyonlaşma hızı (1/s)", title="(e) vakuma tünelleme (WKB, Exp1/4 yöntemi)", ylim=(1e-10, 1e16))
    a.legend(fontsize=7)
    a = ax[1, 2]
    for r in D["mission"]["craft_10kg_100W_2kg"]:
        a.plot(r["burn_years"], r["dv_km_s"], "o")
        a.annotate(f"{r['v_e_km_s']:.0f} km/s, η={r['eta']}", (r["burn_years"], r["dv_km_s"]), fontsize=7)
    a.set(xlabel="yanma süresi (yıl), 100 W", ylabel="Δv (km/s)", title="(f) 10 kg araç, 2 kg In: Δv = v_e ln(M₀/M₁)")
    fig.tight_layout()
    fig.savefig(OUT / "fig_exp28.png", dpi=125)
    print("saved", OUT / "fig_exp28.png")


if __name__ == "__main__":
    {"registry": registry, "rocket": rocket, "cells": cells, "emission": emission, "tunnel": tunnel, "turnstile": turnstile,
     "charge": charge, "mission": mission, "fig": fig}[sys.argv[1]]()
