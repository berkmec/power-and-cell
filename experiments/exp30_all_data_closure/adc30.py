"""Exp30: every Exp0-Exp27 result closed against an independent MEASUREMENT (or marked 'no measurement' / 'method'),
plus measured field-controlled atom-emission propulsion (indium FEEP) against the Exp28 equation.
Criteria were fixed in results/exp30/predictions_exp30.json (14:56) before the anchor sources were read.
Level of reading for every anchor is recorded ('tam metin' = full text read in this session, 'özet' = abstract/search summary).
usage: python adc30.py [fig]  -> results/exp30/closure_all.json, fig_exp30.png"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[2] / "results" / "exp30"
E, EPS0, ME, HBAR, KB, C = 1.602176634e-19, 8.8541878128e-12, 9.1093837015e-31, 1.054571817e-34, 1.380649e-23, 2.99792458e8


def calc():
    r = {}
    # Exp0-2: Coulomb island. Exp1 FEM: C = c_sigma * eps_r * eps0 * a, a = cube edge (cap_cell.py box_bounds), eps_r = 3.9.
    d = np.array([1.58, 2.0])                                 # Shin 2010: 1.58 nm (from C_sigma), ~2 nm (TEM)
    a = np.concatenate([d * (np.pi / 6) ** (1 / 3), d])        # volume-equivalent cube edge ... edge = diameter
    Cm = {rr: c * 3.9 * EPS0 * a * 1e-9 for rr, c in ((0.5, 22.989), (0.2, 16.267))}
    meas_C = 0.42e-18
    r["Exp1_capacitance"] = dict(a_nm=[float(a.min()), float(a.max())], C_model_aF={str(k): [float(v.min() * 1e18), float(v.max() * 1e18)] for k, v in Cm.items()},
                                 C_meas_aF=0.42, ratio_range=[float(min(v.min() for v in Cm.values()) / meas_C), float(max(v.max() for v in Cm.values()) / meas_C)],
                                 Ec_model_eV=[float(E / Cm[0.5].max()), float(E / Cm[0.5].min())], Ec_meas_eV=0.377)
    r["Exp0-2_room_T_size"] = dict(model_a_max_nm_300K=0.4, measured_island_nm=2.0, factor=5.0, verdict="✗",
                                   note="ölçüm: ~2 nm adada 300 K'de dört Coulomb elması; 5 nm adada 200 K altında belirgin blokaj")
    r["Exp2_4K"] = dict(pump_T_e_K=[7, 9], E_add_meV=17, E_over_kT=[17e-3 * E / (KB * 9), 17e-3 * E / (KB * 7)], verdict="✓")
    # Exp4: kappa from measured barrier (3.1 eV) and FN tunnelling mass (~0.4)
    k = lambda m, phi: np.sqrt(2 * m * ME * phi * E) / HBAR * 1e-9
    r["Exp4_kappa"] = dict(kappa_meas_nm=float(k(0.4, 3.1)), kappa_band_nm=[float(k(0.35, 3.1)), float(k(0.45, 3.1))], model_nm=[4.4, 5.7],
                           dev_mid=float(5.05 / k(0.4, 3.1) - 1), verdict="✓" if abs(5.05 / k(0.4, 3.1) - 1) <= 0.2 else "✗")
    r["Exp5_VT"] = dict(model_V=0.513, pdk_V=0.519, dev=0.513 / 0.519 - 1, verdict="✓")
    # Exp6: Pelgrom A_VT ~ 1 mV um per nm oxide (4.15 nm) ... typical 180 nm 5 mV um; device 126 nm x 420 nm
    area = np.sqrt(0.126 * 0.42)
    s = np.array([4.15, 5.0]) / area
    r["Exp6_mismatch"] = dict(sigma_meas_mV=s.tolist(), sigma_MC_mV=[3.9, 5.1], ratio=[float(s.min() / 5.1), float(s.max() / 3.9)], verdict="✗")
    seg = np.array([4.15, 5.0]) / np.sqrt(10 * 0.15)
    r["Exp18_seed_vs_mismatch"] = dict(sigma_seg_mV=seg.tolist(), seed_mV=30, seed_in_sigma=[float(30 / seg.max()), float(30 / seg.min())],
                                       assumption="10 µm segment, L = 0,15 µm (en kötü durum)")
    r["Exp8_fission"] = dict(model_exponent=-1.08, rayleigh=-1.0, dev=0.08, verdict="✓")
    r["Exp22_subband"] = dict(model_single_subband_to=9.5e12, measured_second_onset=6e12, verdict="✗ (önceden kaydedilmemiş ek karşılaştırma)")
    r["Exp25_GHZ"] = dict(model=0.987, measured=[0.925, 0.962], gap_points=[2.5, 6.2], verdict="✗")
    r["Exp26_coherence"] = dict(model=0.86, measured=float(np.exp(-17.5 / 25.0)), verdict="✗")
    r["Exp27_IKAROS"] = dict(F_meas_mN=1.12, F_max_2SA_over_c_1AU_mN=2 * 1361 * 196 / C * 1e3, ratio=1.12e-3 / (2 * 1361 * 196 / C), verdict="✓")
    r["FEEP"] = dict(dev_pct={"V_em sabit, 2–4 mA": 1.6, "3,5 mA, 3–10 kV": 2.8, "V_ex sabit": 4.3}, criterion_pct=5.0, verdict="✓",
                     needles={"37 µA, 6 kV": "28 iğneden 1'i", "3,5 mA": "28/28"}, per_needle_uA=3.5e3 / 28, cell_Ic_uA=[166, 268])
    r["In_field"] = dict(image_hump_V_nm=4.19**2 / 1.44, measured="LMIS uç alanı > 10 V/nm (ölçülmüş tek In değeri bulunamadı)", verdict="kısmen")
    return r


ROWS = [
    # exp, project result, measured anchor, source, reading, verdict, role in the final equation
    ("Exp0", "Coulomb + hapsolma, özyinelemeli boyutlar", "~2 nm Si adada 300 K'de Coulomb + kuantum etkileri", "Shin 2010, arXiv 1001.3724", "tam metin",
     "hapsolma ✓ (nitel) / Coulomb ✗", "turnike adası boyutu"),
    ("Exp1", "c_Σ = 22,989 ε a (FEM)", "C_Σ = 0,42 aF, e²/C_Σ = 0,377 eV (1,58–2 nm)", "Shin 2010", "tam metin", "✗ (model 1,7–3,8× fazla)", "turnike E_C yasası"),
    ("Exp2", "4 K'de yük kuantumu; 77/300 K'de yok", "Si pompası T_e 7–9 K'de ppm; 5 nm adada < 200 K blokaj", "Zhao 2017; Shin 2010", "tam metin",
     "4 K ✓ / 77 K ✗", "turnike sıcaklığı"),
    ("Exp3", "300 K çift nokta ölü", "Si yük kübiti T2* 0,13–2,1 ns, seyreltme soğutucusunda", "Shi 2013, arXiv 1208.0519", "özet", "✓", "kısıt: kriyojenik"),
    ("Exp4", "WKB κ 4,4–5,7 nm⁻¹", "Si/SiO₂ bariyeri 3,1 eV, FN tünel kütlesi ≈ 0,4 → 5,70 nm⁻¹", "Lenzlinger–Snow 1969; arXiv 1512.00255", "özet", "✓ (orta −%11)",
     "tünelleme, GHz turnike"),
    ("Exp5", "V_g(n_s 1e12) 0,513 V", "PDK vth0 0,519 V (ölçüme uydurulmuş kart)", "SkyWater PDK", "birincil dosya", "✓ (%1,2)", "SKY130 alan ölçeği"),
    ("Exp6", "MC σ_VT 3,9–5,1 mV", "Pelgrom A_VT ≈ 4–5 mV·µm → 18–22 mV", "Pelgrom (Berkeley EE240B)", "özet", "✗ (3,5–5,6× küçük)", "değişkenlik"),
    ("Exp7", "≥ 100 aygıt test planı", "—", "—", "—", "yöntem", "ölçüm planı"),
    ("Exp8", "fisyon eşiği üssü −1,08", "yüklü damlalar Rayleigh sınırında bölünüyor", "elektrodinamik terazi ölçümleri", "özet", "✓ (%8)", "bölünme ilkesi"),
    ("Exp9", "araç çapraz denetimi", "—", "—", "—", "yöntem", "doğrulama"),
    ("Exp10", "2B'de K < 0 var, bölünme yok", "düşük yoğunlukta negatif sıkıştırılabilirlik ölçülmüş", "Eisenstein 1992, PRL 68, 674", "özet",
     "kısmen (eşik sayısı yok)", "kısıt: filaman hücre"),
    ("Exp11", "otokataliz gerekli", "kendini çoğaltan noktalar deneyde", "Lee 1994, Nature 369, 215", "özet", "✓", "ilke: otokataliz"),
    ("Exp12", "4 K yüzer gövde kilitlenmesi, kink + histerezis", "180 nm NMOS'ta kink + histerezis, 30 K'de var, 40 K'de yok", "Huang 2020, arXiv 2003.00207",
     "tam metin", "✓ (V_M = 0,80 V ölçülmedi)", "hücre gerilim kenedi"),
    ("Exp13", "balast λ genişliği seçer", "—", "—", "—", "ölçüm yok", "hücre genişliği"),
    ("Exp14", "bölünme boyut penceresi", "—", "—", "—", "ölçüm yok", "kısıt: 1B şerit"),
    ("Exp15", "1B delikler birleşir; B korur", "manyetize Q-makinesinde elektron delikleri ve birleşme", "Saeki 1979, PRL 42, 501", "özet (atıflar)", "✓",
     "iyon bulutu"),
    ("Exp16", "PySINDy", "—", "—", "—", "yöntem", "yöntem"),
    ("Exp17", "j_on, w_s yasası, N ≈ I/(0,75 j_on w_s)", "Si'de yok; benzer: FEEP'te yayan iğne sayısı akımla artıyor", "Mühlich 2022", "tam metin",
     "ölçüm yok (Si)", "I_c, N(t)"),
    ("Exp18", "+30 mV tohum 10/10", "tohum, ölçek uyumsuzluğunun 7–9σ'sı", "Pelgrom ölçeği", "özet", "tutarlı (dolaylı)", "doğum yeri"),
    ("Exp19", "V_DD adımı k = 1, 2, 4", "—", "—", "—", "ölçüm yok", "itki basamakları"),
    ("Exp20", "zorla düşüş 3 → 7", "—", "—", "—", "ölçüm yok", "hızlı artış"),
    ("Exp21", "300 K'de bölünme yok", "kink 40 K'de kayboluyor; n-GaAs filamanları 4,2 K'de", "Huang 2020; Brandl 1989", "tam metin", "✓", "kısıt: ≲ 30 K"),
    ("Exp22", "4 K'de 9,5e12'ye kadar tek alt bant; kaçak fA–pA", "Si(100)'de ikinci alt bant ~6e12 üstünde; kaçak sayısı bulunamadı", "Ray–Martin (IBM)",
     "özet", "✗ alt bant (ek) / kaçak ölçüm yok", "kayıp"),
    ("Exp23", "aktarım 0,985–1,023", "—", "—", "—", "ölçüm yok", "η_aktarım"),
    ("Exp24", "tek atomluk boğaz ~1,2 G₀; tek ³¹P RF hedef 1,000000", "HfO₂'de G₀ basamakları; ³¹P kapı > %98, okuma > %99,8",
     "Long 2013, APL 102, 183505; Pla 2013, Nature 496, 334", "özet", "kısmen (var ✓, sayı iyimser)", "programlanabilir atom"),
    ("Exp25", "GHZ 0,987", "Si GHZ %92,5 / %96,2", "Mądzik 2022; Thorvaldson 2025", "tam metin", "✗ (iyimser)", "kuantum denetim"),
    ("Exp26", "17,5 s'de 0,86", "MREV-16 T2 = 25 s → 0,50", "Ladd 2005", "tam metin", "✗ (iyimser; çöküş sonucu daha da güçlü)", "kısıt: yükle denetim"),
    ("Exp27", "F/P ≤ 1/c (yansımada 2/c); Casimir; DCE", "IKAROS 1,12 mN ≤ 1,78 mN; Casimir %5; DCE fotonları", "JAXA; Lamoreaux 1997; Wilson 2011",
     "özet / tam metin", "✓", "karşılaştırma tabanı"),
]


def fig(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"✓": "#2ca02c", "✗": "#d62728", "kıs": "#ff7f0e", "tut": "#9acd32", "ölç": "#bdbdbd", "yön": "#e0e0e0", "hap": "#ff7f0e", "4 K": "#ff7f0e"}
    fig, ax = plt.subplots(1, 2, figsize=(17, 8.2), gridspec_kw=dict(width_ratios=[1.25, 1]))
    a = ax[0]
    for i, row in enumerate(ROWS):
        v = row[5]
        key = next((k for k in col if v.startswith(k)), "ölç")
        a.barh(i, 1, color=col[key])
        a.text(0.02, i, f"{row[0]}: {v}", va="center", fontsize=7.5)
    a.set(yticks=[], xticks=[], ylim=(len(ROWS) - 0.5, -0.5), title="(a) Exp0–27: her sonuç, bağımsız ölçüme karşı")
    a = ax[1]
    R = res["calc"]
    labels = ["Exp1 C_model/C_ölç", "Exp6 σ_ölç/σ_MC", "Exp4 κ_model/κ_ölç", "Exp5 V_T", "Exp8 üs", "Exp27 F/F_max", "Exp25 GHZ", "Exp26 eşevrelilik",
              "FEEP itki (Exp28 denk.)", "NEXT itki (Exp29)"]
    lo = [R["Exp1_capacitance"]["ratio_range"][0], R["Exp6_mismatch"]["ratio"][0], 4.4 / R["Exp4_kappa"]["kappa_meas_nm"], 0.513 / 0.519, 1.08,
          R["Exp27_IKAROS"]["ratio"], 0.987 / 0.962, 0.86 / R["Exp26_coherence"]["measured"], 1.016, 0.9855]
    hi = [R["Exp1_capacitance"]["ratio_range"][1], R["Exp6_mismatch"]["ratio"][1], 5.7 / R["Exp4_kappa"]["kappa_meas_nm"], 0.513 / 0.519, 1.08,
          R["Exp27_IKAROS"]["ratio"], 0.987 / 0.925, 0.86 / R["Exp26_coherence"]["measured"], 1.043, 1.0145]
    y = np.arange(len(labels))
    a.hlines(y, lo, hi, lw=6, color="C0")
    a.plot(lo, y, "|", color="k"); a.plot(hi, y, "|", color="k")
    a.axvline(1, color="k", lw=0.8)
    a.set(yticks=y, yticklabels=labels, xscale="log", xlabel="model / ölçüm (1 = tam uyum; itkiler: ölçülen/öngörülen)",
          title="(b) nicel karşılaştırmalar", ylim=(len(labels) - 0.5, -0.5))
    fig.tight_layout()
    fig.savefig(OUT / "fig_exp30.png", dpi=125)
    print("saved", OUT / "fig_exp30.png")


if __name__ == "__main__":
    res = dict(calc=calc(), rows=[dict(zip(("exp", "project", "anchor", "source", "reading", "verdict", "role"), r)) for r in ROWS])
    tally = {}
    for r in ROWS:
        v = r[5]
        k = "✓" if v.startswith("✓") else "✗" if v.startswith("✗") else "yöntem" if v == "yöntem" else "ölçüm yok" if v.startswith("ölçüm yok") else "karışık/kısmen"
        tally[k] = tally.get(k, 0) + 1
    res["tally"] = tally
    (OUT / "closure_all.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float), encoding="utf-8")
    for k, v in res["calc"].items():
        print(k, v)
    print("tally", tally, "rows", len(ROWS))
    if "fig" in sys.argv:
        fig(res)
