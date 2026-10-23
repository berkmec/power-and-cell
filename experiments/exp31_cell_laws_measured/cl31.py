"""Exp31: the silicon cell laws (Exp12, 13, 17, 18, 19, 21, 23) against MEASURED silicon current filaments.
Data: W. Mamanee, 'Investigation of current filament behavior in smart power ESD devices', TU Wien thesis (2012, D. Pogany):
300 K silicon NPN ESD devices, TLP pulses, transient interferometric mapping (TIM). Full text read; values below are quoted.
Criteria fixed in results/exp31/predictions_exp31.json (15:14) before this script.
usage: python cl31.py [fig] -> results/exp31/cell_laws.json, fig_exp31.png"""
import json
import sys
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[2] / "results" / "exp31"
# Table 8 (thesis): filament width W_CF (µm) at stress current I_S (A) and filament number N
W = {(1.4, 2): (15.0, 16.6), (1.4, 3): (12.7, 13.4), (2.0, 2): (22.0, 23.0), (2.0, 3): (16.9, 18.6), (2.0, 4): (13.6, 15.2)}
N_RANGE = {0.9: (1, 2), 1.4: (2, 3), 2.0: (3, 4)}         # Fig. 73: N alternates pulse to pulse
BISTABLE = (0.45, 0.9)                                     # single CF below 0.41-0.45 A; N = 1 or 2 for 0.45-0.8 A (and 0.9 A)
L_MAX = 20.0                                               # max W_CF consistent with emitter-sink distance L ≈ 20 µm


def j_of(I, N):
    lo, hi = W[(I, N)]
    return I / (N * hi) * 1e3, I / (N * lo) * 1e3            # mA/µm (low, high)


def run():
    res = {"j_mA_per_um": {f"I={I} A, N={N}": j_of(I, N) for (I, N) in W}}
    mid = lambda I, N: float(np.mean(j_of(I, N)))
    cl1 = {f"N={N}": mid(2.0, N) / mid(1.4, N) - 1 for N in (2, 3)}
    res["CL1"] = dict(ratio_minus_1=cl1, passed=bool(all(abs(v) <= 0.10 for v in cl1.values())))
    cl2 = mid(2.0, 4) / mid(2.0, 2) - 1
    res["CL2"] = dict(ratio_minus_1=cl2, passed=bool(abs(cl2) <= 0.10))
    j_ref = float(np.mean([mid(1.4, 2), mid(2.0, 2)]))
    out3 = {}
    for ws in (L_MAX, 22.5):
        rows = {}
        for I, (a, b) in N_RANGE.items():
            Np = I / (0.75 * j_ref * 1e-3 * ws)
            rows[str(I)] = dict(N_pred=Np, N_meas=[a, b], ok=bool(a - 0.5 <= Np <= b + 0.5), pred_over_mean=Np / ((a + b) / 2))
        out3[f"w_s={ws}"] = rows
    res["CL3"] = dict(j_ref_mA_per_um=j_ref, table=out3, passed=bool(all(r["ok"] for r in out3[f"w_s={L_MAX}"].values())))
    Imax = j_ref * 1e-3 * np.array([L_MAX, 22.5])
    res["CL4"] = dict(single_filament_max_A=Imax.tolist(), bistable_window_A=BISTABLE, upper_edge_dev=float(Imax[0] / BISTABLE[1] - 1),
                      passed=bool(abs(Imax[0] / BISTABLE[1] - 1) <= 0.2),
                      lower_edge_note="0,45 A'deki alt kenar (iki filamanın en erken var olabildiği akım) Exp17'de öngörülmüyor: çok-kararlılık")
    res["CL5"] = dict(observed="I_S = 0,8 A'de yüksek akım yoğunluklu filaman ikiye bölünüyor; iki filaman gücü eşit paylaşıyor (TIM, kıvrımlı aygıt)", passed=True)
    res["CL6"] = dict(observed="tohumsuz: I_S < 0,45 A'de rastgele başlangıç, rastgelelik akımla azalıyor; tohumlu (NPN_Dope, P-gövde ortasında yüksek katkı): başlangıç tekrarlanabilir biçimde ortaya sabit",
                      passed=True)
    res["CL7"] = dict(observed="kısa darbede yüksek akımda daha çok filaman; sabit akımda N darbeden darbeye değişiyor (0,9 A: 1/2; 1,4 A: 2/3; 2 A: 3/4)", passed=True)
    res["CL8"] = dict(observed="I–V'de N başına ayrı dal; 2 filaman dalı 0,45 A'den itibaren, tutma gerilimi ~1 V düşük; kaynak bunu paralel balast R_A ile açıklıyor: ΔV = I R_A (1 − 1/N)",
                      passed=True, caveat="SKY130'un V_M = 0,80 V değeri ölçülmedi; tutma gerilimi bu aygıtta ~18–19 V")
    res["CL9"] = dict(observed="300 K'de bipolar filamanlar çoğalıyor ve bölünüyor (darbeli, ~100 ns–1 µs; yüksek akımda ısıl bozulma)", passed=True)
    res["CL10"] = dict(observed="zorla akım düşüşüyle eşzamanlı bölünme için ölçülmüş benzer bulunamadı", passed=None)
    res["CL11"] = dict(observed="boyut penceresi için ölçülmüş benzer bulunamadı", passed=None)
    res["redistribution"] = "yüksek akımda filamanlar arası akım yeniden dağılımı ölçülmüş (ısıl bozulma süresini kısaltıyor) — Exp23'ün nitel karşılığı"
    return res


UPDATES = {   # Exp30 rows that Exp31 closes with measured silicon / oxide / active-droplet data
    "Exp13": ("TIM'de filaman bölünmesi; balast R_A paralel paylaşım (ΔV = I R_A (1 − 1/N))", "Mamanee 2012 (tam metin)", "✓ (nitel)"),
    "Exp14": ("aktif damlalar kritik hacmin üstünde içi boş sıvı kabuğa dönüşüyor; aktif damla bölünmesi henüz ölçülmemiş",
              "Bergmann 2023, Nat. Commun., PMC10582082 (özet)", "kısmen (kabuk ✓, bölünme penceresi ölçülmedi)"),
    "Exp17": ("sabit N'de j sabit (%0,1 / %5,2); N = I/(0,75 j w_s) aralıkta (0,84–0,90×); tek filaman sınırı 0,889 A ≈ 0,9 A; j N ile %22 düşüyor",
              "Mamanee 2012 (tam metin)", "karışık (3 nicel ✓, 1 ✗)"),
    "Exp18": ("tohumsuz aygıtta rastgele doğum; yerel katkı tohumlu aygıtta doğum yeri ortaya sabit", "Mamanee 2012 (tam metin)", "✓"),
    "Exp19": ("kısa darbede yüksek akımda daha çok filaman; sabit akımda N darbeden darbeye değişiyor", "Mamanee 2012 (tam metin)", "✓ (nitel)"),
    "Exp20": ("V₂O₃'te akım düşürülünce filaman çok sayıda paralel filamana bölünüyor (kaynak: heterojenlik)", "Lange 2020, arXiv 2009.12536 (tam metin)",
              "kısmen (olgu ✓, mekanizma farklı)"),
    "Exp23": ("iki filaman gücü eşit paylaşıyor; yüksek akımda filamanlar arası akım yeniden dağılımı", "Mamanee 2012 (tam metin)", "✓ (nitel)"),
}


def final_tally():
    p = OUT.parents[0] / "exp30" / "closure_all.json"
    rows = json.loads(p.read_text(encoding="utf-8"))["rows"]
    for r in rows:
        if r["exp"] in UPDATES:
            r["anchor"], r["source"], r["verdict"] = UPDATES[r["exp"]]
            r["updated_by"] = "Exp31"
    tally = {}
    for r in rows:
        v = r["verdict"]
        k = "✓" if v.startswith("✓") else "✗" if v.startswith("✗") else "yöntem" if v == "yöntem" else "ölçüm yok" if v.startswith("ölçüm yok") else "karışık/kısmen"
        tally[k] = tally.get(k, 0) + 1
    out = dict(rows=rows, tally=tally)
    (OUT / "closure_final.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    return out


def fig(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    a = ax[0]
    for N, mk in ((2, "o"), (3, "s"), (4, "^")):
        pts = [(I, j_of(I, N)) for (I, n) in W if n == N]
        for I, (lo, hi) in pts:
            a.errorbar(I, (lo + hi) / 2, yerr=(hi - lo) / 2, fmt=mk, color=f"C{N}", capsize=3, label=f"N = {N}" if I == pts[0][0] else None)
    a.axhline(res["CL3"]["j_ref_mA_per_um"], color="k", ls=":", lw=0.8, label="j_ref (N = 2)")
    a.set(xlabel="toplam akım I_S (A)", ylabel="filaman başına j = I/(N·W) (mA/µm)", title="(a) sabit N'de j sabit (CL1), N ile düşüyor (CL2)", xlim=(1.2, 2.2))
    a.legend(fontsize=8)
    a = ax[1]
    I = np.linspace(0.2, 2.2, 200)
    for ws, ls in ((20.0, "-"), (22.5, "--")):
        a.plot(I, I / (0.75 * res["CL3"]["j_ref_mA_per_um"] * 1e-3 * ws), ls, color="C0", label=f"Exp17: N = I/(0,75 j w_s), w_s = {ws} µm")
    for Ii, (lo, hi) in N_RANGE.items():
        a.vlines(Ii, lo, hi, lw=6, color="C3", label="ölçülen N aralığı (TIM)" if Ii == 0.9 else None)
    a.axvspan(*BISTABLE, color="C1", alpha=0.15, label="1/2 iki-kararlılık penceresi")
    a.set(xlabel="toplam akım I_S (A)", ylabel="filaman (hücre) sayısı N", title="(b) sayı yasası (CL3) ve bölünme penceresi (CL4)")
    a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_exp31.png", dpi=125)
    print("saved", OUT / "fig_exp31.png")


if __name__ == "__main__":
    res = run()
    (OUT / "cell_laws.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float), encoding="utf-8")
    for k, v in res.items():
        print(k, v)
    ft = final_tally()
    print("final tally (Exp30 + Exp31)", ft["tally"], "rows", len(ft["rows"]))
    for r in ft["rows"]:
        print(f"  {r['exp']}: {r['verdict']}")
    if "fig" in sys.argv:
        fig(res)
