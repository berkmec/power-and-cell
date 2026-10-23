"""Exp26 figure (moved out of nr26.py to keep it under 500 lines): python nr26_fig.py -> results/exp26/fig_exp26.png"""
import json

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402

from nr26 import OUT, TC_L, drive                                    # noqa: E402


def main():
    L = json.loads((OUT / "law.json").read_text(encoding="utf-8"))
    Rn = json.loads((OUT / "run.json").read_text(encoding="utf-8"))
    Bu = json.loads((OUT / "bus.json").read_text(encoding="utf-8"))
    Me = json.loads((OUT / "measured.json").read_text(encoding="utf-8")) if (OUT / "measured.json").exists() else None
    fig, ax = plt.subplots(2, 4, figsize=(21, 9.5))
    a = ax[0, 0]
    a.loglog(L["W_Hz"], L["single_mean"], "o", ms=3, label="sayım (tek dal)")
    a.loglog(L["W_Hz"], L["double_mean"], "s", ms=3, label="sayım (iki dal)")
    a.loglog(L["W_Hz"], L["law_single"], "-", label="0,63 √(K/W)")
    a.loglog(L["W_Hz"], L["law_double"], "--", label="0,31 √(K/W)")
    a.axhline(L["n_total_mean"], color="k", lw=0.8, ls=":", label=f"6 nm içindeki tüm ²⁹Si ({L['n_total_mean']:.0f})")
    a.axvline(L["W_device_Hz"], color="r", lw=0.8, label=f"cihaz sınırı W = {L['W_device_Hz']:.0f} Hz")
    a.set(xlabel="çözünürlük W (Hz)", ylabel="tek tek seçilebilen çekirdek", title="(a) yalnız aşırı ince anahtar (D0): L1 yasası")
    a.legend(fontsize=7)
    a = ax[0, 1]
    d = np.array(Rn["dgrid_nm"])
    a.semilogx(d, Rn["coverage_mean"], label="yönü seçilebilen gradyan (ortalama)")
    a.semilogx(d, Rn["coverage_min"], ":", label="en kötü yapılandırma")
    a.semilogx(d, Rn["coverage_fixed_mean"], label="tek sabit genel yön (10 rastgele)")
    if "coverage_fixed_lattice_123_mean" in Rn:
        a.semilogx(d, Rn["coverage_fixed_lattice_123_mean"], color="gray", lw=0.8, label="[1,2,3] (örgü yönü: düzlemler çakışık)")
    a.semilogx(d, Rn["law_fixed_mean"], "--", label="L2: ⟨exp(−d* S_k)⟩")
    a.axvline(Rn["d99_median_nm"], color="r", lw=0.8, label=f"d*₉₉ = {Rn['d99_median_nm']:.4f} nm")
    a.set(xlabel="gereken ayrım d* = W/(γG) (nm)", ylabel="kapsanan ²⁹Si oranı", title="(b) D+ kipinde alan kodlaması (6 nm)")
    a.legend(fontsize=7)
    a = ax[0, 2]
    br = Rn["b_vs_r_cfg0"]
    a.semilogy(br["r"], br["b"], ".", ms=2)
    a.axhline(Rn["d99_median_nm"], color="r", lw=0.8)
    a.set(xlabel="donörden uzaklık r (nm)", ylabel="en iyi yönde yalıtım b_k (nm)", title="(c) her çekirdeğin en iyi yalıtımı (yapılandırma 0)")
    a = ax[0, 3]
    for T2 in ("25", "100", "1000"):
        a.hist([x[f"T2={T2}"]["mean_degree"] for x in Bu["configs"]], bins=10, alpha=0.6, label=f"T2 = {T2} s (|J| ≥ 5/T2)")
    a.set(xlabel="ortalama komşu sayısı (dipol ağı)", ylabel="yapılandırma", title="(d) çekirdekler arası veri yolu (20/20 bağlı)")
    a.legend(fontsize=7)
    a = ax[1, 0]
    for key, rows in Rn["feasibility"].items():
        lab = "T2 ∝ t_c⁻² (tavan yok)" if key == "inf" else f"T2 tavanı {key} s"
        a.semilogx([x["tc_ms"] for x in rows], [x["mean_cov"] for x in rows], label=lab)
    a.axvline(TC_L, color="k", lw=0.8, ls=":", label="ölçülen nokta (Ladd 2005)")
    a.set(xlabel="MREV-16 döngüsü t_c (ms)", ylabel="en iyi kapsama (6 nm)", title="(e) katı ölçütle gerçekleşebilirlik (uzatma)")
    a.legend(fontsize=7)
    a = ax[1, 1]
    tcs = np.logspace(-2, np.log10(2), 60)
    a.loglog(tcs, [drive(t)["p_avg_W"] * 1e6 for t in tcs])
    a.axvline(TC_L, color="k", lw=0.8, ls=":")
    for key, fx in Rn["feasible_largest_tc"].items():
        if fx:
            a.plot(fx["tc_ms"], drive(fx["tc_ms"])["p_avg_W"] * 1e6, "o", label=f"en büyük uygun t_c ({key})")
    a.set(xlabel="t_c (ms)", ylabel="ortalama sürüş gücü (µW, 100 nm hat)", title="(f) enerji bütçesi")
    a.legend(fontsize=7)
    if Me:
        a = ax[1, 2]
        for k in (1.5, 2.0, 2.5, 4.0):
            rows = [r for r in Me["rows"] if r["kappa"] == k]
            a.plot([r["inversion"] for r in rows], [r["R_all_nm"] for r in rows], "o-",
                   label=f"κ = {k} (çapraz {Me['crosstalk'][str(k)]:.0e})")
        a.set(xlabel="yazma doğruluğu (hedefte tersinme)", ylabel="bütün ²⁹Si'nin yazılabildiği yarıçap (nm)",
              title="(g) ölçülmüş noktada erişim (T2 = 25 s, ≤ 120 Hz)")
        a.invert_xaxis()
        a.legend(fontsize=7)
        a = ax[1, 3]
        for k in (1.5, 2.0, 2.5, 4.0):
            rows = [r for r in Me["rows"] if r["kappa"] == k]
            a.plot([r["inversion"] for r in rows], [r["cov_6nm"] for r in rows], "o-", label=f"κ = {k}")
        a.set(xlabel="yazma doğruluğu (hedefte tersinme)", ylabel="6 nm yazmacının kapsaması",
              title="(h) ölçülmüş noktada 6 nm (~2100 ²⁹Si)")
        a.invert_xaxis()
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_exp26.png", dpi=125)
    print("saved", OUT / "fig_exp26.png")


if __name__ == "__main__":
    main()
