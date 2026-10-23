"""Exp25 figure: energy routing, Stark tuning, field-coupled atomic NAND, universal central atom (hyperfine register)."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ap25_register as reg  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "exp25"


def main():
    en = json.loads((OUT / "energy.json").read_text(encoding="utf-8"))
    iso = json.loads((OUT / "register_isotopes.json").read_text(encoding="utf-8"))
    lc = OUT / "logic_calibrated.json"
    lg = json.loads(lc.read_text(encoding="utf-8")) if lc.exists() else None
    f, axs = plt.subplots(2, 3, figsize=(18, 10))
    ax = axs[0, 0]
    keys = list(en["routing"])
    sh = np.array([en["routing"][k]["share"] for k in keys])
    x = np.arange(len(keys))
    ax.bar(x - 0.2, sh[:, 0], 0.4, label="uzak hücre (0)")
    ax.bar(x + 0.2, sh[:, 1], 0.4, label="yakın hücre (1)")
    ax.set_xticks(x)
    ax.set_xticklabels(["yönlendirmesiz", "0'a +30 mV", "1'e +30 mV", "0'a +15 mV", "1'e +15 mV"], fontsize=8)
    ax.set_ylabel("ölen hücrenin gücünden pay")
    ax.set_title("(a) enerji yönlendirme: pay yakınlıkla belirleniyor (%73/%27),\nkapı ön gerilimi ±%4 değiştiriyor; verim 0,97–0,98", fontsize=9)
    ax.legend(fontsize=8)
    ax = axs[0, 1]
    E = np.linspace(0.05, 3.5, 200)
    ax.semilogy(E, 117.53 * 2e-2 * E**2 * 1e3, label="³¹P aşırı ince (kuram −2e-2 µm²/V²)")
    ax.semilogy(E, 117.53 * 3.7e-3 * E**2 * 1e3, label="Sb benzeri (ölçülen −3,7e-3)")
    ax.semilogy(E, 27.97e6 * 1e-5 * E**2, label="ESR, g-faktörü, 1 T")
    for wf in (0.5, 1.0, 2.2):
        ax.axvline(0.9 / wf, c="0.6", ls=":", lw=1)
        ax.text(0.9 / wf, 2e0, f" cephe {wf} µm", fontsize=7, rotation=90)
    ax.axvline(3.7, c="C3", lw=1)
    ax.text(3.75, 2e0, "iyonlaşma", fontsize=7, color="C3", rotation=90)
    ax.set_xlabel("hücrenin alanı (V/µm)")
    ax.set_ylabel("enerji düzeyi kayması (kHz)")
    ax.set_title("(b) enerji düzeyi manipülasyonu: hücre alanı atomun\nrezonansını kHz–MHz kaydırıyor (frekansla seçim)", fontsize=9)
    ax.legend(fontsize=7)
    ax = axs[0, 2]
    if lg:
        best = None
        for run in lg["runs"]:
            for t in run["triples"]:
                if t.get("works"):
                    best = (run["seed"], t)
                    break
            if best:
                break
        if best:
            rows = best[1]["rows"]
            ax.bar(range(4), rows, color=["C2", "C2", "C2", "C0"])
            ax.set_xticks(range(4))
            ax.set_xticklabels(["(0,0)→1", "(0,1)→1", "(1,0)→1", "(1,1)→0"])
            ax.set_ylim(0, 1.05)
            ax.axhline(0.9, c="k", ls=":", lw=1)
        s = lg["summary"]
        ax.set_title(f"(c) 3 doğal atomdan NAND (alanla bağlı, 4 K, kalibre eşik):\n{s['working']}/{s['triples']} üçlü, "
                     f"{s['islands_with_NAND']}/{s['islands']} adada; her atomun alanı 30 nm içindekilerin hepsinde > kT", fontsize=9)
    ax.set_ylabel("doğru çıkış olasılığı")
    ax = axs[1, 0]
    rng = np.random.default_rng(3)
    nu = reg.GN * reg.B_FIELD
    for ab, col, lab in ((0.04685, "C0", "doğal Si (%4,7 ²⁹Si)"), (0.001, "C3", "%0,1 ²⁹Si")):
        occ = rng.random(len(reg.SITES)) < ab
        A = reg.A_SITES[occ]
        reg.DMIN = 10e3
        ok, fp, fm = reg.addressable(A, nu)
        ax.vlines((fp - nu) / 1e3, 0, 1 if ab > 0.01 else 0.6, color=col, lw=0.5, alpha=0.5, label=lab)
        ax.vlines((fp[ok] - nu) / 1e3, 0, 1.2, color=col, lw=2)
    ax.set_xlim(-50, 1600)
    ax.set_xlabel("çekirdek rezonansı − ν_n (kHz), elektron ↑")
    ax.set_title("(d) merkez atomun komşuları: her ²⁹Si'nin kendi frekansı;\nkalın = tek tek adreslenebilir (10 kHz seçicilik)", fontsize=9)
    ax.legend(fontsize=8)
    ax.set_yticks([])
    ax = axs[1, 1]
    for dmin, mk in ((100, "o"), (30, "s"), (10, "^")):
        ab = [4.685, 1.0, 0.1, 0.05]
        y = [iso[f"29Si={a:g}%, Dmin={dmin}kHz"]["mean_addressable"] for a in ab]
        z = [iso[f"29Si={a:g}%, Dmin={dmin}kHz"]["frac_strong_within_3nm_addressable"] or 0 for a in ab]
        ax.semilogx(ab, y, "-" + mk, label=f"adreslenebilir sayı, {dmin} kHz")
        ax.semilogx(ab, np.array(z) * 6, "--" + mk, alpha=0.6, label=f"güçlü komşulardan pay ×6, {dmin} kHz")
    ax.set_xlabel("²⁹Si oranı (%)")
    ax.set_ylabel("donör başına adreslenebilir çekirdek")
    ax.set_title("(e) 'tüm komşular' ödünleşimi: seyreltme payı artırır,\nsayıyı azaltır", fontsize=9)
    ax.legend(fontsize=6)
    ax = axs[1, 2]
    keys = [k for k in iso if iso[k]["ghz_mean"]]
    labels = [k.replace("29Si=", "").replace(", Dmin=", "\n") for k in keys]
    for t2, col in (("no_dephasing", "C2"), ("T2e=10ms", "C1"), ("T2e=0.3ms", "C3")):
        ax.plot(range(len(keys)), [iso[k]["ghz_mean"][t2] for k in keys], "o-", c=col, label=t2.replace("no_dephasing", "eşevrelilik kaybı yok"))
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.axhline(0.9, c="k", ls=":", lw=1)
    ax.set_ylabel("GHZ doğruluğu (elektron + 4 çekirdek)")
    ax.set_title("(f) evrensel kontrol: tek RF bobiniyle 5 parçacıklı dolanıklık;\nseyircilerle 0,98–0,99; sınır elektron T2e", fontsize=9)
    ax.legend(fontsize=8)
    f.suptitle("Exp25: atom, alan ve enerji manipülasyonu — evrensel programlanabilir atom", fontsize=12)
    f.tight_layout(rect=(0, 0, 1, 0.95))
    f.savefig(OUT / "fig_exp25.png", dpi=90)


if __name__ == "__main__":
    main()
