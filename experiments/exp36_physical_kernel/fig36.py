"""Exp36 figure: physical surface-kernel, background-aware, noisy single-shot atom control.
Rendered directly at print resolution (4500x3000 px) from preserved Exp36 v1-v5 JSON results.
No simulation or verdict is recomputed; pattern background bars are task-level means from stored records."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np

_here = Path(__file__).resolve()
ROOT = _here.parents[2] if (_here.parents[2] / "results" / "exp36").exists() else Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "exp36"
V = {
    "v1": json.loads((RES / "universal_control36_v1.json").read_text(encoding="utf-8")),
    "v2": json.loads((RES / "universal_control36_v2.json").read_text(encoding="utf-8")),
    "v3": json.loads((RES / "universal_control36_v3.json").read_text(encoding="utf-8")),
    "v4": json.loads((RES / "universal_control36_v4.json").read_text(encoding="utf-8")),
    "v5": json.loads((RES / "universal_control36_v5.json").read_text(encoding="utf-8")),
}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 12,
    "axes.titlesize": 15, "axes.titleweight": "bold", "axes.labelsize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 9.5,
    "axes.grid": True, "grid.alpha": 0.25, "lines.linewidth": 2.4,
    "lines.markersize": 7, "savefig.facecolor": "white",
})
BLUE, SKY, ORANGE, GREEN, RED, GRAY, PURPLE = "#0072B2", "#56B4E9", "#D55E00", "#009E73", "#C1272D", "#666666", "#7B3294"


def agg(version, protocol):
    return V[version]["aggregates"][protocol]


def pattern_bg(version, protocol, pattern):
    vals = [x["off_flip_frac"] for x in V[version]["records"] if x["protocol"] == protocol and x["pattern"] == pattern]
    return 100 * float(np.mean(vals))


fig, ax = plt.subplots(2, 2, figsize=(18, 12))
fig.suptitle("Exp36 — Fiziksel yüzey çekirdeğiyle gürültülü tek-atış atom kontrolü", fontsize=19, fontweight="bold", y=0.995)

# (a) Pareto path: target marginal quality versus overall background disturbance.
a = ax[0, 0]
a.add_patch(Rectangle((1.7, 95), 3.3, 6.5, facecolor=GREEN, alpha=0.10, edgecolor="none", zorder=0))
series = [
    ("v1", "v1 SEG4 μ=.020", GRAY, "SEG-SS1", "SEG-SS5"),
    ("v2", "v2 SEG4 μ=.005", ORANGE, "SEG-SS1", "SEG-SS5"),
    ("v3", "v3 SEG4 μ=.010", RED, "SEG-SS1", "SEG-SS5"),
    ("v4", "v4 aktif güvenlik", PURPLE, "SEG-SS1", "SEG-SS5"),
    ("v5", "v5 SEG6 10 nm", BLUE, "SEG6-SS1", "SEG6-SS5"),
    ("v5", "v5 SEG8 7,5 nm", GREEN, "SEG8-SS1", "SEG8-SS5"),
]
for version, label, color, p1, p5 in series:
    for protocol, marker in ((p1, "o"), (p5, "s")):
        rec = agg(version, protocol)
        x = 100 * rec["all"]["off_flip_frac"]
        y = 100 * rec["nonmono"]["marginal"]
        a.scatter([x], [y], s=125, marker=marker, color=color, edgecolor="black", linewidth=0.8, zorder=5)
    x1 = 100 * agg(version, p1)["all"]["off_flip_frac"]; y1 = 100 * agg(version, p1)["nonmono"]["marginal"]
    x5 = 100 * agg(version, p5)["all"]["off_flip_frac"]; y5 = 100 * agg(version, p5)["nonmono"]["marginal"]
    a.plot([x1, x5], [y1, y5], color=color, alpha=0.45, linewidth=1.4, zorder=2)
a.annotate("SEG6", (100*agg("v5", "SEG6-SS1")["all"]["off_flip_frac"], 100*agg("v5", "SEG6-SS1")["nonmono"]["marginal"]),
           xytext=(-48, -24), textcoords="offset points", fontsize=10, arrowprops=dict(arrowstyle="->", lw=1.0))
a.annotate("SEG8", (100*agg("v5", "SEG8-SS1")["all"]["off_flip_frac"], 100*agg("v5", "SEG8-SS1")["nonmono"]["marginal"]),
           xytext=(10, 8), textcoords="offset points", fontsize=10, arrowprops=dict(arrowstyle="->", lw=1.0))
a.annotate("v2 SEG4", (100*agg("v2", "SEG-SS1")["all"]["off_flip_frac"], 100*agg("v2", "SEG-SS1")["nonmono"]["marginal"]),
           xytext=(8, -4), textcoords="offset points", fontsize=10, arrowprops=dict(arrowstyle="->", lw=1.0))
a.axvline(5, color=RED, linestyle="--", linewidth=1.9)
a.axhline(95, color=RED, linestyle=":", linewidth=1.9)
a.text(1.85, 100.7, "kayıtlı fizibil bölge", color=GREEN, fontsize=10, fontweight="bold")
a.set(xlim=(1.7, 7.2), ylim=(72, 101.5), xlabel="toplam hedef dışı modal çevrilme (%)",
      ylabel="ALT+RAND marjinal bit doğruluğu (%)", title="(a) Pareto yolu: kalite–arka plan dengesi")
handles = [Patch(color=c, label=l) for _, l, c, _, _ in series]
handles += [Line2D([], [], marker="o", color="black", linestyle="none", label="SS1 p=%1"),
            Line2D([], [], marker="s", color="black", linestyle="none", label="SS5 p=%5 çoğunluk")]
a.legend(handles=handles, loc="lower right", ncol=2, framealpha=0.97)

# (b) Final noisy-readout quality by surface-segment pitch.
b = ax[0, 1]
final_series = [
    ("SEG4\nSS1", "v2", "SEG-SS1"), ("SEG4\nSS5", "v2", "SEG-SS5"),
    ("SEG6\nSS1", "v5", "SEG6-SS1"), ("SEG6\nSS5", "v5", "SEG6-SS5"),
    ("SEG8\nSS1", "v5", "SEG8-SS1"), ("SEG8\nSS5", "v5", "SEG8-SS5"),
]
x = np.arange(len(final_series))
marg = [100 * agg(v, p)["nonmono"]["marginal"] for _, v, p in final_series]
exact = [100 * agg(v, p)["nonmono"]["exact_snap"] for _, v, p in final_series]
b1 = b.bar(x - 0.21, marg, width=0.42, color=BLUE, label="marjinal bit doğruluğu")
b2 = b.bar(x + 0.21, exact, width=0.42, color=ORANGE, label="tam desen anlık oranı")
b.bar_label(b1, fmt="%.1f", fontsize=8.5, padding=3)
b.bar_label(b2, fmt="%.1f", fontsize=8.5, padding=3)
for i, (_, v, p) in enumerate(final_series):
    b.text(i, 104.5, agg(v, p)["nonmono"]["modal_count"], ha="center", va="bottom", fontsize=9, fontweight="bold")
b.text(5.45, 108.5, "üstte: modal tam desen", ha="right", fontsize=9)
b.axhline(95, color=RED, linestyle="--", linewidth=1.7, label="marjinal eşiği %95")
b.axhline(65, color=GRAY, linestyle=":", linewidth=1.7, label="tam desen eşiği %65")
b.set(ylim=(0, 112), xticks=x, xticklabels=[s[0] for s in final_series], ylabel="ALT+RAND ortalaması (%)",
      title="(b) Dilim adımına göre gürültülü okuma kalitesi")
b.legend(loc="lower right", framealpha=0.97)

# (c) Background disturbance is concentrated in RAND tasks; the registered budget is overall.
c = ax[1, 0]
patterns = ("MONO", "ALT", "RAND", "TÜM")
bg_series = [
    ("SEG4 SS1", "v2", "SEG-SS1", GRAY),
    ("SEG6 SS1", "v5", "SEG6-SS1", BLUE),
    ("SEG6 SS5", "v5", "SEG6-SS5", SKY),
    ("SEG8 SS1", "v5", "SEG8-SS1", GREEN),
    ("SEG8 SS5", "v5", "SEG8-SS5", PURPLE),
]
x = np.arange(len(patterns)); width = 0.15
for j, (label, version, protocol, color) in enumerate(bg_series):
    vals = [pattern_bg(version, protocol, p) for p in patterns[:3]]
    vals.append(100 * agg(version, protocol)["all"]["off_flip_frac"])
    bars = c.bar(x + (j - 2) * width, vals, width=width, color=color, label=label, edgecolor="black", linewidth=0.3)
    for bar, val in zip(bars, vals):
        if val >= 3:
            c.text(bar.get_x() + bar.get_width()/2, val + 0.35, f"{val:.1f}", ha="center", va="bottom", fontsize=7.5, rotation=90)
c.axhline(5, color=RED, linestyle="--", linewidth=2.0, label="toplam bütçe %5")
c.set(xlim=(-0.55, 3.55), ylim=(0, 15.5), xticks=x, xticklabels=["MONO", "ALT", "RAND", "TÜM\n(kayıtlı ölçüt)"],
      ylabel="hedef dışı modal çevrilme (%)", title="(c) Arka plan bozulması desen sınıfına göre")
c.legend(loc="upper left", ncol=2, framealpha=0.97)

# (d) Depth response and actuator budget for the coarsest passing kernel, SEG6.
d = ax[1, 1]
ds = V["v5"]["depth_stats"]
labels = [f"sığ z<75\nn={ds['shallow_lt75']['n_target_task_bits']}", f"derin z≥75\nn={ds['deep_ge75']['n_target_task_bits']}"]
first = [100 * ds["shallow_lt75"]["first_round_mean"], 100 * ds["deep_ge75"]["first_round_mean"]]
final = [100 * ds["shallow_lt75"]["final_ss1_mean"], 100 * ds["deep_ge75"]["final_ss1_mean"]]
x = np.arange(2)
d1 = d.bar(x - 0.22, first, width=0.44, color=SKY, label="ilk tur (SEG6-ANA tanı değeri)")
d2 = d.bar(x + 0.22, final, width=0.44, color=BLUE, label="son doğrulama (SEG6-SS1)")
d.bar_label(d1, fmt="%.1f", fontsize=10, padding=4)
d.bar_label(d2, fmt="%.1f", fontsize=10, padding=4)
d.axhline(95, color=RED, linestyle="--", linewidth=1.7, label="%95 başvuru çizgisi")
d.set(xlim=(-0.55, 1.55), ylim=(88, 103.5), xticks=x, xticklabels=labels, ylabel="hedef bit marjinali (%)",
      title="(d) SEG6 derinlik yanıtı ve geri besleme kapanışı")
budget = V["v5"]["actuator_budget"]
ss1 = agg("v5", "SEG6-SS1"); ss5 = agg("v5", "SEG6-SS5")
txt = (f"ilk fark {100*ds['first_round_gap_shallow_minus_deep']:.1f} → son fark {100*ds['final_ss1_gap_shallow_minus_deep']:.1f} puan\n"
       f"medyan güncelleme SS1 {ss1['median_updates_nonmono']:.1f}; SS5 {ss5['median_updates_nonmono']:.1f}\n"
       f"≤0,50 V: SEG6 %{100*budget['SEG6']['frac_within_0p50V']:.1f}; SEG8 %{100*budget['SEG8']['frac_within_0p50V']:.1f}")
d.text(0.02, 0.04, txt, transform=d.transAxes, fontsize=9, va="bottom",
       bbox=dict(boxstyle="round,pad=0.45", facecolor="#FFF7E6", edgecolor=ORANGE, alpha=0.97))
d.legend(loc="lower right", framealpha=0.97)

fig.tight_layout(rect=(0, 0, 1, 0.965))
dst = RES / "fig_exp36.png"
try:
    fig.savefig(dst, dpi=250, format="png", pil_kwargs={"compress_level": 6})
    out = dst
except Exception as ex:
    alt = Path.home() / "_exp36_out"
    alt.mkdir(parents=True, exist_ok=True)
    out = alt / "fig_exp36.png"
    fig.savefig(out, dpi=250, format="png", pil_kwargs={"compress_level": 6})
    print(f"WRITE_FAIL {dst} :: {type(ex).__name__}: {ex} -> {out}")
try:
    from PIL import Image
    with Image.open(out) as im:
        print(f"WROTE {out} size={im.size[0]}x{im.size[1]}px bytes={out.stat().st_size}")
except Exception:
    print(f"WROTE {out} bytes={out.stat().st_size}")

