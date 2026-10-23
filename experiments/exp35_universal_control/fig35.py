"""Exp35 figure: universal cell-atom control statistics.
Rendered directly at print resolution (4500x3000 px) from the preserved v2/v3 JSON results.
No simulation or verdict is recomputed; pooled non-monotone values are task-level means from stored tasks."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

_here = Path(__file__).resolve()
ROOT = _here.parents[2] if (_here.parents[2] / "results" / "exp35").exists() else Path(__file__).resolve().parents[2]
V2 = json.loads((ROOT / "results" / "exp35" / "universal_control35.json").read_text(encoding="utf-8"))
V3 = json.loads((ROOT / "results" / "exp35" / "universal_control35_v3.json").read_text(encoding="utf-8"))

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 12,
    "axes.titlesize": 15, "axes.titleweight": "bold", "axes.labelsize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 9.5,
    "axes.grid": True, "grid.alpha": 0.25, "lines.linewidth": 2.4,
    "lines.markersize": 7, "savefig.facecolor": "white",
})

BLUE, ORANGE, GREEN, RED, GRAY, PURPLE = "#0072B2", "#D55E00", "#009E73", "#C1272D", "#666666", "#7B3294"
SERIES = [
    ("G", "v2 G", GRAY, V2),
    ("L20", "v2 L20", BLUE, V2),
    ("CL20", "v2 CL20", "#56B4E9", V2),
    ("L20M", "v3 L20M", ORANGE, V3),
    ("CL20M", "v3 CL20M", GREEN, V3),
    ("CL20A", "v3 CL20A", PURPLE, V3),
]
PATTERNS = ("MONO", "ALT", "RAND")


def tasks(data):
    return [t for island in data["islands"] for t in island["tasks"]]


def pooled(data, protocol, key, patterns=("ALT", "RAND")):
    vals = [t["protocols"][protocol][key] for t in tasks(data) if t["pattern"] in patterns]
    return float(np.mean(vals))


fig, ax = plt.subplots(2, 2, figsize=(18, 12))
fig.suptitle("Exp35 — Evrensel hücre–atom kontrolü: küresel sürüşten kapalı çevrime", fontsize=19, fontweight="bold", y=0.995)

# (a) Modal exact-pattern success by pattern class and controller generation.
a = ax[0, 0]
x = np.arange(len(PATTERNS))
width = 0.125
for j, (protocol, label, color, data) in enumerate(SERIES):
    vals, counts = [], []
    for pat in PATTERNS:
        rec = data["aggregates"][protocol]["by_pattern"][pat]
        vals.append(100 * rec["modal_exact_fraction"])
        counts.append(rec["modal_exact_count"])
    bars = a.bar(x + (j - 2.5) * width, vals, width=width, color=color, label=label, edgecolor="black", linewidth=0.35)
    for bar, cnt in zip(bars, counts):
        a.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.6, cnt,
               ha="center", va="bottom", fontsize=8, rotation=90)
a.set(xlim=(-0.55, 2.55), ylim=(0, 116), xticks=x, xticklabels=["MONO\n(eşik uyumlu)", "ALT\n(değişken)", "RAND\n(rastgele)"],
      ylabel="modal tam desen başarısı (%)", title="(a) Desen sınıfına göre tam yazma başarısı")
# Series legend is placed once below the figure title, outside panel (a), to keep count labels visible.

# (b) Non-monotone quality: marginal bit accuracy versus exact snapshot fraction.
b = ax[0, 1]
labels = [s[1] for s in SERIES]
x = np.arange(len(SERIES))
bit_vals = [100 * pooled(data, protocol, "mean_bit_accuracy") for protocol, _, _, data in SERIES]
exact_vals = [100 * pooled(data, protocol, "exact_snapshot_fraction") for protocol, _, _, data in SERIES]
b1 = b.bar(x - 0.21, bit_vals, width=0.42, color=BLUE, label="ortalama bit doğruluğu")
b2 = b.bar(x + 0.21, exact_vals, width=0.42, color=ORANGE, label="tam desen anlık oranı")
b.bar_label(b1, fmt="%.0f", fontsize=9, padding=3)
b.bar_label(b2, fmt="%.0f", fontsize=9, padding=3)
b.axhline(95, color=RED, linestyle="--", linewidth=1.8, label="kayıtlı marjinal eşiği %95")
b.axhline(65, color=GRAY, linestyle=":", linewidth=1.8, label="kayıtlı tam desen eşiği %65")
b.set(ylim=(0, 112), xticks=x, xticklabels=labels, ylabel="ALT+RAND ortalaması (%)",
      title="(b) Zor desenlerde bit ve tam-desen kalitesi")
b.tick_params(axis="x", rotation=24)
for t in b.get_xticklabels():
    t.set_ha("right")
b.legend(loc="lower right", framealpha=0.97)

# (c) Register scaling under the independent-bit bootstrap reference.
c = ax[1, 0]
boot = V2["bootstrap_scaling"]
Ns = np.array([10, 30, 100], dtype=float)
styles = [
    ("calibrated_60nm", "one_shot_exact_mean", "kalibreli, tek atış", BLUE, "o", "-"),
    ("calibrated_60nm", ("retry3", "exact_completion_mean"), "kalibreli, ≤3 doğrulama denemesi", BLUE, "s", "-"),
    ("uncalibrated_100nm_T4K", "one_shot_exact_mean", "kalibresiz, tek atış", RED, "o", "--"),
    ("uncalibrated_100nm_T4K", ("retry3", "exact_completion_mean"), "kalibresiz, ≤3 deneme", RED, "s", "--"),
]
for pool, key, label, color, marker, ls in styles:
    vals = []
    for N in Ns.astype(int):
        rec = boot[pool]["by_register_size"][str(N)]
        vals.append(rec[key[0]][key[1]] if isinstance(key, tuple) else rec[key])
    c.plot(Ns, vals, marker=marker, linestyle=ls, color=color, label=label, markersize=8)
c.axhline(0.90, color=GREEN, linestyle=":", linewidth=2.0, label="tam kayıt hedefi %90")
c.set(xscale="log", yscale="log", xlim=(8, 125), ylim=(8e-10, 1.25), xticks=Ns, xticklabels=["10", "30", "100"],
      xlabel="kayıt boyutu N (bit)", ylabel="tam kayıt başarı olasılığı",
      title="(c) Ölçekleme: kalibrasyon ve doğrulama denemesi")
c.legend(loc="lower left", framealpha=0.97)

# (d) Target success comes with background disturbance: exactness versus non-target flips.
d = ax[1, 1]
offsets = {"v2 G": (8, 12), "v2 L20": (8, 10), "v2 CL20": (8, 10),
           "v3 L20M": (-86, 10), "v3 CL20M": (-108, -28), "v3 CL20A": (-82, 10)}
for protocol, label, color, data in SERIES:
    xv = 100 * pooled(data, protocol, "nontarget_modal_flip_fraction_vs_G")
    yv = 100 * pooled(data, protocol, "exact_snapshot_fraction")
    d.scatter([xv], [yv], s=130, color=color, edgecolor="black", zorder=5)
    d.annotate(label, (xv, yv), xytext=offsets[label], textcoords="offset points", fontsize=10.5,
               arrowprops=dict(arrowstyle="->", lw=1.1))
d.axvline(5, color=RED, linestyle="--", linewidth=2.0, label="kayıtlı arka-plan bozulma sınırı %5")
d.axhline(65, color=GRAY, linestyle=":", linewidth=1.8, label="tam desen eşiği %65")
d.set(xlim=(-0.4, 7.8), ylim=(-5, 108), xlabel="hedef dışı modal çevrilme oranı (%)",
      ylabel="ALT+RAND tam desen anlık oranı (%)", title="(d) Tam kontrolün bedeli: arka-plan bozulması")
d.legend(loc="lower right", framealpha=0.97)

fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.legend(handles=[Patch(color=color, label=label) for _, label, color, _ in SERIES],
           loc="upper center", bbox_to_anchor=(0.5, 0.972), ncol=6, frameon=False)
dst = ROOT / "results" / "exp35" / "fig_exp35.png"
try:
    fig.savefig(dst, dpi=250, format="png", pil_kwargs={"compress_level": 6})
    out = dst
except Exception as ex:
    alt = Path.home() / "_exp35_out"
    alt.mkdir(parents=True, exist_ok=True)
    out = alt / "fig_exp35.png"
    fig.savefig(out, dpi=250, format="png", pil_kwargs={"compress_level": 6})
    print(f"WRITE_FAIL {dst} :: {type(ex).__name__}: {ex} -> {out}")
try:
    from PIL import Image
    with Image.open(out) as im:
        print(f"WROTE {out} size={im.size[0]}x{im.size[1]}px bytes={out.stat().st_size}")
except Exception:
    print(f"WROTE {out} bytes={out.stat().st_size}")
