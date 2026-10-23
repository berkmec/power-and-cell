"""Exp33 figure: readable statistics for the 4 K SKY130 cell measurements.
Rendered directly at print resolution (4500x1625 px) from results/exp33/cell33.json,
cell33_sensitivity.json (jump diagnostics) and the stored NIST .mdm sweeps (jump scatter).
No verdicts are recomputed here; the jump detector is the registered C33-4 definition."""
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp32_corrections"))
import nist32  # noqa: E402

D = json.loads((ROOT / "results" / "exp33" / "cell33.json").read_text(encoding="utf-8"))
DATA = ROOT / "results" / "exp32" / "nist_sky130"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "axes.titlesize": 15,
    "axes.titleweight": "bold",
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "lines.linewidth": 2.4,
    "lines.markersize": 7,
    "savefig.facecolor": "white",
})
BLUE, ORANGE, GREEN, RED, GRAY = "#0072B2", "#D55E00", "#009E73", "#C1272D", "#666666"
VB_MODEL = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
DVT_MODEL = np.array([0.0, 0.044876592, 0.1025931, 0.183147286, 0.297640108, 0.315595025])

fig, ax = plt.subplots(1, 3, figsize=(18, 6.5))
fig.suptitle("Exp33 — SKY130 hücresi 4 K: gövde, kink ve basamak denetimi", fontsize=18, fontweight="bold")

# (a) Body effect: measured VT shift versus BSIM4 anchor (reverse VB, magnitude proxy)
a = ax[0]
rows = D["C33_1_body_effect"]["per_device"]
for i, r in enumerate(rows):
    vbs = np.array(r["VB_values"])
    dvt = np.array(r["VT_at_VB"]) - r["VT_at_VB"][-1]
    jitter = (i - (len(rows) - 1) / 2) * 0.012
    color = "#9DBBD3" if r["lvt"] else BLUE
    a.plot(vbs + jitter, dvt, "o-", color=color, alpha=0.75, linewidth=1.4, markersize=5)
for lvt, color, ls, lab in ((False, BLUE, "-", "std medyan (5 aygıt)"), (True, GRAY, "--", "lvt medyan (5 aygıt)")):
    sel = [r for r in rows if r["lvt"] == lvt]
    vbs = np.array(sel[0]["VB_values"])
    med = np.median([np.array(r["VT_at_VB"]) - r["VT_at_VB"][-1] for r in sel], axis=0)
    a.plot(vbs, med, marker="o", color=color, linestyle=ls, linewidth=3.0, markersize=8, label=lab, zorder=5)
vb_line = np.linspace(-1.5, 0, 120)
a.plot(vb_line, np.interp(np.abs(vb_line), VB_MODEL, DVT_MODEL), color=RED, linestyle="--",
       linewidth=2.6, label="BSIM4 300 K çapası (Exp12)")
a.annotate("std |γ| medyan = 0,105\nçapa eğimi ≈ 0,32", (-1.08, 0.205), fontsize=10.5,
           bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=GRAY, alpha=0.95))
a.set(xlim=(-1.62, 0.08), ylim=(-0.10, 0.36), xlabel="gövde gerilimi V_B (V)",
      ylabel="ΔV_T = V_T(V_B) − V_T(0) (V)", title="(a) Gövde etkisi: ölçüm BSIM çapasının altında")
a.legend(loc="lower left", framealpha=0.98)

# (b) Kink amplification per device
b = ax[1]
ka = sorted(D["C33_3_kink_amplification"]["rows_VG09"], key=lambda r: r["A_amplification"])
def dev_label(name):
    m = re.search(r"_w([0-9p]+)u_l([0-9p]+)u", name)
    return f"{'lvt' if '_lvt_' in name else 'std'} W{m.group(1).replace('p', ',')}/L{m.group(2).replace('p', ',')}"
labels = [dev_label(r["file"]) for r in ka]
vals = [r["A_amplification"] for r in ka]
colors = [BLUE if r["IB_1p8_A"] > 1e-9 else "#B8C4CC" for r in ka]
bars = b.barh(np.arange(len(ka)), vals, color=colors, edgecolor="black", linewidth=0.5)
for bar, v in zip(bars, vals):
    b.text(v * 1.08, bar.get_y() + bar.get_height() / 2, f"{v:.0f}", va="center", fontsize=9.5)
median_A = D["C33_3_kink_amplification"]["median_A_VG09"]
b.axvline(median_A, color=RED, linestyle="--", linewidth=2.2, label=f"medyan (hüküm kümesi) = {median_A:.0f}")
b.axvline(30, color=GRAY, linestyle=":", linewidth=1.8, label="kayıtlı alt sınır 30")
b.scatter([], [], color="#B8C4CC", edgecolor="black", s=90, label="I_B≤1 nA (hüküm dışı)")
b.set(xscale="log", xlim=(20, 14_000), yticks=np.arange(len(ka)), yticklabels=labels,
      xlabel="kink yükseltimi A = ΔI_D / |I_B| @1,8 V", title="(b) Kink yükseltimi: aygıt aygıt (V_G=0,9 V)")
b.tick_params(axis="y", labelsize=9.5)
b.legend(loc="lower right", framealpha=0.96)

# (c) Jump audit over every stored sweep (registered C33-4 detector)
c = ax[2]
groups = []
for x0, sub, pat, label, color in ((0, "4K", "*IDVD*_4K.mdm", "4 K", BLUE),
                                   (1, "300K", "*IDVD*.mdm", "300 K kontrol", ORANGE)):
    js = []
    for f in sorted((DATA / sub).glob(pat)):
        for var, cc in nist32.read_mdm(f):
            if abs(var.get("VB", 0)) > 1e-9:
                continue
            vd, idd = np.asarray(cc["VD"], float), np.abs(np.asarray(cc["ID"], float))
            ok = (idd[:-1] > 1e-11) & (vd[:-1] >= 0.3) & (idd[1:] > 0)
            if ok.any():
                js.append(float(np.max(np.log10(idd[1:][ok] / idd[:-1][ok]))))
    js = np.array(sorted(js))
    xj = x0 + np.linspace(-0.20, 0.20, len(js))
    c.scatter(xj, js, s=72, color=color, edgecolor="black", linewidth=0.6, label=f"{label} ({len(js)} tarama)", zorder=5)
    groups.append((x0, js))
c.axhline(0.5, color=RED, linestyle="--", linewidth=2.2, label="kayıtlı kilitlenme sınırı 0,5 dex")
c.axhline(0.3, color=GRAY, linestyle=":", linewidth=1.8)
c.annotate("4 K en büyük: 0,386 dex\n33 nA, kink bölgesi", (0.20, 0.386), xytext=(0.28, 0.46),
           fontsize=10.5, arrowprops=dict(arrowstyle="->", lw=1.3))
c.annotate("300 K: 0,641 dex\n0,45 nA taban artefaktı", (1.20, 0.641), xytext=(0.48, 0.64),
           fontsize=10.5, arrowprops=dict(arrowstyle="->", lw=1.3))
c.set(xlim=(-0.55, 1.55), ylim=(-0.03, 0.72), xticks=[0, 1], xticklabels=["4 K", "300 K"],
      ylabel="en büyük bitişik basamak J (dex)", title="(c) Basamak denetimi: bütün taramalar")
c.legend(loc="center", bbox_to_anchor=(0.50, 0.30), framealpha=0.98)

fig.tight_layout(rect=(0, 0, 1, 0.94))
dst = ROOT / "results" / "exp33" / "fig_exp33.png"
try:
    fig.savefig(dst, dpi=250, format="png", pil_kwargs={"compress_level": 6})
    out = dst
except Exception as ex:
    alt = Path.home() / "_cell33_out"
    alt.mkdir(parents=True, exist_ok=True)
    out = alt / "fig_exp33.png"
    fig.savefig(out, dpi=250, format="png", pil_kwargs={"compress_level": 6})
    print(f"WRITE_FAIL {dst} :: {type(ex).__name__}: {ex} -> {out}")
try:
    from PIL import Image
    with Image.open(out) as im:
        print(f"WROTE {out} size={im.size[0]}x{im.size[1]}px bytes={out.stat().st_size}")
except Exception:
    print(f"WROTE {out} bytes={out.stat().st_size}")
