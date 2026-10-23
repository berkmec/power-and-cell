"""Exp34 figure: readable statistics for the ultra-thrust envelope.
Rendered directly at print resolution (4500x3000 px) from results/exp34/thrust34.json;
measured-thruster markers from results/exp28/rocket.json. No new physics is computed."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
T = json.loads((ROOT / "results" / "exp34" / "thrust34.json").read_text(encoding="utf-8"))
R = json.loads((ROOT / "results" / "exp28" / "rocket.json").read_text(encoding="utf-8"))

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
fig, ax = plt.subplots(2, 2, figsize=(18, 12))
fig.suptitle("Exp34 — Ultra iyon itki zarfı (Exp0–Exp33 verisi)", fontsize=19, fontweight="bold")

# (a) Pareto front: thrust/power versus specific impulse
a = ax[0, 0]
for key, label, color in (("pareto_EMI_s1", "EMI⁺ fizibil cephe", BLUE),
                          ("pareto_In_s1", "In⁺ fizibil cephe", ORANGE)):
    pts = T[key]
    a.plot([p["Isp_s"] for p in pts], [p["FP_uN_W"] for p in pts], "o-", color=color, label=label)
for key, color in (("EMI+_precharged", BLUE), ("In+_neutral", ORANGE)):
    o = T["optima_s1"][key]["unconstr"]
    a.scatter([o["Isp_s"]], [o["FP_uN_W"]], marker="*", s=260, color=color, edgecolor="black", zorder=6)
a.scatter([], [], marker="*", s=180, color=GRAY, edgecolor="black", label="kısıtsız tepe (V*)")
meas = R["measured"]
a.scatter([m["Isp_s"] for m in meas], [m["F_over_P_uN_W"] for m in meas],
          marker="^", s=105, color="black", label="ölçülen itkiciler", zorder=7)
a.axhline(1e6 / 299_792_458, color="black", linestyle="--", linewidth=1.8, label="foton sınırı 1/c")
emi_f = T["optima_s1"]["EMI+_precharged"]["constr_beta1"]
a.annotate(f"EMI⁺ fizibil\n{emi_f['FP_uN_W']:.0f} µN/W @ {emi_f['Isp_s']:.0f} s",
           (emi_f["Isp_s"], emi_f["FP_uN_W"]), xytext=(18, 14), textcoords="offset points",
           fontsize=10.5, arrowprops=dict(arrowstyle="->", lw=1.4))
knee = T["knee_In"]["pareto_point"]
a.annotate(f"In diz noktası\n{knee['FP_uN_W']:.0f} µN/W @ {knee['Isp_s']:.0f} s",
           (knee["Isp_s"], knee["FP_uN_W"]), xytext=(-135, 18), textcoords="offset points",
           fontsize=10.5, arrowprops=dict(arrowstyle="->", lw=1.4))
a.set(xscale="log", yscale="log", xlim=(100, 12_000), ylim=(0.001, 2500),
      xlabel="özgül itki Isp (s)", ylabel="itki/güç F/P (µN/W)",
      title="(a) Pareto cephesi: verim–özgül itki dengesi")
a.legend(loc="lower left", framealpha=0.96)

# (b) Species optima: unconstrained versus emission-feasible at beta=1
b = ax[0, 1]
species = [("EMI+_precharged", "EMI⁺"), ("Cs+_neutral", "Cs⁺"), ("In+_neutral", "In⁺"),
           ("Ga+_neutral", "Ga⁺"), ("Xe+_neutral", "Xe⁺"), ("O-_precharged", "O⁻")]
x = np.arange(len(species))
wu = [T["optima_s1"][k]["unconstr"]["FP_uN_W"] for k, _ in species]
wc = [T["optima_s1"][k]["constr_beta1"]["FP_uN_W"] for k, _ in species]
b1 = b.bar(x - 0.19, wu, width=0.38, color=BLUE, label="kısıtsız tepe")
b2 = b.bar(x + 0.19, wc, width=0.38, color=ORANGE, label="fizibil (d≥4,148 nm, β=1)")
b.bar_label(b1, fmt="%.0f", fontsize=9, padding=3)
b.bar_label(b2, fmt="%.0f", fontsize=9, padding=3)
b.axhline(40, color="black", linestyle="--", linewidth=1.8, label="NSTAR ölçülen 40 µN/W")
b.set(yscale="log", ylim=(10, 900), xticks=x, xticklabels=[s[1] for s in species],
      ylabel="F/P (µN/W)", title="(b) Tür optimumları: kısıtsız ve fizibil")
b.legend(loc="lower left", framealpha=0.96)

# (c) Cell supply current across the Exp17 grid and body-effect sensitivity band
c = ax[1, 0]
sets = [("s1_Exp32", "s=1,00 Exp32 (9/9 bölünür)", BLUE),
        ("s0p50_sens", "s=0,50 keşif (6/9)", ORANGE),
        ("s0p33_sens", "s=0,33 keşif (1/9)", GREEN)]
keys = list(T["cell_grids"]["s1_Exp32"].keys())
labels = [k.replace("rD", "").replace("_lam", "/") for k in keys]
x = np.arange(len(keys))
for j, (skey, lab, color) in enumerate(sets):
    vals = []
    for k in keys:
        v = T["cell_grids"][skey][k]["I_c_mA"]
        vals.append(np.nan if v is None else v)
    bars = c.bar(x + (j - 1) * 0.27, vals, width=0.27, color=color, label=lab)
    for xi, v in zip(x + (j - 1) * 0.27, vals):
        if np.isnan(v):
            c.scatter([xi], [0.08], marker="x", s=95, color=RED, linewidths=2.2, zorder=6)
c.scatter([], [], marker="x", s=95, color=RED, linewidths=2.2, label="bölünme yok")
c.set(xlim=(-0.6, len(keys) - 0.4), ylim=(0, 2.9), xticks=x, xticklabels=labels,
      xlabel="balast ızgarası  rD / λ (µm)", ylabel="hücre akımı I_c (mA)",
      title="(c) Hücre beslemesi: I_c ve gövde-etkisi bandı")
c.tick_params(axis="x", rotation=32)
for t in c.get_xticklabels():
    t.set_ha("right")
c.legend(loc="upper right", framealpha=0.96)

# (d) Reference mission operating points
d = ax[1, 1]
missions = T["missions_10kg_100W_2kg"]
names = ["EMI⁺ 8,4 V", "In⁺ diz 1,20 kV", "In⁺ 2 kV", "In⁺ 5 kV", "Xe⁺ 1,8 kV"]
colors = [BLUE, ORANGE, ORANGE, RED, GREEN]
offsets = [(10, 10), (-105, -34), (12, 10), (-95, 4), (14, -30)]
for m, name, color, off in zip(missions, names, colors, offsets):
    d.scatter(m["burn_y_eta06"], m["dv_km_s"], s=115, color=color, edgecolor="black", zorder=5)
    d.annotate(f"{name}\n{m['FP_uN_W']:.0f} µN/W", (m["burn_y_eta06"], m["dv_km_s"]),
               xytext=off, textcoords="offset points", fontsize=10,
               arrowprops=dict(arrowstyle="->", lw=1.1))
d.axvline(5, color=GRAY, linestyle=":", linewidth=1.8)
d.axhline(10, color=GRAY, linestyle=":", linewidth=1.8)
d.text(5.12, 1.0, "5 yıl", fontsize=10.5, color=GRAY)
d.text(0.009, 10.45, "Δv = 10 km/s", fontsize=10.5, color=GRAY)
d.set(xscale="log", xlim=(0.005, 7), ylim=(0, 22.5),
      xlabel="yanma süresi (yıl, 100 W; η=0,6 referansı)", ylabel="Δv (km/s)",
      title="(d) Referans görev: 10 kg, 100 W, 2 kg yakıt")

fig.tight_layout(rect=(0, 0, 1, 0.965))
dst = ROOT / "results" / "exp34" / "fig_exp34.png"
try:
    fig.savefig(dst, dpi=250, format="png", pil_kwargs={"compress_level": 6})
    out = dst
except Exception as ex:
    alt = Path.home() / "_exp34_out"
    alt.mkdir(parents=True, exist_ok=True)
    out = alt / "fig_exp34.png"
    fig.savefig(out, dpi=250, format="png", pil_kwargs={"compress_level": 6})
    print(f"WRITE_FAIL {dst} :: {type(ex).__name__}: {ex} -> {out}")
try:
    from PIL import Image
    with Image.open(out) as im:
        print(f"WROTE {out} size={im.size[0]}x{im.size[1]}px bytes={out.stat().st_size}")
except Exception:
    print(f"WROTE {out} bytes={out.stat().st_size}")
