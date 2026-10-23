import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = json.loads(Path(__file__).resolve().parents[2].joinpath("results", "exp39", "ship39.json").read_text(encoding="utf-8"))
cont = D["continuous_ladder"]
rungs = {int(r["n"]): r for r in D["rungs"]}
n = np.array([r["n"] for r in cont])
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
fig.suptitle("Exp39 — Denklem 1 (tam bütün): $\\varepsilon=0{,}1^n$, $\\gamma=10^n$ nonilyon merdiveni", fontsize=16, fontweight="bold")

ax = axes[0, 0]
ax.plot(n, [r["log10_gamma"] for r in cont], color="#1f4e79", lw=2.4, label="$\\log_{10}\\gamma=n$")
ax.set_xlabel("n")
ax.set_ylabel("$\\log_{10}\\gamma$", color="#1f4e79")
ax2 = ax.twinx()
ax2.plot(n, [r["log10_one_minus_beta"] for r in cont], color="#c00000", lw=2.2, label="$\\log_{10}(1-\\beta)$")
ax2.set_ylabel("$\\log_{10}(1-\\beta)$", color="#c00000")
ax.scatter([1, 30], [1, 30], color="#1f4e79", zorder=5)
ax.annotate("n=1: $\\gamma=10$\n$\\beta=0{,}994987$", (1, 1), xytext=(3, 4), arrowprops=dict(arrowstyle="->"), fontsize=9)
ax.annotate("n=30: $\\gamma=10^{30}$\n= 1 nonilyon\n$1-\\beta=5\\times10^{-61}$", (30, 30), xytext=(17, 20), arrowprops=dict(arrowstyle="->"), fontsize=9)
ax.set_title("(a) Hız c'ye asimptotik yaklaşır; katsayı $10^n$ büyür")
ax.grid(alpha=.25)

ax = axes[0, 1]
lanes = [
    ("EMI_R1_feasible", "Exp34 EMI fizibil", "#9467bd"),
    ("In_knee_R2", "Exp34 In diz (44,8 km/s)", "#1f77b4"),
    ("In_5kV_formula", "In 5 kV (91,7 km/s)", "#2ca02c"),
    ("photon_lower_bound", "Foton alt sınırı", "#ff7f0e"),
]
for key, label, color in lanes:
    vals = np.array([r["mass_ratio_log10"][key] for r in cont])
    ax.plot(n, np.log10(1.0 + vals), label=label, color=color, lw=2.1)
ax.axhline(math.log10(1.0 + math.log10(10/8)), color="black", ls="--", lw=1.5, label="Referans gemi 10/8")
ax.set_xlabel("n")
ax.set_ylabel("$\\log_{10}(M_{ıslak}/M_{kuru})$ (çift-log eksen)")
_ticks = [0, 1, 10, 100, 1000, 10000, 100000, 1000000]
ax.set_yticks([math.log10(1.0 + v) for v in _ticks])
ax.set_yticklabels([f"{v:g}" for v in _ticks])
ax.set_title("(b) Relativistik roket kütle-oranı bedeli")
ax.legend(fontsize=8, loc="upper left")
ax.grid(alpha=.25)
ax.text(.98, .08, "n=1:\nIn diz log10R=8 696\nfoton R=19,95", transform=ax.transAxes, ha="right", fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=.85))

ax = axes[1, 0]
mask = np.array([r["log10_kinetic_energy_J"] is not None for r in cont])
ax.plot(n[mask], [r["log10_kinetic_energy_J"] for r, m in zip(cont, mask) if m], color="#1f4e79", lw=2.3, label="$\\log_{10}E_{kin}$ [J]")
ax.set_xlabel("n")
ax.set_ylabel("$\\log_{10}E_{kin}$ [J]", color="#1f4e79")
ax2 = ax.twinx()
ax2.plot(n, [math.log10(r["in_wavelength_over_planck"]) for r in cont], color="#c00000", lw=2.1, label="$\\log_{10}(\\lambda_{In}/\\ell_P)$")
ax2.axhline(0, color="#c00000", ls=":", lw=1.2)
ax2.set_ylabel("$\\log_{10}(\\lambda_{In}/\\ell_P)$", color="#c00000")
ax.axvline(D["planck_crossing"]["n"], color="gray", ls="--", lw=1.4)
ax.text(D["planck_crossing"]["n"] + .35, 24, "Planck geçişi\nn=17,856", fontsize=9)
ax.set_title("(c) 8 kg kuru gemi enerjisi ve tek In iyonu $h\\nu$ boyu")
ax.grid(alpha=.25)

ax = axes[1, 1]
sel = [r for r in cont if r["n"] >= 0.1]
x = np.array([r["n"] for r in sel])
ax.plot(x, [math.log10(r["ideal_100W_supply_years"]) for r in sel], color="#c00000", lw=2.3, label="İdeal 100 W enerji tedarik süresi")
ax.plot(x, [math.log10(r["alpha_centauri_ship_proper_years"]) for r in sel], color="#1f77b4", lw=2.3, label="Alfa Centauri gemi öz süresi (seyir)")
ax.axhline(math.log10(4.37), color="#2ca02c", ls="--", lw=1.7, label="Dünya çerçevesi seyir ≈4,37 yıl")
ax.set_xlabel("n")
ax.set_ylabel("$\\log_{10}$(yıl)")
ax.set_title("(d) Zaman ölçekleri: tedarik büyür, öz süre küçülür")
ax.legend(fontsize=8, loc="center left")
ax.grid(alpha=.25)
ax.text(.03, .05, "Gerçek Exp34 yanışları (10→8 kg):\nEMI Δv=0,85 km/s; In diz 10,0 km/s; In 5 kV 20,4 km/s\n→ hiçbiri n=1 basamağına ulaşmıyor.",
        transform=ax.transAxes, fontsize=9, bbox=dict(boxstyle="round", facecolor="white", alpha=.88))

fig.tight_layout(rect=[0, 0, 1, .94])
out = Path(__file__).resolve().parents[2].joinpath("results", "exp39", "fig_exp39.png")
try:
    fig.savefig(out, dpi=300)
except PermissionError:
    out = Path.home() / "_fig39.png"
    fig.savefig(out, dpi=300)
print(out, out.stat().st_size)
