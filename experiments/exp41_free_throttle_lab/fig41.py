import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
D = json.loads((ROOT / "results" / "exp41" / "lab41.json").read_text(encoding="utf-8"))
OUT = ROOT / "results" / "exp41"
FINITE = D["finite_ship"]
EXT = D["external_bench"]
EVENTS = D["external_events"]["1g_FULL"]
LANES = {
    "EMI_R1_feasible": ("EMI fizibil", "#9467bd"),
    "In_knee_R2": ("In diz", "#1f77b4"),
    "In_5kV_formula": ("In 5 kV", "#2ca02c"),
}
PROFILES = {
    "FULL": ("Tam gaz", "#d62728"),
    "ECO25": ("%25 gaz", "#1f77b4"),
    "RAMP": ("Rampa", "#2ca02c"),
    "PULSE50": ("Darbe %50", "#ff7f0e"),
}


def save(fig, name):
    path = OUT / name
    try:
        fig.savefig(path, dpi=300)
    except PermissionError:
        path = Path.home() / name
        fig.savefig(path, dpi=300)
    plt.close(fig)
    print(path, path.stat().st_size)


def farr(lane, profile, field):
    return np.array([r[field] for r in FINITE[lane][profile]["samples"]], dtype=float)


def earr(run, field):
    return np.array([r[field] for r in EXT[run]["samples"]], dtype=float)


def fig_a():
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle("Exp41A — Sonlu yakıtlı laboratuvar: hız serbest, gaz profili yalnız saati değiştirir", fontsize=16, fontweight="bold")
    ax = axes[0, 0]
    for p, (lab, col) in PROFILES.items():
        ax.plot(farr("In_knee_R2", p, "earth_years"), farr("In_knee_R2", p, "delta_v_km_s"), label=lab, color=col, lw=2.1)
    ax.set(xlabel="Dünya zamanı [yıl]", ylabel="Δv [km/s]", title="(a) In diz: serbest hızın zamanla evrimi")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=.25)
    ax.text(.03, .92, "Bütün profiller aynı 2 kg yakıtı yakar\n→ son Δv aynı: 10,001 km/s", transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=.9))

    ax = axes[0, 1]
    for p, (lab, col) in PROFILES.items():
        ax.plot(farr("In_knee_R2", p, "earth_years"), farr("In_knee_R2", p, "mass_kg"), label=lab, color=col, lw=2.1)
    ax.set(xlabel="Dünya zamanı [yıl]", ylabel="Gemi kütlesi [kg]", title="(b) In diz: yakıt tüketimi ve gaz profili")
    ax.legend(fontsize=8)
    ax.grid(alpha=.25)

    ax = axes[1, 0]
    for lane, (lab, col) in LANES.items():
        ax.plot(farr(lane, "FULL", "mass_kg"), farr(lane, "FULL", "delta_v_km_s"), label=lab, color=col, lw=2.2)
    ax.invert_xaxis()
    ax.set(xlabel="Kalan kütle [kg] (10 → 8)", ylabel="Δv [km/s]", title="(c) Tam gazda kütle–hız yörüngesi")
    ax.legend(fontsize=8)
    ax.grid(alpha=.25)

    ax = axes[1, 1]
    lanes = list(LANES)
    x = np.arange(len(lanes))
    width = .19
    for j, profile in enumerate(PROFILES):
        vals = [FINITE[l][profile]["final"]["earth_years"] for l in lanes]
        bars = ax.bar(x + (j - 1.5) * width, vals, width=width, label=PROFILES[profile][0], color=PROFILES[profile][1])
        ax.bar_label(bars, fmt="%.3g", fontsize=7, padding=1)
    ax.set_yscale("log")
    ax.set_xticks(x, [LANES[l][0] for l in lanes])
    ax.set(ylabel="Yanış süresi [yıl] (log)", title="(d) Gaz profili son hızı değil, yanış süresini değiştirir")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=.25)
    fig.tight_layout(rect=[0, 0, 1, .94])
    save(fig, "fig_exp41a_finite_lab.png")


def fig_b():
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle("Exp41B — Dış beslemeli gaz tezgâhı: hız serbest, c asimptotik sınır", fontsize=16, fontweight="bold")
    runs = [("0.1g_FULL", "0,1g", "#1f77b4"), ("1g_FULL", "1g", "#d62728"), ("10g_FULL", "10g", "#2ca02c")]
    ax = axes[0, 0]
    for run, lab, col in runs:
        ax.plot(earr(run, "proper_years"), earr(run, "log10_gamma"), label=lab, color=col, lw=2.2)
    ax.set(xlabel="Gemi öz zamanı [yıl]", ylabel="$\\log_{10}\\gamma$", title="(a) Serbest hızın γ okuması (hedef değil)", xlim=(0, 85))
    ax.legend(fontsize=9)
    ax.grid(alpha=.25)

    ax = axes[0, 1]
    for run, lab, col in runs:
        ax.plot(earr(run, "proper_years"), earr(run, "log10_one_minus_beta"), label=lab, color=col, lw=2.2)
    ax.set(xlabel="Gemi öz zamanı [yıl]", ylabel="$\\log_{10}(1-\\beta)$", title="(b) c'ye kalan hız açığı", xlim=(0, 85))
    ax.legend(fontsize=9)
    ax.grid(alpha=.25)

    ax = axes[1, 0]
    for run, lab, col in runs:
        ax.plot(earr(run, "proper_years"), np.log10(1 + earr(run, "earth_years")), label=lab, color=col, lw=2.2)
    ax.set(xlabel="Gemi öz zamanı [yıl]", ylabel="$\\log_{10}(1+t_{Dünya}/yıl)$", title="(c) Öz zaman – Dünya zamanı ayrışması", xlim=(0, 85))
    ax.legend(fontsize=9)
    ax.grid(alpha=.25)

    ax = axes[1, 1]
    for run, lab, col in runs:
        ax.plot(np.log10(1 + earr(run, "earth_years")), np.log10(1 + earr(run, "distance_ly")), label=lab, color=col, lw=2.2)
    lim = [0, 35]
    ax.plot(lim, lim, color="black", ls="--", lw=1.4, label="ışık konisi x=ct")
    ax.set(xlabel="$\\log_{10}(1+t_{Dünya}/yıl)$", ylabel="$\\log_{10}(1+x/[ışık yılı])$", title="(d) Dünya çizgisi ışık konisinin altında kalır", xlim=(0, 35), ylim=(0, 35))
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=.25)
    fig.tight_layout(rect=[0, 0, 1, .94])
    save(fig, "fig_exp41b_external_bench.png")


def fig_c():
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle("Exp41C — Serbest koşuda olaylar, enerji faturası ve In `hν` göstergesi", fontsize=16, fontweight="bold")
    ax = axes[0, 0]
    run = "1g_FULL"
    ax.plot(earr(run, "proper_years"), earr(run, "log10_gamma"), color="#1f4e79", lw=2.3, label="$\\log_{10}\\gamma$")
    ax.set(xlabel="Gemi öz zamanı [yıl]", ylabel="$\\log_{10}\\gamma$", title="(a) 1g tam gaz: γ ve In dalga boyu göstergesi")
    ax2 = ax.twinx()
    ax2.plot(earr(run, "proper_years"), np.log10(earr(run, "in_lambda_over_planck")), color="#c00000", lw=2.1, label="$\\log_{10}(\\lambda_{In}/\\ell_P)$")
    ax2.axhline(0, color="#c00000", ls=":", lw=1.2)
    ax2.set_ylabel("$\\log_{10}(\\lambda_{In}/\\ell_P)$", color="#c00000")
    for key, txt in [("gamma_10", "γ=10"), ("in_lambda_planck_length", "In λ=ℓ_P"), ("gamma_1e30_readout", "γ=10$^{30}$ okuması")]:
        x = EVENTS[key]["proper_years"]
        ax.axvline(x, color="gray", ls="--", lw=1.1)
        ax.text(x + .7, 5, txt, rotation=90, fontsize=8, va="bottom")
    ax.grid(alpha=.25)

    ax = axes[0, 1]
    labels = [
        ("gamma_10", "γ=10"),
        ("in_lambda_1e-18_m", "In λ=10$^{-18}$ m"),
        ("in_planck_energy", "In E=E_P"),
        ("in_lambda_equals_schwarzschild", "In λ=r_s"),
        ("in_lambda_planck_length", "In λ=ℓ_P"),
        ("gamma_1e30_readout", "γ=10$^{30}$"),
    ]
    y = np.arange(len(labels))
    proper = [math.log10(EVENTS[k]["proper_years"]) for k, _ in labels]
    earth = [math.log10(EVENTS[k]["earth_years"]) for k, _ in labels]
    ax.barh(y + .19, proper, height=.36, label="gemi öz zamanı", color="#1f77b4")
    ax.barh(y - .19, earth, height=.36, label="Dünya zamanı", color="#d62728")
    ax.set_yticks(y, [lab for _, lab in labels])
    ax.invert_yaxis()
    ax.set(xlabel="$\\log_{10}$(yıl)", title="(b) 1g olay zaman çizelgesi (olaylar hedef değil, okuma)")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(axis="x", alpha=.25)

    ax = axes[1, 0]
    x = earr(run, "log10_gamma")
    ekin = np.log10(np.maximum(earr(run, "kinetic_energy_J"), 1e-300))
    power = np.log10(np.maximum(earr(run, "instantaneous_external_power_W"), 1e-300))
    ax.plot(x[1:], ekin[1:], color="#1f4e79", lw=2.3, label="kinetik enerji")
    ax.set(xlabel="$\\log_{10}\\gamma$", ylabel="$\\log_{10}E_{kin}$ [J]", title="(c) 8 kg dış-besleme enerji faturası")
    ax2 = ax.twinx()
    ax2.plot(x[1:], power[1:], color="#c00000", lw=2.1, label="anlık dış güç")
    ax2.set_ylabel("$\\log_{10}P_{dış}$ [W]", color="#c00000")
    ax2.axhline(2, color="#c00000", ls=":", lw=1.2)
    ax2.text(33, 2.3, "100 W", ha="right", color="#c00000", fontsize=9)
    ax.grid(alpha=.25)

    ax = axes[1, 1]
    for run2, lab, col in [("1g_FULL", "1g tam gaz", "#d62728"), ("1g_RAMP", "1g rampa", "#2ca02c"), ("1g_PULSE50", "1g darbe %50", "#ff7f0e")]:
        ax.plot(earr(run2, "proper_years"), earr(run2, "log10_gamma"), label=lab, color=col, lw=2.2)
    ax.set(xlabel="Gemi öz zamanı [yıl]", ylabel="$\\log_{10}\\gamma$", title="(d) 1g gaz profilleri: serbest hız geçmişe bağlıdır")
    ax.legend(fontsize=8)
    ax.grid(alpha=.25)
    ax.text(.03, .06, "Son okumalar (85 öz-yıl):\nrampa γ≈10$^{37,58}$; darbe γ≈10$^{18,75}$", transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=.9))
    fig.tight_layout(rect=[0, 0, 1, .94])
    save(fig, "fig_exp41c_events_energy.png")


if __name__ == "__main__":
    fig_a()
    fig_b()
    fig_c()

