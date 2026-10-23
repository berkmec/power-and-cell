import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
D = json.loads((ROOT / "results" / "exp40" / "h40.json").read_text(encoding="utf-8"))
OUT = ROOT / "results" / "exp40"
SP = D["species"]
LABELS = {
    "electron": ("Elektron", "#1f77b4"),
    "proton": ("Proton", "#2ca02c"),
    "indium_atom": ("İndiyum", "#d62728"),
    "planck_mass_particle": ("Planck kütlesi", "#9467bd"),
    "ship_8kg_single_quantum_counterfactual": ("8 kg tek-kuantum varsayımı", "#7f7f7f"),
}
LP = D["constants"]["planck_length_m"]
MP = D["constants"]["planck_mass_kg"]
EP_EV = D["constants"]["planck_energy_J"] / 1.602176634e-19


def arr(name, field):
    return np.array([r[field] for r in SP[name]["rows"]], dtype=float)


def save(fig, name):
    path = OUT / name
    try:
        fig.savefig(path, dpi=300)
    except PermissionError:
        path = Path.home() / name
        fig.savefig(path, dpi=300)
    plt.close(fig)
    print(path, path.stat().st_size)


def fig_a():
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle("Exp40A — Denklem 2 (`hν` hattı): frekans, enerji ve uzunluk ölçeği merdiveni", fontsize=16, fontweight="bold")
    names = ("electron", "proton", "indium_atom", "planck_mass_particle")
    ax = axes[0, 0]
    for name in names[:3]:
        lab, col = LABELS[name]
        ax.plot(arr(name, "n"), arr(name, "log10_lambda_h_m"), label=lab, color=col, lw=2.1)
    for y, txt in [(-10, "atom 10$^{-10}$ m"), (-14, "çekirdek 10$^{-14}$ m"), (-18, "elektrozayıf 10$^{-18}$ m"), (math.log10(LP), "Planck boyu")]:
        ax.axhline(y, color="black", ls=":", lw=1.1)
        ax.text(29.7, y + .35, txt, ha="right", fontsize=8)
    ax.set(xlabel="n  ($\\gamma=10^n$)", ylabel="$\\log_{10}\\lambda_h$ [m]", title="(a) $\\lambda_h=h/(\\gamma m_0c)$ ölçek inişi", ylim=(-37, -9))
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=.25)
    ax.text(.36, .06, "8 kg gemi tek kuantum sayılsa\n$\\lambda_h(n=0)=2,76\\times10^{-43}$ m: grafiğin altında", transform=ax.transAxes, fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", alpha=.88))

    ax = axes[0, 1]
    for name in names[:3]:
        lab, col = LABELS[name]
        ax.plot(arr(name, "n"), arr(name, "log10_frequency_Hz"), label=lab, color=col, lw=2.1)
    ax.axhline(math.log10(D["constants"]["c_m_s"] / LP), color="black", ls=":", lw=1.2)
    ax.text(29.5, math.log10(D["constants"]["c_m_s"] / LP) + .4, "c/$\\ell_P$", ha="right", fontsize=9)
    ax.set(xlabel="n", ylabel="$\\log_{10}\\nu$ [Hz]", title="(b) Frekans merdiveni $\\nu=\\gamma m_0c^2/h$")
    ax.legend(fontsize=8)
    ax.grid(alpha=.25)

    ax = axes[1, 0]
    for name in names[:3]:
        lab, col = LABELS[name]
        ax.plot(arr(name, "n"), arr(name, "log10_energy_eV"), label=lab, color=col, lw=2.1)
    ax.axhline(math.log10(EP_EV), color="black", ls=":", lw=1.2)
    ax.text(29.5, math.log10(EP_EV) + .5, "Planck enerjisi", ha="right", fontsize=9)
    ax.set(xlabel="n", ylabel="$\\log_{10}E$ [eV]", title="(c) Enerji merdiveni $E=h\\nu=\\gamma m_0c^2$")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=.25)

    ax = axes[1, 1]
    rows = [r for r in SP["indium_atom"]["rows"] if 0 < r["n"] <= 5]
    x = np.array([r["n"] for r in rows])
    beta = np.array([r["beta"] for r in rows])
    excess = np.array([100 * (r["photon_to_particle_momentum"] - 1) for r in rows])
    ax.plot(x, beta, color="#1f4e79", lw=2.3, label="$\\lambda_h/\\lambda_{dB}=\\beta$")
    ax.set(xlabel="n (0–5 yakın bakış)", ylabel="$\\lambda_h/\\lambda_{dB}=\\beta$", title="(d) hν boyu, parçacık de Broglie boyuyla aynı mı?", ylim=(0, 1.03))
    ax2 = ax.twinx()
    ax2.semilogy(x, excess, color="#c00000", lw=2.1, label="foton momentum fazlası [%]")
    ax2.set_ylabel("Foton momentum fazlası [%] (log)", color="#c00000")
    ax.scatter([1], [0.994987], color="#1f4e79", zorder=5)
    ax.annotate("n=1: β=0,99499\nfazla=%0,503", (1, .995), xytext=(2.0, .68), arrowprops=dict(arrowstyle="->"), fontsize=9)
    ax.grid(alpha=.25)
    fig.tight_layout(rect=[0, 0, 1, .94])
    save(fig, "fig_exp40a_scale.png")


def fig_b():
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle("Exp40B — Boyut denetimi: birim değişmez, ölçek Planck rejimlerine iner", fontsize=16, fontweight="bold")
    ax = axes[0, 0]
    logm = np.linspace(-31, 1, 320)
    nn = np.linspace(0, 30, 300)
    X, Y = np.meshgrid(logm, nn)
    log_lambda_ratio = math.log10(2 * math.pi * MP) - X - Y
    log_schwarz_ratio = math.log10(math.pi) + 2 * (math.log10(MP) - X - Y)
    log_energy_ratio = Y + X - math.log10(MP)
    im = ax.pcolormesh(X, Y, log_lambda_ratio, cmap="RdBu_r", vmin=-30, vmax=30, shading="auto")
    cs1 = ax.contour(X, Y, log_energy_ratio, levels=[0], colors=["#2ca02c"], linewidths=[1.8])
    cs2 = ax.contour(X, Y, log_schwarz_ratio, levels=[0], colors=["black"], linewidths=[1.8])
    cs3 = ax.contour(X, Y, log_lambda_ratio, levels=[0], colors=["white"], linewidths=[2.0])
    ax.clabel(cs1, fmt={0.0: "E=E_P"}, fontsize=8)
    ax.clabel(cs2, fmt={0.0: "λ=r_s"}, fontsize=8)
    ax.clabel(cs3, fmt={0.0: "λ=ℓ_P"}, fontsize=8)
    for name in ("electron", "proton", "indium_atom", "ship_8kg_single_quantum_counterfactual"):
        x = math.log10(SP[name]["mass_kg"])
        ax.axvline(x, color=LABELS[name][1], ls="--", lw=1.2)
        ax.text(x, 29.2, LABELS[name][0].split()[0], rotation=90, va="top", fontsize=7, color=LABELS[name][1])
    ax.set(xlabel="$\\log_{10}(m_0/[kg])$", ylabel="n", title="(a) Kütle–γ rejim haritası; renk: $\\log_{10}(\\lambda_h/\\ell_P)$")
    fig.colorbar(im, ax=ax, label="$\\log_{10}(\\lambda_h/\\ell_P)$")

    ax = axes[0, 1]
    for name in LABELS:
        lab, col = LABELS[name]
        ax.plot(arr(name, "n"), np.log10(arr(name, "lambda_over_planck")), label=lab, color=col, lw=2, ls="--" if "ship" in name else "-")
    ax.axhline(0, color="black", lw=1.3)
    ax.set(xlabel="n", ylabel="$\\log_{10}(\\lambda_h/\\ell_P)$", title="(b) Planck boyuna iniş", ylim=(-40, 26))
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=.25)

    ax = axes[1, 0]
    for name in LABELS:
        lab, col = LABELS[name]
        ax.plot(arr(name, "n"), np.log10(arr(name, "lambda_over_schwarzschild")), label=lab, color=col, lw=2, ls="--" if "ship" in name else "-")
    ax.axhline(0, color="black", lw=1.3)
    ax.set(xlabel="n", ylabel="$\\log_{10}(\\lambda_h/r_s(E))$", title="(c) Kendi enerjisinin Schwarzschild yarıçapıyla karşılaştırma", ylim=(-80, 50))
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=.25)

    ax = axes[1, 1]
    species = ("electron", "proton", "indium_atom")
    keys = [
        ("lambda_crosses_electroweak_1e-18_m", "10$^{-18}$ m", "#ff7f0e"),
        ("energy_reaches_planck_energy", "E=E_P", "#2ca02c"),
        ("lambda_equals_schwarzschild", "λ=r_s", "#1f77b4"),
        ("lambda_crosses_planck_length_m", "λ=ℓ_P", "#d62728"),
    ]
    x = np.arange(len(species))
    width = .19
    for j, (key, lab, col) in enumerate(keys):
        vals = [SP[s]["thresholds_n"][key] for s in species]
        bars = ax.bar(x + (j - 1.5) * width, vals, width=width, label=lab, color=col)
        ax.bar_label(bars, fmt="%.2f", fontsize=7, padding=1)
    ax.set_xticks(x, ["Elektron", "Proton", "İndiyum"])
    ax.set(ylabel="eşik n değeri", title="(d) Eşik sırası her türde aynı: E_P → r_s → ℓ_P", ylim=(0, 26))
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=.25)
    fig.tight_layout(rect=[0, 0, 1, .94])
    save(fig, "fig_exp40b_dimension_audit.png")


def fig_c():
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle("Exp40C — Varsayımsal ek boyut: denklemde R yok; erişim yalnız ek modelle tanımlanır", fontsize=16, fontweight="bold")
    curves = D["hypothetical_extra_dimension"]["curves"]
    ax = axes[0, 0]
    for name in ("electron", "proton", "indium_atom"):
        lab, col = LABELS[name]
        rows = curves[name]
        x = np.array([r["log10_radius_m"] for r in rows])
        ax.plot(x, [r["geometric_threshold_n"] for r in rows], color=col, ls="--", lw=1.8, label=f"{lab}: λ=R")
        ax.plot(x, [r["first_kk_energy_threshold_n"] for r in rows], color=col, lw=2.1, label=f"{lab}: KK E≥ħc/R")
    ax.axhline(30, color="black", ls=":", lw=1.3)
    ax.text(-10.5, 30.5, "nonilyon n=30", ha="right", fontsize=9)
    ax.set(xlabel="$\\log_{10}(R/[m])$", ylabel="eşik n", title="(a) Geometrik ve KK enerji eşikleri (0,798 dex kayık)", ylim=(-8, 38))
    ax.legend(fontsize=7, ncol=2, loc="lower left")
    ax.grid(alpha=.25)

    ax = axes[0, 1]
    rows = curves["indium_atom"]
    ax.plot([r["log10_radius_m"] for r in rows], [math.log10(r["first_kk_gap_eV"]) for r in rows], color="black", lw=2.3, label="KK aralığı ħc/R")
    base = math.log10(SP["indium_atom"]["rest_energy_eV"])
    for n in (1, 10, 20, 30):
        ax.axhline(base + n, ls="--", lw=1.3, label=f"In E (n={n})")
    ax.set(xlabel="$\\log_{10}(R/[m])$", ylabel="$\\log_{10}E$ [eV]", title="(b) İndiyum enerjisi ilk KK aralığını ne zaman aşar?")
    ax.legend(fontsize=8)
    ax.grid(alpha=.25)

    ax = axes[1, 0]
    scenarios = D["hypothetical_extra_dimension"]["scenarios"]
    for key, vals_by_species in scenarios.items():
        rows = [r for r in vals_by_species["indium_atom"] if r["n"] >= .1 and r["kmax_continuous"] > 0]
        ax.plot([r["n"] for r in rows], [math.log10(r["kmax_continuous"]) for r in rows], lw=2, label=key)
    ax.axhline(0, color="black", lw=1.4)
    ax.text(29.5, .5, "ilk KK modu", ha="right", fontsize=9)
    ax.set(xlabel="n", ylabel="$\\log_{10}k_{max}$ (sürekli kestirim)", title="(c) İndiyum için kinematik KK mod kapasitesi")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=.25)

    ax = axes[1, 1]
    m_in = SP["indium_atom"]["mass_kg"]
    nn = np.linspace(0, 30, 301)
    logr = np.linspace(-35, -10, 260)
    X, Y = np.meshgrid(nn, logr)
    beta = np.sqrt(np.maximum(0.0, 1.0 - 10.0 ** (-2 * X)))
    logk = np.full_like(X, -80.0)
    positive = beta > 0
    logk[positive] = X[positive] + np.log10(beta[positive]) + math.log10(m_in * D["constants"]["c_m_s"] / D["constants"]["hbar_J_s"]) + Y[positive]
    im = ax.pcolormesh(X, Y, logk, cmap="viridis", vmin=-15, vmax=15, shading="auto")
    cs = ax.contour(X, Y, logk, levels=[0], colors=["white"], linewidths=[2])
    ax.clabel(cs, fmt={0.0: "kmax=1"}, fontsize=8)
    lambda_c = SP["indium_atom"]["compton_wavelength_h_over_mc_m"]
    ax.plot(nn, np.log10(lambda_c) - nn, color="#ff7f0e", lw=2, label="λ_h=R geometrik çizgisi")
    ax.axhline(math.log10(LP), color="white", ls=":", lw=1.4)
    ax.text(29.5, math.log10(LP) + .6, "R=Planck boyu", ha="right", color="white", fontsize=8)
    ax.set(xlabel="n", ylabel="$\\log_{10}(R/[m])$", title="(d) İndiyum: R–n erişim haritası; renk $\\log_{10}k_{max}$", ylim=(-35, -10))
    ax.legend(fontsize=8, loc="upper right")
    fig.colorbar(im, ax=ax, label="$\\log_{10}k_{max}$")
    fig.tight_layout(rect=[0, 0, 1, .94])
    save(fig, "fig_exp40c_extra_dimension.png")


if __name__ == "__main__":
    fig_a()
    fig_b()
    fig_c()

