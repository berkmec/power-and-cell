"""Exp22b: is the quantum confinement of the inversion layer visible in the foundry's own SKY130 data, and how large is
the oxide tunnelling current (both room-temperature quantum effects)?
1. ngspice, official nfet_01v8 tt card (capMod = 2: BSIM4 charge-thickness model), 27 C: gate capacitance C_gg(V_G) of
   W = 10 um devices with L = 25 and 50 um (same bin); the difference removes overlap and fringe capacitance.
2. Same stack (t_ox 4.148 nm, N_A 1.7e17, n+ poly) self-consistent Poisson-Schroedinger (src/qcell/mos1d.py, 300 K) and the
   classical (Boltzmann, no quantisation) analytic C-V. The curves are aligned at the inflection of the inversion edge.
3. Direct tunnelling through 4.148 nm SiO2 (Schuegraf-Hu WKB form, m_ox = 0.4-0.5 m0, barrier 3.1 eV) vs thermionic emission.
usage: python cv22.py  ->  results/exp22/cv_tunnel.json"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))
import spice_tools as st  # noqa: E402
from qcell.mos1d import Mos, Q, KB, HBAR, M0  # noqa: E402

OUT = ROOT / "results" / "exp22"
EPS0, EPS_SI, EPS_OX = 8.8541878128e-12, 11.7, 3.9
TOX, NA, NI = 4.148e-9, 1.7e23, 1.0e16
COX = EPS_OX * EPS0 / TOX                                     # F/m^2
VG = np.round(np.arange(-1.0, 2.0001, 0.05), 3)
LS, WD, F_AC = (25.0, 50.0), 10.0, 1e5


def spice_cv():
    lines = [st.header(1.0, 1.0)]
    for L in LS:
        for k, v in enumerate(VG):
            lines.append(f"Vg{k}_{L:g} g{k}_{L:g} 0 dc {v} ac 1\nX{k}_{L:g} 0 g{k}_{L:g} 0 0 sky130_fd_pr__nfet_01v8 "
                         f"l={L:g} w={WD:g}\n")
    lines.append(f".ac lin 1 {F_AC:g} {F_AC:g}\n.control\nrun\n")
    names = []
    for L in LS:
        for k in range(len(VG)):
            lines.append(f"let c{k}_{L:g} = imag(i(vg{k}_{L:g}))\n")
            names.append(f"c{k}_{L:g}")
    lines.append(f"wrdata {st.wsl(st.SP / 'cv22.dat')} " + " ".join(names) + "\n.endc\n.end\n")
    log, bad = st.run("cv22", "* exp22 gate C-V, official SKY130 nfet_01v8 tt card\n" + "".join(lines))
    if bad:
        print("ngspice:", bad[:5])
    row = np.atleast_2d(st.load("cv22"))[0]
    vals = np.abs(row[1::2]) / (2 * np.pi * F_AC)                 # F per device
    c25, c50 = vals[:len(VG)], vals[len(VG):]
    c_int = (c50 - c25) / ((LS[1] - LS[0]) * WD * 1e-12)          # F/m^2
    return dict(VG=VG.tolist(), C_L25_F=c25.tolist(), C_L50_F=c50.tolist(), C_area_F_m2=c_int.tolist())


def classical_cv():
    vt = KB * 300.0 / Q
    phif = vt * np.log(NA / NI)
    vfb = -(1.12 / 2 + phif)
    psi = np.linspace(-0.3, 1.3, 3201)
    b = psi / vt
    F2 = (np.exp(-b) + b - 1) + (NI / NA) ** 2 * (np.exp(b) - b - 1)
    Qs = np.sign(psi) * np.sqrt(2 * EPS_SI * EPS0 * KB * 300.0 * NA * np.maximum(F2, 0))
    vg = vfb + psi + Qs / COX
    return vg, Qs


def quantum_cv():
    mos = Mos(T=300.0)
    rows = []
    for ps in np.round(np.arange(0.30, 1.301, 0.02), 3):
        s = mos.solve(float(ps))
        rows.append((s["Vg"], EPS_SI * EPS0 * s["Es"], s["ns_cm2"]))
        print(f"  P-S psi_s {ps:.2f}: Vg {s['Vg']:.3f} ns {s['ns_cm2']:.2e}", flush=True)
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2]


def deriv(x, y):
    o = np.argsort(x)
    x, y = x[o], y[o]
    return 0.5 * (x[1:] + x[:-1]), np.diff(y) / np.diff(x)


def align(vg, c):
    """inflection of the inversion edge (largest dC/dV_G for V_G above the C minimum)"""
    i0 = int(np.argmin(c))
    xm, dc = deriv(vg[i0:], c[i0:])
    return float(xm[np.argmax(dc)])


def tunnelling():
    out = {}
    for m_ox in (0.4, 0.5):
        phi = 3.1 * Q
        A = Q**3 / (16 * np.pi**2 * HBAR * phi)
        B = 4 * np.sqrt(2 * m_ox * M0) * phi**1.5 / (3 * HBAR * Q)
        for vox in (0.5, 1.0, 1.5, 1.8):
            E = vox / TOX
            r = 1 - vox / 3.1
            J = A * E**2 / (1 - np.sqrt(r)) ** 2 * np.exp(-B * (1 - r**1.5) / E)   # A/m^2
            out[f"m{m_ox}_V{vox}"] = J * 1e-4                                        # A/cm^2
    out["thermionic_A_cm2_300K"] = 120.0 * 300.0**2 * np.exp(-3.1 / (KB * 300.0 / Q))
    out["area_100x100um_cm2"] = 1e-4
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sp = spice_cv()
    vgs = np.array(sp["VG"])
    cs = np.array(sp["C_area_F_m2"])
    vgc, qc = classical_cv()
    xc, ccl = deriv(vgc, qc)
    vgq, qq, nsq = quantum_cv()
    xq, cqu = deriv(vgq, qq)
    res = dict(spice=sp, classical=dict(VG=xc[::20].tolist(), C=ccl[::20].tolist()),
               quantum=dict(VG=xq.tolist(), C=cqu.tolist(), ns_cm2=nsq.tolist()), C_ox=COX, tunnelling=tunnelling())
    comp = {}
    for name, x, c in (("spice_card", vgs, cs), ("classical", xc, ccl), ("quantum_PS", xq, cqu)):
        v0 = align(x, c)
        comp[name] = dict(V_inflection=v0, **{f"C/Cox at +{d:.1f} V": float(np.interp(v0 + d, x, c) / COX)
                                              for d in (0.3, 0.6, 0.9)})
    res["comparison"] = comp
    (OUT / "cv_tunnel.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print("C_ox = %.4f uF/cm2" % (COX * 100))
    for k, v in comp.items():
        print(f"{k:12s}", {kk: round(vv, 4) for kk, vv in v.items()})
    print("tunnelling (A/cm2):", {k: f"{v:.2e}" for k, v in res["tunnelling"].items()})


def fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    qs = json.loads((OUT / "quantum_signatures.json").read_text(encoding="utf-8"))
    cv = json.loads((OUT / "cv_tunnel.json").read_text(encoding="utf-8"))
    f, axs = plt.subplots(1, 3, figsize=(16, 5))
    ax = axs[0]
    ns = np.array([r["ns_cm2"] for r in qs["rows"]])
    ax.plot(ns, [r["dE10_meV"] for r in qs["rows"]], "o-", label="E₁ − E₀ (ikinci alt bant)")
    ax.plot(ns, [r["EF_minus_E0_4K_meV"] for r in qs["rows"]], "s-", label="E_F − E₀ (4 K)")
    ax.plot(ns, [100 * r["ground_fraction_300K"] for r in qs["rows"]], "^--", c="C3", label="taban alt bantta % (300 K)")
    ax.set_xscale("log")
    ax.set_xlabel("n_s (cm⁻²)")
    ax.set_ylabel("meV  /  %")
    ax.set_title("(a) inversiyon katmanı alt bantları: 4 K'de tüm aralıkta\ntek alt bant (elektrik kuantum limiti)", fontsize=10)
    ax.legend(fontsize=8)
    ax = axs[1]
    vgs, cs = np.array(cv["spice"]["VG"]), np.array(cv["spice"]["C_area_F_m2"])
    for (x, c), lab, key in (((vgs, cs), "SKY130 resmî kart (BSIM4 capMod 2)", "spice_card"),
                             ((np.array(cv["classical"]["VG"]), np.array(cv["classical"]["C"])), "klasik (kuantumsuz)", "classical"),
                             ((np.array(cv["quantum"]["VG"]), np.array(cv["quantum"]["C"])), "kuantum Poisson–Schrödinger", "quantum_PS")):
        v0 = cv["comparison"][key]["V_inflection"]
        ax.plot(x - v0, c / cv["C_ox"], label=lab)
    ax.set_xlim(-0.6, 1.2)
    ax.set_ylim(0.0, 1.05)
    ax.axhline(1.0, c="k", lw=0.8, ls=":")
    ax.set_xlabel("V_G − V_bükülme (V)")
    ax.set_ylabel("C / C_ox")
    ax.set_title("(b) inversiyon kapasitesi: kart klasiğin %2,4 altında,\nkuantum P-S'nin %5 üstünde (+0,9 V'ta)", fontsize=10)
    ax.legend(fontsize=8)
    ax = axs[2]
    t = cv["tunnelling"]
    for m, col in (("0.4", "C0"), ("0.5", "C1")):
        v = [0.5, 1.0, 1.5, 1.8]
        ax.semilogy(v, [t[f"m{m}_V{x}"] * 1e-2 for x in v], "o-", c=col, label=f"doğrudan tünelleme, 1 mm², m_ox = {m} m₀")
    ax.axhline(1e-15, c="k", ls="--", lw=1, label="elektrometre tabanı ~1 fA")
    ax.axhline(t["thermionic_A_cm2_300K"] * 1e-2, c="C3", ls=":", lw=1, label="ısıl emisyon (300 K): ~1e-47 A (eksen dışında)")
    ax.set_ylim(1e-20, 1e-9)
    ax.set_xlabel("oksit gerilimi (V)")
    ax.set_ylabel("kapı akımı (A), 1 mm² kapasitör")
    ax.set_title("(c) 4,148 nm SiO₂'den tünelleme: oda sıcaklığında\nölçülebilir kuantum etkisi (fA–pA, 1,8 V)", fontsize=10)
    ax.legend(fontsize=7)
    f.suptitle("Exp22: SKY130 silikonunda kuantum durumları — veriden ve modelden ne görülüyor, nasıl ölçülür", fontsize=12)
    f.tight_layout(rect=(0, 0, 1, 0.94))
    f.savefig(OUT / "fig_exp22.png", dpi=95)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "fig":
        fig()
    else:
        main()
