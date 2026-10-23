"""Exp32-C4: SKY130's OWN cryogenic measurement against the cell model (Exp12/17).
Data (copied to results/exp32/nist_sky130/): NIST mds2-2997 'CMOS Sky130 Primitives measured at cryogenic temperatures'
(ETest tile, Lakeshore closed-cycle probe station, chuck and surroundings at 4 K, Agilent 4156B; Akturk, Hoskins, Shrestha et al.)
and the room-temperature measurement of the same structures (google/skywater-pdk-sky130-raw-data).
Definitions and criteria fixed in results/exp32/predictions_exp32_nist.json before this script.
usage: python nist32.py [fig] -> results/exp32/nist_sky130.json, fig_exp32_nist.png"""
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp32"
DATA = OUT / "nist_sky130"


def read_mdm(path):
    t = Path(path).read_text(encoding="latin-1")
    blocks = []
    for b in t.split("BEGIN_DB")[1:]:
        var = {k: float(v) for k, v in re.findall(r"ICCAP_VAR (\w+)\s+(\S+)", b)}
        body = b.split("END_DB")[0].splitlines()
        head = next(l for l in body if l.strip().startswith("#")).replace("#", "").split()
        rows = np.array([[float(x) for x in l.split()] for l in body if re.match(r"^\s*-?[0-9.]", l)])
        blocks.append((var, {h: rows[:, i] for i, h in enumerate(head)}))
    return blocks


def slope(x, y, lo, hi):
    s = (x >= lo - 1e-9) & (x <= hi + 1e-9)
    return float(np.polyfit(x[s], y[s], 1)[0])


def idvd_metrics(blocks, vg, vb=0.0):
    b = next(c for v, c in blocks if abs(v.get("VG", -9) - vg) < 1e-6 and abs(v.get("VB", -9) - vb) < 1e-6)
    vd, idd = b["VD"], b["ID"]
    ib = b.get("IB", np.zeros_like(vd))
    g0 = slope(vd, idd, 0.5, 1.0)
    g1 = slope(vd, idd, 1.6, 1.8)
    g = np.gradient(idd, vd)
    onset = None
    for i in range(len(vd) - 2):
        if vd[i] >= 0.5 and all(g[i + k] > 1.5 * g0 for k in range(3)):
            onset = float(vd[i]); break
    big = np.abs(ib) > 1e-9
    ib_on = float(vd[np.argmax(big)]) if big.any() else None
    M = np.where(big, np.abs(ib) / np.maximum(idd, 1e-30), np.nan)
    return dict(VG=vg, kink_ratio=g1 / g0 if g0 > 0 else None, kink_onset_V=onset, IB_1nA_onset_V=ib_on,
                ID_at_1p2=float(np.interp(1.2, vd, idd)), ID_at_1p8=float(np.interp(1.8, vd, idd)),
                M_at={f"{v:.1f}": float(np.interp(v, vd, np.nan_to_num(M, nan=0))) for v in (1.4, 1.6, 1.8)}, vd=vd, idd=idd, ib=ib)


def idvg_metrics(blocks, W, L):
    out = {}
    for var, c in blocks:
        vg, idd = c["VG"], np.abs(c["ID"])
        lg = np.log10(np.clip(idd, 1e-13, None))
        i0 = np.argmax(lg > -11)
        vt = float(np.interp(np.log10(1e-7 * W / L), lg[i0:], vg[i0:]))
        ss = np.diff(vg) / np.maximum(np.diff(lg), 1e-9) * 1e3
        ok = (idd[1:] > 1e-10) & (idd[1:] < 1e-6)
        out[f"VD={var['VD']:g},VB={var['VB']:g}"] = dict(VT=vt, SSmin_mV_dec=float(np.min(ss[ok])) if ok.any() else None)
    return out


def size_of(name):
    m = re.search(r"_w([0-9p]+)u_l([0-9p]+)u", name)
    return float(m.group(1).replace("p", ".")), float(m.group(2).replace("p", "."))


def survey():
    rows = []
    for f in sorted((DATA / "4K").glob("*IDVD*_4K.mdm")):
        W, L = size_of(f.name)
        lvt = "_lvt_" in f.name
        r = dict(file=f.name, W=W, L=L, lvt=lvt)
        bl = read_mdm(f)
        vgs = sorted({round(v["VG"], 3) for v, _ in bl})
        for vg in vgs:
            if vg < 0.85:
                continue
            m = idvd_metrics(bl, vg)
            r[f"VG{vg:g}"] = {k: v for k, v in m.items() if k not in ("vd", "idd", "ib")}
        g = f.with_name(f.name.replace("IDVD", "IDVG"))
        if g.exists():
            r["IDVG"] = idvg_metrics(read_mdm(g), W, L)
        rows.append(r)
    return rows


def room_temperature():
    out = {}
    for f in sorted((DATA / "300K").glob("*IDVD*.mdm")):
        bl = read_mdm(f)
        out[f.name] = {f"VG{vg:g}": {k: v for k, v in idvd_metrics(bl, vg).items() if k not in ("vd", "idd", "ib")} for vg in (0.72, 1.08)}
        g = f.with_name(f.name.replace("IDVD", "IDVG"))
        if g.exists():
            W, L = size_of(f.name)
            out[f.name]["IDVG"] = idvg_metrics(read_mdm(g), W, L)
    return out


def recalibrate():
    """Exp17 constants with SKY130-measured 4 K multiplication (W7/L0.15, V_G = 0.9 V) instead of BSIM4-300K x fii."""
    bl = read_mdm(DATA / "4K" / "sky130_fd_pr__nfet_01v8_w7u_l0p15u_m1(8392_7_8_IDVD)_4K.mdm")
    m = idvd_metrics(bl, 0.9)
    vd, idd, ib = m["vd"], m["idd"], m["ib"]
    s = np.abs(ib) > 1e-9
    a1, a0 = np.polyfit(vd[s], np.log10(np.abs(ib[s]) / idd[s]), 1)
    sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
    import law
    import reduced_model as rm
    bsim_M = rm.M_of
    base = law.constants()                                     # Exp17 as published (BSIM4 300 K M, fii = 10, ss = 20 mV/dec)
    rm.M_of = lambda v: np.minimum(10 ** (a0 + a1 * np.asarray(v, dtype=float)), 1.0)
    res = dict(fit=dict(log10M_intercept=float(a0), slope_dec_per_V=float(a1), n_points=int(s.sum()),
                        M_meas_at={f"{v}": float(10 ** (a0 + a1 * v)) for v in (0.8, 1.2, 1.4, 1.6, 1.8)},
                        M_bsim_x10_at={f"{v}": float(10 * bsim_M(v)) for v in (0.8, 1.2, 1.4, 1.6, 1.8)}),
               exp17_published=base, corrected={})
    for ss in (0.006, 0.010, 0.020):
        p = dict(law.BASE, fii=1, ss=ss)
        try:
            res["corrected"][f"ss={ss}"] = law.constants(p)
        except Exception as e:                                 # no bistability inside the search window
            res["corrected"][f"ss={ss}"] = dict(error=repr(e))
    # measured 4 K transfer curve (W7/L0.15, V_D = 1.8 V, V_B = 0): j_d(V_ov) = I_D / W, V_ov = V_G - V_T(V_D = 1.8 V)
    g = read_mdm(DATA / "4K" / "sky130_fd_pr__nfet_01v8_w7u_l0p15u_m1(8392_7_8_IDVG)_4K.mdm")
    c = next(c for v, c in g if abs(v["VD"] - 1.8) < 1e-6 and abs(v["VB"]) < 1e-6)
    vt18 = idvg_metrics([(dict(VD=1.8, VB=0.0), c)], 7.0, 0.15)["VD=1.8,VB=0"]["VT"]
    vov_m, jd_m = c["VG"] - vt18, np.clip(np.abs(c["ID"]) / 7.0, 1e-30, None)
    jd_meas = lambda vov: 10 ** np.interp(vov, vov_m, np.log10(jd_m))
    jd_orig = law._jd
    law._jd = lambda V, p: jd_meas(p["vov0"] + rm.dvt_of(V))
    try:
        res["corrected"]["measured_transfer"] = law.constants(dict(law.BASE, fii=1))
    except Exception as e:
        res["corrected"]["measured_transfer"] = dict(error=repr(e))
    res["transfer"] = dict(VT_VD1p8=vt18, jd_at_vov={f"{v}": float(jd_meas(v)) for v in (-0.1, 0.0, 0.1, 0.2, 0.3)})
    law._jd = jd_orig
    rm.M_of = bsim_M
    return res


def fig(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
    a = ax[0]
    b4 = read_mdm(DATA / "4K" / "sky130_fd_pr__nfet_01v8_w7u_l0p15u_m1(8392_7_8_IDVD)_4K.mdm")
    b3 = read_mdm(DATA / "300K" / "sky130_fd_pr__nfet_01v8_w7u_l0p15u_m1(8392_7_8_IDVD).mdm")
    for bl, vg, lab, ls in ((b4, 0.9, "4 K, V_G = 0,9 V", "-"), (b4, 1.2, "4 K, V_G = 1,2 V", "-"), (b3, 0.72, "300 K, V_G = 0,72 V", "--"),
                            (b3, 1.08, "300 K, V_G = 1,08 V", "--")):
        m = idvd_metrics(bl, vg)
        a.plot(m["vd"], m["idd"] / np.interp(1.0, m["vd"], m["idd"]), ls, label=lab)
    a.set(xlabel="V_D (V)", ylabel="I_D / I_D(1 V)", title="(a) SKY130 nfet_01v8 W7/L0,15: 4 K'de kink, 300 K'de yok")
    a.legend(fontsize=8)
    a = ax[1]
    m = idvd_metrics(b4, 0.9)
    s = np.abs(m["ib"]) > 1e-9
    a.semilogy(m["vd"][s], np.abs(m["ib"][s]) / m["idd"][s], "o", label="ölçülen 4 K |I_B|/I_D (V_G = 0,9 V)")
    v = np.linspace(0.6, 1.9, 100)
    f = res["recal"]["fit"]
    a.semilogy(v, 10 ** (f["log10M_intercept"] + f["slope_dec_per_V"] * v), "k-", lw=0.8, label="log-doğrusal uydurma")
    sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))
    import reduced_model as rm
    a.semilogy(v, 10 * rm.M_of(v), "r--", label="Exp12/17: BSIM4 300 K × F_ii = 10")
    a.axvline(0.8018, color="gray", ls=":", label="Exp17 V_M = 0,80 V")
    a.set(xlabel="V_D (V)", ylabel="çoğaltma M", title="(b) darbe iyonlaşması: ölçülen vs modelin varsayımı")
    a.legend(fontsize=7)
    a = ax[2]
    rows = [(k, d.get("V_M"), d.get("V_down")) for k, d in res["recal"]["corrected"].items() if isinstance(d, dict) and "V_M" in d]
    labels = ["Exp17 (yayımlanan)"] + [f"ölçülen M, {k}" for k, _, _ in rows]
    vm = [res["recal"]["exp17_published"]["V_M"]] + [x[1] for x in rows]
    a.barh(range(len(vm)), vm, color=["C3"] + ["C0"] * len(rows))
    a.set(yticks=range(len(vm)), yticklabels=labels, xlabel="Maxwell (kenet) gerilimi V_M (V)", title="(c) hücre kenedi: SKY130 4 K girdileriyle")
    fig.tight_layout()
    fig.savefig(OUT / "fig_exp32_nist.png", dpi=125)
    print("saved", OUT / "fig_exp32_nist.png")


if __name__ == "__main__":
    res = dict(survey_4K=survey(), room_T=room_temperature(), recal=recalibrate())
    (OUT / "nist_sky130.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float), encoding="utf-8")
    for r in res["survey_4K"]:
        k9 = r.get("VG0.9", {})
        vt = r.get("IDVG", {}).get("VD=0.1,VB=0", {})
        print(f"{'lvt ' if r['lvt'] else ''}W{r['W']:g}/L{r['L']:g}: VT {vt.get('VT', float('nan')):.3f} SS {vt.get('SSmin_mV_dec')}  "
              f"VG0.9 kink {k9.get('kink_ratio')} onset {k9.get('kink_onset_V')} IB-onset {k9.get('IB_1nA_onset_V')} M1.8 {k9.get('M_at', {}).get('1.8')}")
    for f, d in res["room_T"].items():
        print("300K", f[22:50], {k: (round(v["kink_ratio"], 2) if isinstance(v, dict) and v.get("kink_ratio") else v) for k, v in d.items() if k != "IDVG"},
              {k: round(v["VT"], 3) for k, v in d.get("IDVG", {}).items()})
    print("recal fit", res["recal"]["fit"])
    print("Exp17 published", res["recal"]["exp17_published"])
    for k, v in res["recal"]["corrected"].items():
        print("corrected", k, v)
    if "fig" in sys.argv:
        fig(res)


def cells(res):
    """Exp17 split law with the published and the SKY130-measured-input constants: lambda_c, w_split, I_c = 0.75 j_on w_s,
    and the In+ thrust per cell at 5 kV (Exp28/29 equation, gamma = 0.9)."""
    sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
    import law
    sets = {"Exp17 (yayımlanan)": res["recal"]["exp17_published"], "SKY130 4 K ölçülmüş girdiler": res["recal"]["corrected"]["measured_transfer"]}
    f_per_A = 0.9 * np.sqrt(2 * 114.818 * 1.66053906660e-27 * 5e3 / 1.602176634e-19)
    out = {}
    for name, c in sets.items():
        rows = {}
        for r_D in (5.0, 10.0, 20.0):
            for lam in (40.0, 60.0, 120.0):
                ws = law.w_split(lam, r_D, c)
                ic = 0.75 * c["j_on"] * ws if np.isfinite(ws) else None
                rows[f"rD{r_D:g}_lam{lam:g}"] = dict(lambda_c_um=law.lambda_c(r_D, c), w_split_um=ws, I_c_mA=ic * 1e3 if ic else None,
                                                    F_cell_uN_5kV_In=ic * f_per_A * 1e6 if ic else None)
        out[name] = rows
    return out


if __name__ == "__main__" and "cells" in sys.argv:
    res = json.loads((OUT / "nist_sky130.json").read_text(encoding="utf-8"))
    res["cells"] = cells(res)
    (OUT / "nist_sky130.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float), encoding="utf-8")
    for name, rows in res["cells"].items():
        for k, v in rows.items():
            print(name, k, {a: (round(b, 3) if isinstance(b, float) else b) for a, b in v.items()})
