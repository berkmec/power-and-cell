"""Exp33: the SKY130 cell at 4 K - what the existing NIST 4 K raw data can still measure + lab registration.
Criteria fixed in results/exp33/predictions_exp33.json BEFORE this script computed any target observable.
Part A (existing data): C33-1 body effect VT(VB) vs Exp12 BSIM4 anchors; C33-2 body diode IB(VB);
C33-3 kink amplification A = dID_excess/|IB|; C33-4 single-step jump audit; C33-5 sweep-direction audit.
Part B (model-only): Exp17 split law with Exp32 measured constants -> registered lab table (C33-6/7);
the rD x lambda sub-grid is recomputed with law.w_split/lambda_c and compared to the stored Exp32 grid.
No physical lab measurement is performed. Stored E: files are read, never overwritten by this script;
output goes to results/exp33/cell33.json, falling back to ~/_cell33_out if the results folder is not writable."""
import json
import sys
import time
from pathlib import Path

import numpy as np

t0 = time.time()
_here = Path(__file__).resolve()
ROOT = _here.parents[2] if (_here.parents[2] / "results" / "exp32").exists() else Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp32_corrections"))
sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
import nist32  # noqa: E402  (read_mdm / idvd_metrics / idvg_metrics / size_of / slope, Exp32 definitions)
import law  # noqa: E402  (Exp17 split law)

DATA = ROOT / "results" / "exp32" / "nist_sky130"
VB_MODEL = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
DVT_MODEL = np.array([0.0, 0.044876592, 0.1025931, 0.183147286, 0.297640108, 0.315595025])


def monotonic(a):
    d = np.diff(np.asarray(a, dtype=float))
    return bool(np.all(d > 0) or np.all(d < 0))


def inventory_and_direction():
    files = {}
    nonmono = []
    for sub in ("4K", "300K"):
        for f in sorted((DATA / sub).glob("*.mdm")):
            bl = nist32.read_mdm(f)
            sweep_var = "VD" if "IDVD" in f.name else "VG"
            ok = True
            for var, c in bl:
                if not monotonic(c[sweep_var]):
                    ok = False
                    nonmono.append(f"{sub}/{f.name} block {var}")
            files[f"{sub}/{f.name}"] = dict(
                n_blocks=len(bl), heads=sorted(bl[0][1].keys()) if bl else [],
                VG=sorted({round(v["VG"], 4) for v, _ in bl if "VG" in v}),
                VB=sorted({round(v["VB"], 4) for v, _ in bl if "VB" in v}),
                VD=sorted({round(v["VD"], 4) for v, _ in bl if "VD" in v}),
                all_sweeps_monotonic=ok)
    return files, nonmono


def body_effect():
    per_dev = {}
    for f in sorted((DATA / "4K").glob("*IDVG*_4K.mdm")):
        W, L = nist32.size_of(f.name)
        met = nist32.idvg_metrics(nist32.read_mdm(f), W, L)
        tab = {}
        for k, v in met.items():
            vd_s, vb_s = k.split(",")
            tab.setdefault(float(vd_s.split("=")[1]), {})[float(vb_s.split("=")[1])] = v["VT"]
        per_dev[f.name] = dict(W=W, L=L, lvt="_lvt_" in f.name, VT_by_VD_VB={str(k): {str(a): b for a, b in sorted(d.items())} for k, d in tab.items()})
    # verdict at VD = 0.1, standard devices, largest |VB|
    rows = []
    for name, d in per_dev.items():
        t = d["VT_by_VD_VB"].get("0.1")
        if not t or "0.0" not in t:
            continue
        vbs = np.array(sorted(float(k) for k in t))
        vts = np.array([t[f"{v:g}"] if f"{v:g}" in t else t[str(v)] for v in vbs])
        vb_max = vbs[np.argmax(np.abs(vbs))]
        dvt = float(vts[np.argmax(np.abs(vbs))] - t["0.0"])
        gamma = float(np.polyfit(vbs, vts, 1)[0]) if len(vbs) > 1 else None
        model = float(np.interp(abs(vb_max), VB_MODEL, DVT_MODEL))
        rows.append(dict(file=name, lvt=d["lvt"], VB_values=[float(x) for x in vbs], VT_at_VB=[float(x) for x in vts],
                         VB_max=float(vb_max), dVT_meas=dvt, dVT_model_BSIM4=model,
                         ratio_abs=float(abs(dvt) / model) if model > 0 else None, gamma_dVT_dVB=gamma,
                         sign_ok=bool(dvt * vb_max < 0) if vb_max != 0 else None))
    std = [r for r in rows if not r["lvt"] and abs(r["VB_max"]) >= 0.3 and r["ratio_abs"]]
    if std:
        verdict = "pass" if all(0.5 <= r["ratio_abs"] <= 2.0 and r["sign_ok"] and 0.15 <= abs(r["gamma_dVT_dVB"]) <= 0.60 for r in std) else "fail"
    else:
        verdict = "unmeasurable-in-dataset"
    return dict(per_device=rows, verdict_C33_1=verdict,
                note="VT by nist32.idvg_metrics (1e-7*W/L criterion); model = Exp12 BSIM4 300 K anchors interpolated in |VB| (magnitude proxy, sign checked separately)")


def body_diode():
    rows = []
    has_ib = False
    for f in sorted((DATA / "4K").glob("*IDVG*_4K.mdm")):
        W, L = nist32.size_of(f.name)
        for var, c in nist32.read_mdm(f):
            if "IB" not in c or abs(var.get("VD", 9) - 0.1) > 1e-9:
                continue
            has_ib = True
            idd = np.abs(c["ID"])
            sel = idd < 1e-9
            ib = np.abs(c["IB"])[sel] if sel.any() else np.abs(c["IB"])
            rows.append(dict(file=f.name, W=W, lvt="_lvt_" in f.name, VB=float(var["VB"]),
                             IB_median_subthreshold_A=float(np.median(ib)) if len(ib) else None,
                             IB_median_per_um=float(np.median(ib) / W) if len(ib) else None))
    if not has_ib or not rows:
        return dict(rows=[], verdict_C33_2="unmeasurable-in-dataset", note="no IB column in IDVG blocks at VD=0.1")
    fwd = sorted({r["VB"] for r in rows if r["VB"] > 0})
    von = None
    for vb in fwd:
        meds = [r["IB_median_subthreshold_A"] for r in rows if r["VB"] == vb]
        if meds and max(meds) > 1e-9:
            von = vb
            break
    rev = [r for r in rows if r["VB"] == min(r["VB"] for r in rows)]
    rev_leak_max_per_um = max((r["IB_median_per_um"] for r in rev if r["IB_median_per_um"] is not None), default=None)
    if fwd:
        verdict = "pass" if (von is not None and 0.85 <= von <= 1.25) else ("fail (no 1 nA crossing up to max forward VB)" if von is None else "fail")
    else:
        verdict = "pass (reverse-only branch)" if (rev_leak_max_per_um is not None and rev_leak_max_per_um < 1e-11) else "fail (reverse-only branch)"
    return dict(rows=rows, forward_VB_values=fwd, Von_first_1nA=von, reverse_leak_max_per_um=rev_leak_max_per_um,
                verdict_C33_2=verdict, note="IB median over subthreshold VG points (ID<1 nA) at VD=0.1 V per (device,VB) block")


def kink_amplification():
    rows = []
    for f in sorted((DATA / "4K").glob("*IDVD*_4K.mdm")):
        bl = nist32.read_mdm(f)
        for vg in sorted({round(v["VG"], 3) for v, _ in bl if abs(v.get("VB", 0)) < 1e-9}):
            if vg < 0.85:
                continue
            m = nist32.idvd_metrics(bl, vg)
            vd, idd, ib = m["vd"], np.abs(m["idd"]), np.abs(m["ib"])
            g0 = nist32.slope(vd, idd, 0.5, 1.0)
            id10 = float(np.interp(1.0, vd, idd))
            id18 = float(np.interp(1.8, vd, idd))
            ib18 = float(np.interp(1.8, vd, ib))
            rows.append(dict(file=f.name, VG=vg, A_amplification=float((id18 - (id10 + g0 * 0.8)) / ib18) if ib18 > 0 else None,
                             IB_1p8_A=ib18, ID_1p8_A=id18))
    sel = [r for r in rows if abs(r["VG"] - 0.9) < 1e-9 and r["IB_1p8_A"] > 1e-9 and r["A_amplification"] is not None]
    med = float(np.median([r["A_amplification"] for r in sel])) if sel else None
    verdict = "pass" if (sel and 200 <= med <= 20000 and all(r["A_amplification"] > 30 for r in sel)) else ("fail" if sel else "unmeasurable-in-dataset")
    return dict(rows_VG09=[r for r in rows if abs(r["VG"] - 0.9) < 1e-9], n_in_verdict=len(sel), median_A_VG09=med, verdict_C33_3=verdict)


def jump_audit():
    per_temp = {}
    for sub, pat in (("4K", "*IDVD*_4K.mdm"), ("300K", "*IDVD*.mdm")):
        sweeps = []
        for f in sorted((DATA / sub).glob(pat)):
            for var, c in nist32.read_mdm(f):
                if abs(var.get("VB", 0)) > 1e-9:
                    continue
                vd, idd = np.asarray(c["VD"], dtype=float), np.abs(np.asarray(c["ID"], dtype=float))
                ok = (idd[:-1] > 1e-11) & (vd[:-1] >= 0.3) & (idd[1:] > 0)
                if not ok.any():
                    continue
                j = np.log10(idd[1:][ok] / idd[:-1][ok])
                sweeps.append(dict(file=f.name, VG=float(var.get("VG", float("nan"))), J_max_dex=float(np.max(j)),
                                   J_at_VD=float(vd[:-1][ok][int(np.argmax(j))])))
        per_temp[sub] = dict(n_sweeps=len(sweeps), max_J_dex=max((s["J_max_dex"] for s in sweeps), default=None),
                             n_over_0p5=sum(s["J_max_dex"] > 0.5 for s in sweeps),
                             n_in_0p3_0p5=sum(0.3 < s["J_max_dex"] <= 0.5 for s in sweeps),
                             worst=sorted(sweeps, key=lambda s: -s["J_max_dex"])[:5])
    ok = all(v["n_over_0p5"] == 0 and v["n_in_0p3_0p5"] <= 1 for v in per_temp.values())
    return dict(per_temperature=per_temp, verdict_C33_4="pass" if ok else "fail")


def max_rel_dev(a, b):
    devs = []
    def walk(x, y):
        if isinstance(y, dict) and isinstance(x, dict):
            for k in y:
                if k in x:
                    walk(x[k], y[k])
        elif isinstance(y, list) and isinstance(x, list):
            for u, v in zip(x, y):
                walk(u, v)
        elif isinstance(y, (int, float)) and isinstance(x, (int, float)):
            devs.append(abs(x - y) if y == 0 else abs(x - y) / abs(y))
    walk(a, b)
    return max(devs) if devs else None


def cell_law_lab_table():
    """Part B: Exp17 split law with the stored Exp32 measured constants (recompute grid; no root-finding redo)."""
    stored = json.loads((ROOT / "results" / "exp32" / "nist_sky130.json").read_text(encoding="utf-8"))
    c = stored["recal"]["corrected"]["measured_transfer"]
    f_per_A = 0.9 * np.sqrt(2 * 114.818 * 1.66053906660e-27 * 5e3 / 1.602176634e-19)
    grid, lab = {}, {}
    for r_D in (5.0, 10.0, 20.0):
        for lam in (40.0, 60.0, 120.0):
            ws = law.w_split(lam, r_D, c)
            ic = 0.75 * c["j_on"] * ws if np.isfinite(ws) else None
            key = f"rD{r_D:g}_lam{lam:g}"
            grid[key] = dict(lambda_c_um=law.lambda_c(r_D, c), w_split_um=ws, I_c_mA=ic * 1e3 if ic else None,
                             F_cell_uN_5kV_In=ic * f_per_A * 1e6 if ic else None)
            if ic:
                lab[key] = dict(W_min_um_for_4_cells=4 * ws,
                                I_steps_mA_k1_2_4_7=[k * ic * 1e3 for k in (1, 2, 4, 7)],
                                P_per_cell_mW_at_VM=ic * c["V_M"] * 1e3, P_per_cell_mW_at_1p8V=ic * 1.8 * 1e3,
                                P_7cells_mW_at_VM=7 * ic * c["V_M"] * 1e3)
    dev = max_rel_dev(grid, stored["cells"]["SKY130 4 K ölçülmüş girdiler"]) if "SKY130 4 K ölçülmüş girdiler" in stored.get("cells", {}) else \
        max_rel_dev(grid, list(stored["cells"].values())[1])
    return dict(constants_measured_transfer=c, recomputed_grid=grid, vs_stored_Exp32_cells_max_rel_dev=dev,
                consistency_C33_6="pass" if (dev is not None and dev < 1e-6) else "fail", lab_table=lab,
                design_point="rD10_lam60", note="V_D window for lab sweeps: [V_down, V_up] from constants; fronts at V_M")


if __name__ == "__main__":
    files, nonmono = inventory_and_direction()
    print(f"EXP33 inventory files={len(files)} nonmonotonic_sweeps={nonmono}")
    be = body_effect()
    print(f"EXP33 C33-1 body_effect verdict={be['verdict_C33_1']}")
    for r in be["per_device"]:
        print(f"  BE {r['file'][:52]} VB={r['VB_values']} dVT={r['dVT_meas']} model={r['dVT_model_BSIM4']} ratio={r['ratio_abs']} gamma={r['gamma_dVT_dVB']} sign_ok={r['sign_ok']}")
    bd = body_diode()
    print(f"EXP33 C33-2 body_diode verdict={bd['verdict_C33_2']} Von={bd.get('Von_first_1nA')} rev_leak_max_per_um={bd.get('reverse_leak_max_per_um')}")
    ka = kink_amplification()
    print(f"EXP33 C33-3 kink_amplification verdict={ka['verdict_C33_3']} median_A={ka['median_A_VG09']} n={ka['n_in_verdict']}")
    for r in ka["rows_VG09"]:
        print(f"  KA {r['file'][:52]} A={r['A_amplification']} IB18={r['IB_1p8_A']:.3e}")
    ja = jump_audit()
    print(f"EXP33 C33-4 jumps verdict={ja['verdict_C33_4']} " +
          " ".join(f"{k}:n={v['n_sweeps']},maxJ={v['max_J_dex']},>0.5={v['n_over_0p5']},0.3-0.5={v['n_in_0p3_0p5']}" for k, v in ja["per_temperature"].items()))
    cl = cell_law_lab_table()
    cm = cl["constants_measured_transfer"]
    print(f"EXP33 C33-6 consistency={cl['consistency_C33_6']} dev={cl['vs_stored_Exp32_cells_max_rel_dev']} "
          f"VM={cm['V_M']:.6f} Vdown={cm['V_down']:.6f} Vup={cm['V_up']:.6f} jon={cm['j_on']:.6e} dVc={cm['dV_c']:.6f}")
    dp = cl["lab_table"]["rD10_lam60"]
    print(f"EXP33 design rD10_lam60: Wmin4={dp['W_min_um_for_4_cells']:.1f}um Isteps_mA={[round(x,3) for x in dp['I_steps_mA_k1_2_4_7']]} "
          f"Pcell_VM={dp['P_per_cell_mW_at_VM']:.3f}mW P7_VM={dp['P_7cells_mW_at_VM']:.2f}mW")
    out = dict(inventory=files, C33_5_sweep_direction=dict(nonmonotonic=nonmono, verdict="confirmed (no hysteresis measurable)" if not nonmono else "refuted"),
               C33_1_body_effect=be, C33_2_body_diode=bd, C33_3_kink_amplification=ka, C33_4_jump_audit=ja,
               PartB_cell_law=cl,
               lab_status="NOT MEASURED this session: no cryostat / ballasted wide test chip / instruments; C33-6/7 are registered predictions + protocol inputs only")
    txt = json.dumps(out, indent=1, ensure_ascii=False, default=float)
    dst = ROOT / "results" / "exp33" / "cell33.json"
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(txt, encoding="utf-8")
        print("WROTE", dst)
    except Exception as ex:
        alt = Path.home() / "_cell33_out"
        alt.mkdir(parents=True, exist_ok=True)
        (alt / "cell33.json").write_text(txt, encoding="utf-8")
        print(f"WRITE_FAIL {dst} :: {type(ex).__name__}: {ex} -> WROTE {alt / 'cell33.json'}")
    print(f"ELAPSED_S={time.time()-t0:.1f}")
