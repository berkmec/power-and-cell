"""Exp33 POST-REGISTERED EXPLORATORY checks (not part of the C33 verdicts in predictions_exp33.json):
1) diagnostics of the two worst jump-audit sweeps (actual ID values at the flagged step);
2) sensitivity of the Exp17/Exp32 cell constants to the C33-1 finding: the 4 K measured body effect
   (reverse VB, |dVT/dVB| std median 0.105; |dVT| ratio vs clamped BSIM4 anchor at |VB|=1.5 median 0.498)
   is weaker than the BSIM4 300 K anchor (anchor mean slope 0.316 to VB=1). Recompute law.constants on the
   exact nist32 measured_transfer code path with dvt_of scaled by s = 1.0 (baseline, must match stored
   Exp32 constants), 0.50 and 0.33. Reverse-bias magnitude is used as a forward scale proxy - the
   forward/reverse asymmetry is unknown, so this is a sensitivity bound, not a corrected constant set."""
import json
import sys
import time
from pathlib import Path

import numpy as np

t0 = time.time()
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp32_corrections"))
sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
import nist32
import law

DATA = ROOT / "results" / "exp32" / "nist_sky130"
out = dict(kind="POST-REGISTERED EXPLORATORY (not a C33 verdict)")

# ---- 1) worst-jump diagnostics ----
diag = []
for sub, fname, vg_t in (("300K", "sky130_fd_pr__nfet_01v8_w7u_l8u_m1(8392_5_6_IDVD).mdm", 0.0),
                         ("4K", "sky130_fd_pr__nfet_01v8_lvt_w0p42u_l0p15u_m1(8391_9_10_IDVD)_4K.mdm", 0.6)):
    for var, c in nist32.read_mdm(DATA / sub / fname):
        if abs(var.get("VB", 0)) > 1e-9 or abs(var.get("VG", -9) - vg_t) > 1e-9:
            continue
        vd, idd = np.asarray(c["VD"], float), np.abs(np.asarray(c["ID"], float))
        ok = (idd[:-1] > 1e-11) & (vd[:-1] >= 0.3) & (idd[1:] > 0)
        idx = np.nonzero(ok)[0]
        j = np.log10(idd[1:][ok] / idd[:-1][ok])
        i = int(idx[int(np.argmax(j))])
        diag.append(dict(temp=sub, file=fname, VG=vg_t, VD_step=[float(vd[i]), float(vd[i + 1])],
                         ID_step_A=[float(idd[i]), float(idd[i + 1])], J_dex=float(np.max(j)),
                         ID_median_sweep_A=float(np.median(idd))))
        print(f"EXP33DIAG {sub} VG={vg_t} VD {vd[i]:.3f}->{vd[i+1]:.3f} ID {idd[i]:.3e}->{idd[i+1]:.3e} A J={np.max(j):.3f} dex")
out["worst_jump_diagnostics"] = diag

# ---- 2) dvt sensitivity on the nist32 measured_transfer code path ----
bl = nist32.read_mdm(DATA / "4K" / "sky130_fd_pr__nfet_01v8_w7u_l0p15u_m1(8392_7_8_IDVD)_4K.mdm")
m = nist32.idvd_metrics(bl, 0.9)
vd_a, idd_a, ib_a = m["vd"], m["idd"], m["ib"]
s_mask = np.abs(ib_a) > 1e-9
a1, a0 = np.polyfit(vd_a[s_mask], np.log10(np.abs(ib_a[s_mask]) / idd_a[s_mask]), 1)
g = nist32.read_mdm(DATA / "4K" / "sky130_fd_pr__nfet_01v8_w7u_l0p15u_m1(8392_7_8_IDVG)_4K.mdm")
c18 = next(c for v, c in g if abs(v["VD"] - 1.8) < 1e-6 and abs(v["VB"]) < 1e-6)
vt18 = nist32.idvg_metrics([(dict(VD=1.8, VB=0.0), c18)], 7.0, 0.15)["VD=1.8,VB=0"]["VT"]
vov_m, jd_m = c18["VG"] - vt18, np.clip(np.abs(c18["ID"]) / 7.0, 1e-30, None)
jd_meas = lambda vov: 10 ** np.interp(vov, vov_m, np.log10(jd_m))

rm = law.rm
orig_M, orig_dvt, orig_jd = rm.M_of, rm.dvt_of, law._jd
rm.M_of = lambda v: np.minimum(10 ** (a0 + a1 * np.asarray(v, dtype=float)), 1.0)
sens = {}
for s in (1.0, 0.50, 0.33):
    law._jd = lambda V, p, s=s: jd_meas(p["vov0"] + s * orig_dvt(V))
    try:
        cst = law.constants(dict(law.BASE, fii=1))
        ws = law.w_split(60.0, 10.0, cst)
        sens[f"dvt_scale_{s:g}"] = dict(constants=cst, ws_rD10_lam60_um=ws,
                                        Ic_rD10_lam60_mA=0.75 * cst["j_on"] * ws * 1e3 if np.isfinite(ws) else None)
        print(f"EXP33SENS s={s:g} VM={cst['V_M']:.6f} Vdown={cst['V_down']:.6f} Vup={cst['V_up']:.6f} "
              f"jon={cst['j_on']:.6e} dVc={cst['dV_c']:.6f} regime={cst['regime']} ws60/10={ws:.2f} Ic={sens[f'dvt_scale_{s:g}']['Ic_rD10_lam60_mA']}")
    except Exception as e:
        sens[f"dvt_scale_{s:g}"] = dict(error=repr(e))
        print(f"EXP33SENS s={s:g} ERROR {e!r}")
law._jd, rm.M_of = orig_jd, orig_M
stored = json.loads((ROOT / "results" / "exp32" / "nist_sky130.json").read_text(encoding="utf-8"))["recal"]["corrected"]["measured_transfer"]
base = sens.get("dvt_scale_1", {}).get("constants", {})
out["baseline_vs_stored_Exp32"] = {k: [base.get(k), stored.get(k)] for k in ("V_M", "V_down", "V_up", "j_on", "dV_c")}
out["dvt_sensitivity"] = sens
out["scale_basis"] = "s=0.50: median |dVT_meas|/dVT_BSIM4-clamped at VB=-1.5 over 5 std devices (C33-1); s=0.33: median |gamma| 0.105 / anchor mean slope 0.316"
out["caveat"] = "reverse-bias magnitude used as forward scale proxy; forward body effect at 4 K remains unmeasured (no VB>0 in dataset)"

txt = json.dumps(out, indent=1, ensure_ascii=False, default=float)
dst = ROOT / "results" / "exp33" / "cell33_sensitivity.json"
try:
    dst.write_text(txt, encoding="utf-8")
    print("WROTE", dst)
except Exception as ex:
    alt = Path.home() / "_cell33_out"
    alt.mkdir(parents=True, exist_ok=True)
    (alt / "cell33_sensitivity.json").write_text(txt, encoding="utf-8")
    print(f"WRITE_FAIL {dst} :: {type(ex).__name__}: {ex} -> WROTE {alt / 'cell33_sensitivity.json'}")
print(f"ELAPSED_S={time.time()-t0:.1f}")
