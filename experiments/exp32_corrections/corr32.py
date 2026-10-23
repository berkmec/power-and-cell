"""Exp32: correct the models that failed against measurement and re-test each correction on data NOT used to make it.
 C1 artificial-atom capacitance: weakly gated island -> self-capacitance C = 2 pi eps0 eps_ox d (+ lever arm), tested on
    Shin 2010 SET_A (C_sigma = 0.42 aF, TEM d ~ 2 nm) and SET_B (d ~ 5 nm: one peak at 300 K, clear blockade below 200 K)
 C2 mismatch: Pelgrom A_VT from the measured ~1 mV um / nm rule vs the SKY130 PDK mismatch coefficient (independent)
 C3 second subband vs doping: results of sub32.py (own solver) vs Ray-Martin ~6e12 cm^-2
 C5 sign of dj_high/dV_D in the Exp17 SKY130 model -> registered prediction for how j per cell changes with N
Criteria: results/exp32/predictions_exp32.json (15:26).
usage: python corr32.py -> results/exp32/corrections.json"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp32"
E, EPS0, KB = 1.602176634e-19, 8.8541878128e-12, 1.380649e-23


def c1():
    C = lambda d_nm: 2 * np.pi * EPS0 * 3.9 * d_nm * 1e-9
    seta = dict(d_TEM_nm=2.0, C_pred_aF=C(2.0) * 1e18, C_meas_aF=0.42, dev=C(2.0) / 0.42e-18 - 1)
    seta["passed"] = bool(abs(seta["dev"]) <= 0.15)
    ec = E / C(5.0)                                  # e^2/C in eV (E/C gives volts = eV per e)
    setb = dict(d_nm=5.0, e2_over_C_eV=ec, at300K_kT=ec * E / (KB * 300), at200K_kT=ec * E / (KB * 200))
    setb["passed"] = bool(setb["at300K_kT"] < 7 <= setb["at200K_kT"])
    lever = dict(alpha_G_meas=0.094 / 0.42, alpha_G_exp1=17.549 / 22.989)
    # corrected thresholds (self-capacitance, eps_ox 3.9): e^2/(2C) >= 10 kT
    dmax = lambda T: E**2 / (2 * 10 * KB * T) / (2 * np.pi * EPS0 * 3.9) * 1e9
    thr = {f"{T} K": dmax(T) for T in (300, 77, 4)}
    return dict(SET_A=seta, SET_B=setb, lever_arm=lever, d_max_nm_e2_over_2C_10kT=thr,
                exp1_threshold_300K_nm=0.4, passed=bool(seta["passed"] and setb["passed"]))


def c2():
    area = np.sqrt(0.126 * 0.42)
    rule = np.array([4.15, 5.0])
    pdk = 3.356
    return dict(A_VT_rule_mV_um=rule.tolist(), A_VT_pdk_mV_um=pdk, ratio=(rule / pdk).tolist(), sigma_mV_rule=(rule / area).tolist(),
                sigma_mV_pdk=pdk / area, passed=bool(np.all(rule / pdk <= 1.5)),
                seed_30mV_in_sigma_10um_x_0p15um=[30 / (a / np.sqrt(1.5)) for a in (pdk, 5.0)],
                caveat="PDK vth0_slope biriminin V·µm olduğu varsayıldı (Exp5'te kesinleştirilemedi)")


def c3():
    rows = {}
    for f in sorted(OUT.glob("subband_NA*.json")):
        d = json.loads(f.read_text(encoding="utf-8"))
        rows[f"{d['NA_cm3']:.1e}"] = d["onset_bracket_cm2"]
    return dict(onset_brackets=rows, measured_cm2=6e12, window=[3.75e12, 9.6e12])


def c5():
    sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
    import law
    c = law.constants()
    vd, jh = law.j_high_curve(c=c, n=120)
    slope = np.polyfit(vd, jh, 1)[0]
    return dict(V_down=c["V_down"], V_M=c["V_M"], j_on=c["j_on"], j_high_at_V_down=float(jh[0]), j_high_at_V_M=float(np.interp(c["V_M"], vd, jh)),
                j_high_at_1p2V=float(jh[-1]), dj_dVD_sign=int(np.sign(slope)),
                prediction=("SKY130 modelinde j_high V_D ile artıyor: aynı toplam akımda daha çok (daha dar) hücre, ortalama daha yüksek V_D ve "
                            "daha yüksek j demek → j N ile artmalı. Ölçülen ESD (ısıl, 300 K) aygıtında düşüyor: aygıt sınıfı farkı; SKY130 4 K ölçümü ayırt eder.")
                if slope > 0 else "SKY130 modelinde j_high V_D ile azalıyor: j N ile düşmeli (ESD ölçümüyle aynı yön).")


if __name__ == "__main__":
    res = dict(C1=c1(), C2=c2(), C3=c3(), C5=c5())
    (OUT / "corrections.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float), encoding="utf-8")
    for k, v in res.items():
        print(k, json.dumps(v, ensure_ascii=False, default=float)[:900])
