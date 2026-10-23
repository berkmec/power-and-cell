"""Exp21: can the filament division work at room temperature?
The 4 K model (exp12/17) uses 300 K BSIM4 tables (I_d, I_sub, V_T(V_bs)) plus 4 K knobs: subthreshold slope SS, body-diode
onset V_on and e-folding V_s, impact-ionisation factor F_ii, well leakage g_sub. At 300 K the knobs take room-temperature
values (SS ~ 85 mV/dec, V_s ~ kT/q, F_ii = 1); a floating body needs an isolated p-well (deep n-well, no body tie).
scan    : bistability (three roots of the body equation) and the Maxwell point over 300 K parameter sets, with the high
          state defined by the root structure (not by V_b > 0.5 V as in law.constants)
run     : 1D device runs with a 300 K set that is bistable inside the table range, split test against the law
usage: python rt21.py scan | predict | run [nproc] | analyze   ->  results/exp21/"""
import itertools
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import integrate, optimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
import run as r17  # noqa: E402  (low priority, split_events)
from run import law, sp  # noqa: E402

OUT = ROOT / "results" / "exp21"
VG = np.linspace(0.0, 1.3, 13001)


def roots(vds, p):
    f = law._F(VG, vds, p)
    idx = np.nonzero(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    return [optimize.brentq(law._F, VG[i], VG[i + 1], args=(vds, p)) for i in idx]


def constants_T(p, vd_grid=np.arange(0.40, 2.601, 0.005)):
    """bistable window [V_down, V_up] where F(V_b) has three roots (low, middle, high), Maxwell voltage by equal area"""
    rows = [(v, roots(v, p)) for v in vd_grid]
    bi = [(v, rr) for v, rr in rows if len(rr) >= 3]
    if not bi:
        return dict(bistable=False)
    v_down, v_up = bi[0][0], bi[-1][0]

    def area(v):
        rr = roots(v, p)
        return integrate.quad(lambda x: law._F(x, v, p), rr[0], rr[-1], limit=200)[0] if len(rr) >= 3 else np.nan
    a_lo, a_hi = area(v_down + 1e-4), area(v_up - 1e-4)
    out = dict(bistable=True, V_down=float(v_down), V_up=float(v_up), area_lo=float(a_lo), area_hi=float(a_hi))
    if not (a_lo < 0 < a_hi):
        out.update(regime="no-Maxwell")
        return out
    v_m = optimize.brentq(area, v_down + 1e-4, v_up - 1e-4, xtol=1e-6)
    rr = roots(v_m, p)
    out.update(regime="Maxwell", V_M=float(v_m), dV_c=float(v_m - v_down), j_on=float(law._jd(rr[-1], p)),
               j_low=float(law._jd(rr[0], p)), vb_low=float(rr[0]), vb_mid=float(rr[1]), vb_high=float(rr[-1]))
    return out


# 300 K sets from the scan (Maxwell regime): RT1 = realistic SS 85 mV/dec, RT2 = steeper 75 mV/dec
SETS = dict(RT1=dict(ss=0.085, vs=0.026, fii=1.0, vov0=-0.1, g_sub=1e-11, von=0.70),
            RT2=dict(ss=0.075, vs=0.026, fii=1.0, vov0=-0.1, g_sub=1e-11, von=0.70))
R_D, R_L, W, V_HIGH = 300.0, 160e3, 600.0, 0.22


def jobs():
    out = []
    for lam in (45, 60, 80, 120):
        out.append(dict(name=f"RT1_lam{lam}", set="RT1", lam=float(lam), seed=0, sigma_scale=1.0))
    for s in (1, 2):
        out.append(dict(name=f"RT1_lam60_s{s}", set="RT1", lam=60.0, seed=s, sigma_scale=1.0))
    out.append(dict(name="RT1_lowmm_lam60", set="RT1", lam=60.0, seed=0, sigma_scale=0.1))
    for lam in (60, 80):
        out.append(dict(name=f"RT2_lam{lam}", set="RT2", lam=float(lam), seed=0, sigma_scale=1.0))
    return out


def law_T(c, lam):
    """split law with a conducting low state: the uniform j_low lowers the whole rail and leaves the centre-to-edge
    difference to the contrast dj = j_on - j_low; total width w = (I - j_low W) / dj"""
    dj = c["j_on"] - c["j_low"]
    lc = float(np.sqrt(2 * c["dV_c"] / (R_D * dj)))
    lc0 = float(np.sqrt(2 * c["dV_c"] / (R_D * c["j_on"])))
    ws = lambda l_c: float(2 * lam * np.log(lam / (lam - l_c))) if lam > l_c else float("inf")
    w = ws(lc)
    r_v = R_D * lam**2
    a = w / (2 * lam)
    vdd = c["V_M"] + r_v * c["j_on"] * (1 - np.exp(-2 * a)) / 2 + R_L * (c["j_on"] * w + c["j_low"] * (W - w))
    return dict(dj=dj, lambda_c=lc, w_split=w, lambda_c_jon_only=lc0, w_split_jon_only=ws(lc0), vdd_split=float(vdd))


def predict():
    f = OUT / "predictions.json"
    if f.exists():
        print("predictions already registered (write-once):", f)
        return
    cs = {k: constants_T(p) for k, p in SETS.items()}
    rows = {j["name"]: dict(job=j, **law_T(cs[j["set"]], j["lam"])) for j in jobs()}
    text = ["P1: in every run a filament forms at 300 K and at least one persists over the upper half of the sweep",
            "P2: total filament width = (I - j_low W)/(j_on - j_low) within 15 % (median over steps with filaments)",
            "P3: first split width within 20 % of w_split (dj law) for lambda >= 60 um (7 runs); lambda = 45 um "
            "(1.17 lambda_c, near-critical) within 35 %",
            "P4: centre-to-edge drop at the splits / dV_c (dj in the rail) = 1.0 +- 0.15 in the mean over all splits",
            "P5: the dj law predicts the first split better than the j_on-only law in the mean absolute error"]
    out = dict(registered=time.strftime("%Y-%m-%d %H:%M"), constants=cs, R_D=R_D, R_L=R_L, W=W, v_high=V_HIGH,
               predictions=text, rows=rows)
    f.write_text(json.dumps(out, indent=1), encoding="utf-8")
    for k, v in cs.items():
        print(k, {kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in v.items()})
    for k, v in rows.items():
        print(f"{k:18s} lambda_c {v['lambda_c']:.1f} um  w_split {v['w_split']:.1f} um (j_on-only {v['w_split_jon_only']:.1f})"
              f"  V_DD {v['vdd_split']:.2f} V")


def refine(p):
    """constants_T with V_down and V_up bisected to 1e-5 V (the scan grid is 5 mV)"""
    c = constants_T(p)
    has3 = lambda v: len(roots(v, p)) >= 3
    a, b = c["V_down"] - 0.005, c["V_down"]
    for _ in range(20):
        m = 0.5 * (a + b)
        a, b = (a, m) if has3(m) else (m, b)
    v_down = b
    a, b = c["V_up"], c["V_up"] + 0.005
    for _ in range(20):
        m = 0.5 * (a + b)
        a, b = (m, b) if has3(m) else (a, m)
    v_up = a
    c.update(V_down=v_down, V_up=v_up, dV_c=c["V_M"] - v_down, margin_up=v_up - c["V_M"])
    return c


def jobs2(cs):
    """design-window test: D = r_v (j_on - j_low) / 2 inside (dV_c, V_up - V_M) -> division; above -> nucleation of new
    domains first; below -> no division. r_D = 2 D / (dj lambda^2)."""
    out = []
    for s in ("RT1", "RT2"):
        c = cs[s]
        dj = c["j_on"] - c["j_low"]
        D = dict(inside=0.5 * (c["dV_c"] + c["margin_up"]), above=1.5 * c["margin_up"], below=0.7 * c["dV_c"])
        plan = [("inside", 30, 0.1), ("inside", 50, 0.1), ("above", 30, 0.1)]
        if s == "RT1":
            plan += [("inside", 30, 0.3), ("below", 30, 0.1)]
        for kind, lam, sig in plan:
            r_d = 2 * D[kind] / (dj * lam**2)
            out.append(dict(name=f"W_{s}_{kind}_lam{lam}_mm{sig}", set=s, lam=float(lam), r_D=float(r_d), D=D[kind],
                            kind=kind, seed=0, sigma_scale=sig))
    return out


def predict2():
    f = OUT / "predictions2.json"
    if f.exists():
        print("already registered:", f)
        return
    cs = {k: refine(p) for k, p in SETS.items()}
    rows = {}
    for j in jobs2(cs):
        c = cs[j["set"]]
        lc = j["lam"] * float(np.sqrt(c["dV_c"] / j["D"]))
        w = float(2 * j["lam"] * np.log(j["lam"] / (j["lam"] - lc))) if j["lam"] > lc else float("inf")
        rows[j["name"]] = dict(job=j, lambda_c=lc, w_split=w)
    text = ["Q1 inside (D in the window): domains divide; the first split happens before the 2nd nucleation and at "
            "w_split within 25 %",
            "Q2 above (D > V_up - V_M): new domains nucleate before any split; no split before the 3rd nucleation",
            "Q3 below (D < dV_c): no split over the whole sweep",
            "Q4 w_split(lambda 50)/w_split(lambda 30) inside, mm 0.1: = law ratio within 25 % (both sets)",
            "Q5 inside with mismatch x0.3 (sigma 3.7 mV per 0.5 um): reported, no prediction (tolerance probe)"]
    out = dict(registered=time.strftime("%Y-%m-%d %H:%M"), constants=cs, R_L=R_L, W=W, v_high=V_HIGH, predictions=text,
               rows=rows)
    f.write_text(json.dumps(out, indent=1), encoding="utf-8")
    for k, v in cs.items():
        print(k, {kk: (round(vv, 6) if isinstance(vv, float) else vv) for kk, vv in v.items()})
    for k, v in rows.items():
        print(f"{k:30s} D {v['job']['D']*1e3:5.1f} mV r_D {v['job']['r_D']:7.1f} lambda_c {v['lambda_c']:5.1f} "
              f"w_split {v['w_split']:7.1f}")


def analyze1():
    """round 1: P2 (total high width vs (I - j_low W)/(j_on - j_low), median over steps with filaments)"""
    pr = json.loads((OUT / "predictions.json").read_text(encoding="utf-8"))
    out = {}
    for name, row in pr["rows"].items():
        c = pr["constants"][row["job"]["set"]]
        r = json.loads((OUT / "runs" / f"{name}.json").read_text(encoding="utf-8"))
        err = [(x["frac_high"] * W) / ((x["I_uA"] * 1e-6 - c["j_low"] * W) / (c["j_on"] - c["j_low"])) - 1
               for x in r["hist"] if x["n_fil"] > 0]
        out[name] = dict(P2_median_err=float(np.median(err)), max_n_fil=max(x["n_fil"] for x in r["hist"]),
                         widest=max(max(x["widths_um"] or [0]) for x in r["hist"]), n_splits_counted=len(r["splits"]))
        print(f"{name:18s} P2 median err {out[name]['P2_median_err']:+.3f} | max n {out[name]['max_n_fil']} widest "
              f"{out[name]['widest']:.1f} um | counted splits {out[name]['n_splits_counted']}")
    (OUT / "analysis1.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    fig_, axs = plt.subplots(2, 2, figsize=(13, 9))
    ax = axs[0, 0]
    for p, lab, col in ((law.BASE, "4 K modeli (donmuş kuyu)", "C0"), (SETS["RT1"], "300 K, RT1 (SS 85 mV/dek, g 1e-11)", "C3"),
                        (SETS["RT2"], "300 K, RT2 (SS 75 mV/dek)", "C1")):
        vd = np.arange(0.4, 2.6, 0.002)
        lo, hi = [], []
        for v in vd:
            rr = roots(v, p)
            lo.append(law._jd(rr[0], p) if rr else np.nan)
            hi.append(law._jd(rr[-1], p) if len(rr) >= 3 else np.nan)
        ax.semilogy(vd, np.array(lo) * 1e6, c=col, lw=1.2, label=f"{lab}: alçak dal")
        ax.semilogy(vd, np.array(hi) * 1e6, c=col, lw=2.5, ls="--", label="yüksek dal (bistabil pencere)")
    ax.set_xlabel("V_D (V)")
    ax.set_ylabel("j (µA/µm)")
    ax.set_title("(a) bistabil pencere: 4 K'de ~770 mV, 300 K'de 43–68 mV; akım oranı 10⁴'e karşı 5–12", fontsize=10)
    ax.legend(fontsize=7)
    for ax, name, title in ((axs[0, 1], "RT1_lam60", "(b) bağımsız 0,5 µm hücreler, SKY130 uyumsuzluğu: tuz-biber, bölünme yok"),
                            (axs[1, 0], "W_RT2_inside_lam30_mm0.1", "(c) uyumsuzluk × 0,1, r_v pencere içinde: alanlar çekirdeklenerek çoğalıyor"),
                            (axs[1, 1], "B_RT2_Rsq1e+08_lam50", "(d) dirençli gövde (1e8 Ω/□): tek alan büyüyor, yasa genişliğinde bölünmüyor")):
        r = json.loads((OUT / "runs" / f"{name}.json").read_text(encoding="utf-8"))
        segs = [[(c - w / 2, x["vdd"]), (c + w / 2, x["vdd"])] for x in r["hist"] for c, w in zip(x["centres_um"], x["widths_um"])]
        ax.add_collection(LineCollection(segs, colors="k", linewidths=1.0))
        ax.set_xlim(0, W)
        ax.set_ylim(1.8, 7.5)
        ax.set_xlabel("konum (µm)")
        ax.set_ylabel("V_DD (V)")
        ax.set_title(title, fontsize=9)
    pr2 = json.loads((OUT / "predictions2.json").read_text(encoding="utf-8"))
    for ax, key in ((axs[1, 0], "W_RT2_inside_lam30_mm0.1"), (axs[1, 1], "W_RT2_inside_lam50_mm0.1")):
        ax.text(5, 7.2, f"yasanın bölünme genişliği {pr2['rows'][key]['w_split']:.0f} µm", fontsize=8, color="C3")
    fig_.suptitle("Exp21: standart bölünme oda sıcaklığına inmiyor (SKY130 indirgenmiş model, 300 K)", fontsize=12)
    fig_.tight_layout(rect=(0, 0, 1, 0.96))
    fig_.savefig(OUT / "fig_exp21.png", dpi=95)


def predict3():
    """body coupling: a real 300 K p-well (R_sq ~ 3 kOhm/sq, body coupling length ~2-5 mm) and a hypothetical resistive
    body (1e8 Ohm/sq, ~10-30 um) that averages the mismatch; RT2, D inside the window, full mismatch"""
    f = OUT / "predictions3.json"
    if f.exists():
        print("already registered:", f)
        return
    p2 = json.loads((OUT / "predictions2.json").read_text(encoding="utf-8"))
    rows = {}
    for rsq in (3e3, 1e8):
        for lam in (30, 50):
            base = dict(p2["rows"][f"W_RT2_inside_lam{lam}_mm0.1"])
            j = dict(base["job"], name=f"B_RT2_Rsq{rsq:.0e}_lam{lam}", sigma_scale=1.0, R_sq=rsq)
            rows[j["name"]] = dict(job=j, lambda_c=base["lambda_c"], w_split=base["w_split"])
    text = ["R1 (registered): R_sq = 3e3 (real p-well): the body acts as one node, no spatial domains; the whole device "
            "switches at once (frac_high jumps from <0.05 to >0.9 within one 20 mV step) and n_fil <= 1 throughout",
            "R2 (exploration, no prediction): R_sq = 1e8: do domains form and divide when the body averages the mismatch?"]
    f.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), predictions=text, rows=rows), indent=1),
                 encoding="utf-8")
    print("\n".join(text))


def analyze2(pfile="predictions2.json", ofile="analysis2.json"):
    """per run: order of nucleations (births) and splits, first split width vs the registered law value"""
    pr = json.loads((OUT / pfile).read_text(encoding="utf-8"))
    out = {}
    for name, row in pr["rows"].items():
        r = json.loads((OUT / "runs" / f"{name}.json").read_text(encoding="utf-8"))
        h = r["hist"]
        seq = []
        for x in h:
            seq += ["b"] * x["births"] + ["s"] * x["splits"]
        first = r["splits"][0] if r["splits"] else None
        nb_before = seq.index("s") if "s" in seq else seq.count("b")
        top = max(h, key=lambda x: x["vdd"])
        out[name] = dict(kind=row["job"]["kind"], w_pred=row["w_split"], n_splits=len(r["splits"]),
                         births_before_first_split=nb_before, total_births=seq.count("b"),
                         first_split=dict(vdd=first["vdd"], parent_w=first["parent_w"], kids_w=first["kids_w"],
                                          err=(first["parent_w"] / row["w_split"] - 1) if np.isfinite(row["w_split"]) else None)
                         if first else None,
                         end=dict(vdd=top["vdd"], n_fil=top["n_fil"], widest=max(top["widths_um"] or [0]),
                                  frac_high=top["frac_high"]), sequence="".join(seq)[:60],
                         max_n_fil=max(x["n_fil"] for x in h),
                         max_jump_frac_high=max(abs(b["frac_high"] - a["frac_high"]) for a, b in zip(h, h[1:])))
    (OUT / ofile).write_text(json.dumps(out, indent=1), encoding="utf-8")
    for k, v in out.items():
        fs = v["first_split"]
        fst = (f"first split {fs['vdd']:.2f} V parent {fs['parent_w']:.1f} um -> {[round(w, 1) for w in fs['kids_w']]}"
               + (f" err {fs['err']:+.3f}" if fs["err"] is not None else "")) if fs else "no split"
        print(f"{k:28s} pred {v['w_pred']:6.1f} | births before 1st split {v['births_before_first_split']}, splits "
              f"{v['n_splits']}, births {v['total_births']} | {fst} | end n {v['end']['n_fil']} widest "
              f"{v['end']['widest']:.0f} | max n {v['max_n_fil']} max jump frac_high {v['max_jump_frac_high']:.2f} | "
              f"{v['sequence']}")
    return out


def run(job):
    t0 = time.time()
    p = SETS[job["set"]]
    r_d = job.get("r_D", R_D)
    kw = dict(W=W, dy=0.5, R_sq=job.get("R_sq", 1e12), R_L=R_L, vdd_up=(1.0, 7.5), dvdd=0.02, seed=job["seed"], r_D=r_d,
              contacts="distributed", r_v=r_d * job["lam"] ** 2, down=False, sigma_scale=job["sigma_scale"],
              t_relax=1.0, implicit_load=True, v_high=V_HIGH, **p)
    r = sp.simulate(**kw)
    res = dict(job=job, splits=r17.split_events(r["hist"]), hist=r["hist"], wall=time.time() - t0)
    (OUT / "runs").mkdir(parents=True, exist_ok=True)
    (OUT / "runs" / f"{job['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    return f"{job['name']}: {len(res['splits'])} bölünme, n_max {max(h['n_fil'] for h in r['hist'])} ({res['wall']:.0f} s)"


def scan():
    grid = dict(ss=(0.075, 0.085, 0.095), vs=(0.026, 0.030), von=(0.55, 0.60, 0.65, 0.70),
                g_sub=(1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11), vov0=(-0.1, 0.0), fii=(1.0,))
    keys = list(grid)
    res = []
    t0 = time.time()
    for vals in itertools.product(*[grid[k] for k in keys]):
        p = dict(zip(keys, vals))
        c = constants_T(p)
        res.append(dict(p=p, c=c))
    (OUT / "scan.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    nb = [x for x in res if x["c"].get("bistable")]
    nm = [x for x in nb if x["c"].get("regime") == "Maxwell"]
    print(f"{len(res)} sets, bistable {len(nb)}, Maxwell regime {len(nm)} ({time.time() - t0:.0f} s)")
    for x in sorted(nm, key=lambda x: x["c"]["V_M"])[:12]:
        print(" ", x["p"], {k: round(v, 4) for k, v in x["c"].items() if isinstance(v, float)})
    return res


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[1] == "scan":
        scan()
    elif sys.argv[1] == "predict":
        predict()
    elif sys.argv[1] == "predict2":
        predict2()
    elif sys.argv[1] == "analyze1":
        analyze1()
    elif sys.argv[1] == "fig":
        fig()
    elif sys.argv[1] == "analyze2":
        analyze2()
    elif sys.argv[1] == "predict3":
        predict3()
    elif sys.argv[1] == "analyze3":
        analyze2("predictions3.json", "analysis3.json")
    elif sys.argv[1] in ("run", "run2", "run3"):
        if sys.argv[1] in ("run2", "run3"):
            pr = json.loads((OUT / f"predictions{sys.argv[1][-1]}.json").read_text(encoding="utf-8"))
            js = [v["job"] for v in pr["rows"].values()]
        else:
            js = jobs()
        js = [j for j in js if not (OUT / "runs" / f"{j['name']}.json").exists()]
        with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 3) as pool:
            for line in pool.imap_unordered(run, js):
                print(line, flush=True)
        print("TAMAM", flush=True)
