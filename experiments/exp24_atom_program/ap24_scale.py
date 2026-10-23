"""Exp24 (scale): can a programmable silicon filament cell shrink toward atomic size?
Derivation used here (all inputs from the exp12/17 4 K model and the SKY130 card):
  cell size        w_cell >= 2 lambda_c,  lambda_c = sqrt(2 dV_c / (r_D j_on))         (exp17 law, lambda >> lambda_c)
  programming      a gate seed of amplitude a over the cell region; birth at the unit with the largest overdrive (exp18),
                   units of width dy carry Pelgrom mismatch sigma_u = A_VT / sqrt(dy L)  (A_VT = 3.356 mV um, L = 0.15 um)
                   P_prog = P( max_{m seed units}(a + sigma_u z) > max_{N-m other units}(sigma_u z) )   (extreme values)
  physical floors  independent body unit >= depletion width W_dep; noise: Kramers barrier of the floating body
                   dPhi = (c_b w / q) * integral ln(i_out / i_in) dV  (birth-death quasi-potential), needs dPhi >~ 26 for
                   10 years; atomic limit: N_A a^3 dopants per atomic cell.
Simulation test of P_prog: lambda_c = 36, 3.6, 0.36, 0.036 um (r_D = 10, 1e3, 1e5, 1e7 Ohm/um), dy = lambda_c / 10,
W = 30 lambda_c, lambda = 2 lambda_c, seed a = 30 mV over lambda_c at 0.3 W, 10 mismatch fields each; plus lambda_c = 0.36 um
with a = 150 mV. Only the first birth is used.
usage: python ap24_scale.py derive | predict | run [nproc] | analyze   ->  results/exp24/"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import integrate, optimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp17_silicon_split"))
import run as r17  # noqa: E402  (low priority)
from run import law, sp  # noqa: E402

OUT = ROOT / "results" / "exp24"
Q, KB, EPS0 = 1.602176634e-19, 1.380649e-23, 8.8541878128e-12
A_VT, L_CH, NA = 3.356e-3, 0.15, 1.7e17          # V um, um, cm^-3
CASES = [(36.0, 0.030), (3.6, 0.030), (0.36, 0.030), (0.036, 0.030), (0.36, 0.150)]


def p_prog(a, sigma, m, n, samples=200000, seed=0):
    rng = np.random.default_rng(seed)
    best_seed = a + sigma * rng.standard_normal((samples, m)).max(axis=1)
    rest = sigma * rng.standard_normal((samples, n - m)).max(axis=1) if n - m < 400 else \
        sigma * np.sqrt(2 * np.log(n - m)) * np.ones(samples)
    return float(np.mean(best_seed > rest))


def derive():
    p, c = law.BASE, law.constants()
    VG = np.linspace(0.0, 1.3, 26001)
    vd = c["V_M"]
    f = law._F(VG, vd, p)
    idx = np.nonzero(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    rr = [optimize.brentq(law._F, VG[i], VG[i + 1], args=(vd, p)) for i in idx]
    lo, mid, hi = (rr[0], rr[1], rr[-1]) if len(rr) >= 3 else (0.0, rr[0], rr[-1])
    i_in = lambda V: p["fii"] * law.rm.M_of(vd) * law._jd(V, p)                       # holes generated (A/um)
    i_out = lambda V: 1e-6 * np.exp((V - p["von"]) / p["vs"]) + p["g_sub"] * V          # holes removed (A/um)
    eps = 1e-6
    up = integrate.quad(lambda V: np.log(i_out(V) / i_in(V)), lo + eps, mid, limit=400)[0]     # low -> high barrier
    dn = integrate.quad(lambda V: np.log(i_in(V) / i_out(V)), mid, hi, limit=400)[0]           # high -> low barrier
    c_b = 1e-15                                                                      # F/um (model)
    phi_req = float(np.log(3.15e8 / 1e-3))                                           # 10 years / ~1 ms attempt time
    w_noise = {k: phi_req * Q / (c_b * v) for k, v in (("low_to_high", up), ("high_to_low", dn))}
    eps_si = 11.7 * EPS0
    w_dep = np.sqrt(2 * eps_si * 0.9 / (Q * NA * 1e6)) * 1e6                          # um, body swing 0.9 V
    l_debye = {f"{T} K": np.sqrt(eps_si * KB * T / (Q**2 * NA * 1e6)) * 1e9 for T in (4, 300)}   # nm
    spacing = (NA * 1e6) ** (-1 / 3) * 1e9                                           # nm
    res = dict(
        roots_at_VM=dict(V_M=vd, low=lo, mid=mid, high=hi), barrier_integrals_V=dict(low_to_high=up, high_to_low=dn),
        phi_required_10yr=phi_req, w_noise_um=w_noise, W_dep_um=float(w_dep), debye_length_nm=l_debye,
        acceptor_spacing_nm=float(spacing),
        dopants_per_unit=dict(atomic_cell_0p5nm=NA * (0.5e-7) ** 3, unit_Wdep=NA * (w_dep * 1e-4) * (L_CH * 1e-4) * (w_dep * 1e-4)),
        lambda_c_formula="lambda_c = sqrt(2 dV_c/(r_D j_on))", dV_c=c["dV_c"], j_on=c["j_on"],
        lambda_c_um={f"r_D={r:g}": float(np.sqrt(2 * c["dV_c"] / (r * c["j_on"]))) for r in (10, 1e3, 1e5, 1e7, 1e9)},
        sigma_unit_mV={f"dy={d:g}um": A_VT / np.sqrt(d * L_CH) * 1e3 for d in (3.6, 0.36, 0.08, 0.036, 0.0036, 0.0005)},
        seed_needed_mV_for_P90={}, poly_pitch_um=0.36)
    for d in (3.6, 0.36, 0.08, 0.036, 0.0005):
        s = A_VT / np.sqrt(d * L_CH)
        n, m = 300, 10
        a = optimize.brentq(lambda a_: p_prog(a_, s, m, n, 20000) - 0.9, 0.0, 20 * s)
        res["seed_needed_mV_for_P90"][f"dy={d:g}um"] = a * 1e3
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "derived_scale.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


def jobs():
    out = []
    for lc, a in CASES:
        for f in range(10):
            out.append(dict(name=f"S_lc{lc:g}_a{a * 1e3:g}_f{f}", lc=lc, a=a, field=f))
    return out


def setup(job):
    c = law.constants()
    lc = job["lc"]
    r_D = 2 * c["dV_c"] / (lc**2 * c["j_on"])
    W, dy, lam = 30 * lc, lc / 10, 2 * lc
    R_L = 0.2 / (c["j_on"] * lc)
    return dict(r_D=r_D, W=W, dy=dy, lam=lam, R_L=R_L, seed_c=0.3 * W)


def predict():
    f = OUT / "predictions_scale.json"
    if f.exists():
        print("already registered:", f)
        return
    rows = {}
    for lc, a in CASES:
        s = setup(dict(lc=lc))
        sig = A_VT / np.sqrt(s["dy"] * L_CH)
        P = p_prog(a, sig, 10, 300)
        lo, hi = [int(x) for x in __import__("scipy").stats.binom.interval(0.9, 10, P)]
        rows[f"lc{lc:g}_a{a * 1e3:g}"] = dict(lambda_c=lc, a_mV=a * 1e3, dy_um=s["dy"], sigma_unit_mV=sig * 1e3, P_prog=P,
                                             births_in_seed_90pct_interval=[lo, hi], r_D=s["r_D"])
        print(f"lambda_c {lc:7.3f} um a {a*1e3:4.0f} mV: dy {s['dy']:.4f} um sigma {sig*1e3:6.1f} mV -> P_prog {P:.3f}, "
              f"expected births in seed of 10: {lo}-{hi}")
    text = ["S1: for every (lambda_c, a) case the number of fields (of 10) whose first cell is born on the seed lies in the "
            "registered 90 % binomial interval of P_prog (extreme-value formula with Pelgrom units of width lambda_c/10)",
            "S2: the first cell width scales with lambda_c (first cell / lambda_c the same within a factor 2 for 36, 3.6 "
            "and 0.36 um)",
            "validity: units below the depletion width (~80 nm) and below ~100 stored holes are outside the physics of the "
            "model; lambda_c = 0.036 um is run only as a test of the formula"]
    f.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), rows=rows, predictions=text), indent=1),
                 encoding="utf-8")


def run(job):
    t0 = time.time()
    s = setup(job)
    n = int(round(s["W"] / s["dy"]))
    y = (np.arange(n) + 0.5) * s["dy"]
    prog = np.where(np.abs(y - s["seed_c"]) <= job["lc"] / 2, job["a"], 0.0)
    r = sp.simulate(W=s["W"], dy=s["dy"], R_sq=1e12, R_L=s["R_L"], vdd_up=(1.0, 1.8), dvdd=0.005, seed=100 + job["field"],
                    r_D=s["r_D"], contacts="distributed", r_v=s["r_D"] * s["lam"] ** 2, down=False, prog=prog,
                    t_relax=1.0, implicit_load=True, **law.BASE)
    first = next((h for h in r["hist"] if h["n_fil"] > 0), None)
    res = dict(job=job, setup=s, first=dict(vdd=first["vdd"], centres=first["centres_um"], widths=first["widths_um"])
               if first else None, wall=time.time() - t0)
    if first:
        k = int(np.argmax(first["widths_um"]))
        res["first"]["on_seed"] = bool(abs(first["centres_um"][k] - s["seed_c"]) <= job["lc"])
    (OUT / "scale_runs").mkdir(parents=True, exist_ok=True)
    (OUT / "scale_runs" / f"{job['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    return f"{job['name']}: {res['first']} ({res['wall']:.0f} s)"


def analyze():
    pr = json.loads((OUT / "predictions_scale.json").read_text(encoding="utf-8"))
    out = {}
    for key, row in pr["rows"].items():
        rs = [json.loads(f.read_text(encoding="utf-8")) for f in sorted((OUT / "scale_runs").glob(f"S_{key}_f*.json"))]
        hits = sum(1 for r in rs if r["first"] and r["first"].get("on_seed"))
        w = [max(r["first"]["widths"]) / row["lambda_c"] for r in rs if r["first"]]
        lo, hi = row["births_in_seed_90pct_interval"]
        out[key] = dict(n=len(rs), births_on_seed=hits, predicted_P=row["P_prog"], interval=[lo, hi], S1=lo <= hits <= hi,
                        first_width_over_lambda_c=float(np.median(w)) if w else None, no_birth=sum(1 for r in rs if not r["first"]))
        print(f"{key:14s} on seed {hits}/{len(rs)} (pred P {row['P_prog']:.2f}, interval {lo}-{hi}) S1 "
              f"{'OK' if out[key]['S1'] else 'FAIL'} | first width / lambda_c {out[key]['first_width_over_lambda_c']}")
    (OUT / "analysis_scale.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    cmd = sys.argv[1]
    if cmd == "derive":
        derive()
    elif cmd == "predict":
        predict()
    elif cmd == "analyze":
        analyze()
    elif cmd == "run":
        js = [j for j in jobs() if not (OUT / "scale_runs" / f"{j['name']}.json").exists()]
        with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 3) as pool:
            for line in pool.imap_unordered(run, js):
                print(line, flush=True)
        print("TAMAM", flush=True)
