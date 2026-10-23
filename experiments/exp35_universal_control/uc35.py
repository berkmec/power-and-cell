"""Exp35: universal cell-atom control stack and full-control limits.

Main computation reuses the Exp24 Coulomb-glass bodies and their stored equilibrium-calibrated addresses.
It asks the next question Exp24 did not test: can several qualified atoms be written simultaneously to an
arbitrary bit pattern? Controllers compared under predictions_exp35.json:
  G       best single global body voltage (threshold rule on calibrated addresses)
  L10/20/40  G plus simultaneous target-centred Gaussian local energy offsets (+/-15 meV start)
  CL20    L20 plus up to 5 read-update rounds, then a fresh-RNG verification segment

A second, explicitly independent-bit bootstrap uses stored Exp24 single-atom probabilities to expose register
scaling; its independence assumption is compared with the correlated Coulomb-glass pattern runs above.
No physical lab control is performed. Local offsets model segmented-gate/Stark field modes, not per-atom SKY130 gates.
Output: results/exp35/universal_control35.json (fallback ~/_exp35_out)."""
import itertools
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

_here = Path(__file__).resolve()
ROOT = _here.parents[2] if (_here.parents[2] / "results" / "exp24").exists() else Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp24_atom_program"))
from ap24_atoms import EA0, build, relax  # noqa: E402

KT4 = 8.617333262e-5 * 4.0
WINDOW = (0.805, 1.006)
A0_EV = 0.015
AMP_CLIP_EV = 0.060
GAIN_EV = 0.025
FINAL_SWEEPS = 1200
FEEDBACK_SWEEPS = 600
PROTOCOLS = ("G", "L10", "L20", "L40", "CL20")


def load(rel):
    return json.loads((ROOT / rel).read_text(encoding="utf-8"))


def kernel(r, target_sites, ell):
    if ell is None:
        return np.zeros((len(r), len(target_sites)))
    d2 = ((r[:, None, :] - r[None, target_sites, :]) ** 2).sum(-1)
    return np.exp(-d2 / (2 * ell * ell))


def mc_sample(nocc0, eps, W, targets, sweeps, rng, track_all=False):
    """Exp24-faithful proposal distribution: 1/2 random single flip, 1/2 uniform occupied-to-empty hop.

    Maintained occupied/empty lists make the hop proposal distributionally identical to Exp24's
    rng.choice(np.nonzero(...)) implementation while avoiding repeated full scans."""
    nocc = np.array(nocc0, dtype=int, copy=True)
    Wloc = W
    h = Wloc @ nocc
    N = len(nocc)
    targets = np.array(targets, dtype=int)
    occ = [int(x) for x in np.flatnonzero(nocc == 1)]
    emp = [int(x) for x in np.flatnonzero(nocc == 0)]
    pos_occ = np.full(N, -1, dtype=int)
    pos_emp = np.full(N, -1, dtype=int)
    for q, site in enumerate(occ):
        pos_occ[site] = q
    for q, site in enumerate(emp):
        pos_emp[site] = q

    def remove(lst, pos, site):
        idx = int(pos[site])
        last = lst[-1]
        lst[idx] = last
        pos[last] = idx
        lst.pop()
        pos[site] = -1

    def add(lst, pos, site):
        pos[site] = len(lst)
        lst.append(site)

    def set_site(site, new):
        if nocc[site] == new:
            return
        nocc[site] = new
        if new:
            remove(emp, pos_emp, site)
            add(occ, pos_occ, site)
        else:
            remove(occ, pos_occ, site)
            add(emp, pos_emp, site)

    total = sweeps * N
    kinds = rng.random(total)
    ia = rng.integers(0, N, total)
    hop_u1 = rng.random(total)
    hop_u2 = rng.random(total)
    uacc = rng.random(total)
    burn = sweeps // 3
    target_acc = np.zeros(len(targets), dtype=float)
    all_acc = np.zeros(N, dtype=float) if track_all else None
    snapshots = []
    cnt = 0
    p = 0
    for sweep in range(sweeps):
        for _ in range(N):
            if kinds[p] < 0.5:
                i = int(ia[p])
                dO = (1 - 2 * nocc[i]) * (eps[i] + h[i])
                if dO <= 0 or uacc[p] < np.exp(-dO / KT4):
                    sgn = 1 - 2 * int(nocc[i])
                    new = int(nocc[i]) + sgn
                    h += sgn * Wloc[:, i]
                    set_site(i, new)
            elif occ and emp:
                i = occ[min(int(hop_u1[p] * len(occ)), len(occ) - 1)]
                j = emp[min(int(hop_u2[p] * len(emp)), len(emp) - 1)]
                dO = (eps[j] + h[j] - Wloc[i, j]) - (eps[i] + h[i])
                if dO <= 0 or uacc[p] < np.exp(-dO / KT4):
                    h += Wloc[:, j] - Wloc[:, i]
                    set_site(i, 0)
                    set_site(j, 1)
            p += 1
        if sweep >= burn:
            target_acc += nocc[targets]
            if track_all:
                all_acc += nocc
            cnt += 1
            if (sweep - burn) % 5 == 0:
                snapshots.append(nocc[targets].copy())
    return dict(state=nocc, mean_n=target_acc / cnt, snapshots=np.array(snapshots, dtype=int),
                all_mean=(all_acc / cnt if track_all else None), samples=int(cnt))


def metrics(sample, target_bits):
    target_n = 1 - np.array(target_bits, dtype=int)
    mean_n = sample["mean_n"]
    marginals = np.where(target_n == 1, mean_n, 1.0 - mean_n)
    snaps = sample["snapshots"]
    exact = float(np.mean(np.all(snaps == target_n[None, :], axis=1))) if len(snaps) else 0.0
    modal_n = (mean_n >= 0.5).astype(int)
    return dict(mean_n=[float(x) for x in mean_n], marginals=[float(x) for x in marginals],
                mean_bit_accuracy=float(np.mean(marginals)), marginal_product=float(np.prod(marginals)),
                exact_snapshot_fraction=exact, modal_exact=bool(np.array_equal(modal_n, target_n)),
                modal_bits=[int(1 - x) for x in modal_n])


def best_global(addresses, bits):
    addresses = np.array(addresses, float)
    cand = [WINDOW[0], WINDOW[1]] + [float((a + b) / 2) for a, b in zip(addresses[:-1], addresses[1:])]
    best = None
    for vb in cand:
        vb = float(np.clip(vb, *WINDOW))
        pred = np.where(vb < addresses, 0, 1)
        score = (int(np.sum(pred == np.array(bits))), float(np.min(np.abs(addresses - vb))), -abs(vb - np.mean(addresses)))
        if best is None or score > best[0]:
            best = (score, vb, pred)
    return float(best[1])


def make_patterns(m, island_seed):
    mono = [1] * (m // 2) + [0] * (m - m // 2)  # ascending address crosses B- (0) -> B0 (1) as Vb rises
    alt = [i % 2 for i in range(m)]
    rng = np.random.default_rng(35000 + island_seed)
    rand = alt[:]
    for _ in range(20):
        cand = [int(x) for x in rng.integers(0, 2, m)]
        if 0 < sum(cand) < m and any(cand[i] < cand[i + 1] for i in range(m - 1)):
            rand = cand
            break
    return {"MONO": mono, "ALT": alt, "RAND": rand}


def evaluate(state, eps, W, targets, bits, sweeps, seed, track_all=False):
    sample = mc_sample(state, eps, W, targets, sweeps, np.random.default_rng(seed), track_all=track_all)
    out = metrics(sample, bits)
    out["state"] = sample["state"]
    out["all_mean"] = sample["all_mean"]
    return out


def island_job(run):
    seed = int(run["seed"])
    rows = sorted([x for x in run["rows"] if x["in_cell_window"]], key=lambda x: x["address"])[:4]
    sites = [int(x["site"]) for x in rows]
    addresses = [float(x["address"]) for x in rows]
    r, W, D = build(60.0, False, seed)
    kernels = {ell: kernel(r, sites, ell) for ell in (10.0, 20.0, 40.0)}
    patterns = make_patterns(len(sites), seed)
    pair_d = []
    for i, j in itertools.combinations(sites, 2):
        pair_d.append(float(np.linalg.norm(r[i] - r[j])))
    relaxed_cache = {}
    tasks = []
    for p_idx, (pname, bits) in enumerate(patterns.items()):
        bits = [int(x) for x in bits]
        target_n = 1 - np.array(bits, dtype=int)
        vb = best_global(addresses, bits)
        key = round(vb, 10)
        if key not in relaxed_cache:
            n0, h0 = np.zeros(len(r), int), np.zeros(len(r))
            relax(n0, h0, EA0 + vb + D, W)
            relaxed_cache[key] = n0
        state0 = relaxed_cache[key]
        base_eps = EA0 + vb + D
        seed_base = 450000 + seed * 10000 + p_idx * 1000
        results = {}

        g = evaluate(state0, base_eps, W, sites, bits, FINAL_SWEEPS, seed_base + 900, track_all=True)
        g_modal_all = (g["all_mean"] >= 0.5).astype(int)
        results["G"] = {k: g[k] for k in ("mean_n", "marginals", "mean_bit_accuracy", "marginal_product",
                                          "exact_snapshot_fraction", "modal_exact", "modal_bits")}
        results["G"]["nontarget_modal_flip_fraction_vs_G"] = 0.0

        amps0 = np.where(target_n == 0, A0_EV, -A0_EV)
        for ell in (10.0, 20.0, 40.0):
            name = f"L{int(ell)}"
            ev = evaluate(state0, base_eps + kernels[ell] @ amps0, W, sites, bits,
                          FINAL_SWEEPS, seed_base + 900, track_all=True)
            rec = {k: ev[k] for k in ("mean_n", "marginals", "mean_bit_accuracy", "marginal_product",
                                      "exact_snapshot_fraction", "modal_exact", "modal_bits")}
            modal_all = (ev["all_mean"] >= 0.5).astype(int)
            mask = np.ones(len(r), bool)
            mask[sites] = False
            rec["nontarget_modal_flip_fraction_vs_G"] = float(np.mean(modal_all[mask] != g_modal_all[mask]))
            rec["amplitudes_eV"] = [float(x) for x in amps0]
            results[name] = rec

        # Closed loop at ell=20 nm: read mean occupations, update each target-centred bump, then verify fresh.
        amps = amps0.copy()
        state = np.array(state0, copy=True)
        rounds, updates = [], 0
        for rnd in range(5):
            ev = evaluate(state, base_eps + kernels[20.0] @ amps, W, sites, bits,
                          FEEDBACK_SWEEPS, seed_base + 500 + rnd, track_all=False)
            state = ev["state"]
            rounds.append(dict(round=rnd, mean_n=ev["mean_n"], marginals=ev["marginals"],
                               modal_exact=ev["modal_exact"], amplitudes_eV=[float(x) for x in amps]))
            if ev["modal_exact"] and min(ev["marginals"]) >= 0.90:
                break
            if rnd == 4:
                break
            amps = np.clip(amps + GAIN_EV * np.clip(np.array(ev["mean_n"]) - target_n, -1, 1),
                           -AMP_CLIP_EV, AMP_CLIP_EV)
            updates += 1
        ev = evaluate(state, base_eps + kernels[20.0] @ amps, W, sites, bits,
                      FINAL_SWEEPS, seed_base + 950, track_all=True)
        rec = {k: ev[k] for k in ("mean_n", "marginals", "mean_bit_accuracy", "marginal_product",
                                  "exact_snapshot_fraction", "modal_exact", "modal_bits")}
        modal_all = (ev["all_mean"] >= 0.5).astype(int)
        mask = np.ones(len(r), bool)
        mask[sites] = False
        rec["nontarget_modal_flip_fraction_vs_G"] = float(np.mean(modal_all[mask] != g_modal_all[mask]))
        rec["amplitudes_eV"] = [float(x) for x in amps]
        rec["feedback_updates"] = int(updates)
        rec["feedback_rounds"] = rounds
        results["CL20"] = rec

        tasks.append(dict(island_seed=seed, pattern=pname, bits=bits, target_sites=sites,
                          addresses_V=addresses, global_Vb_V=vb,
                          min_target_distance_nm=min(pair_d) if pair_d else None,
                          protocols=results))
    print(f"island seed {seed}: {len(tasks)} tasks done", flush=True)
    return dict(island_seed=seed, n_atoms=int(len(r)), target_sites=sites, addresses_V=addresses,
                min_target_distance_nm=min(pair_d) if pair_d else None, tasks=tasks)


def agg(tasks, protocol, key, pattern_filter=None):
    vals = [t["protocols"][protocol][key] for t in tasks if pattern_filter is None or t["pattern"] in pattern_filter]
    return float(np.mean(vals)) if vals else None


def modal_frac(tasks, protocol, pattern_filter=None):
    vals = [bool(t["protocols"][protocol]["modal_exact"]) for t in tasks
            if pattern_filter is None or t["pattern"] in pattern_filter]
    return float(np.mean(vals)) if vals else None, int(sum(vals)), int(len(vals))


def bootstrap_scaling():
    small = load("results/exp24/small_bodies.json")
    thermal = load("results/exp24/thermal.json")
    pools = {
        "calibrated_60nm": np.array([row[f"P_target_{b}"] for run in small["runs"] if run["side"] == 60.0
                                      for row in run["rows"] for b in ("0", "1")], dtype=float),
        "uncalibrated_100nm_T4K": np.array([row[f"P_target_{b}"] for run in thermal["runs"] if run["T"] == 4.0
                                             for row in run["rows"] for b in ("0", "1")], dtype=float),
    }
    out = {}
    for name, pool in pools.items():
        rng = np.random.default_rng(46000)
        by_n = {}
        for N in (1, 3, 5, 10, 30, 100, 300):
            p = pool[rng.integers(0, len(pool), size=(20000, N))]
            rec = dict(one_shot_exact_mean=float(np.mean(np.prod(p, axis=1))),
                       modal_exact_fraction=float(np.mean(np.all(p > 0.5, axis=1))))
            for R in (2, 3, 5):
                q = 1 - (1 - p) ** R
                attempts = sum((1 - p) ** j for j in range(R))
                rec[f"retry{R}"] = dict(exact_completion_mean=float(np.mean(np.prod(q, axis=1))),
                                       mean_total_attempts=float(np.mean(np.sum(attempts, axis=1))))
            by_n[str(N)] = rec
        out[name] = dict(pool_size=int(len(pool)), pool_mean=float(np.mean(pool)), pool_min=float(np.min(pool)),
                         by_register_size=by_n)
    return out


def verdicts(tasks, boot):
    mono, nm = ("MONO",), ("ALT", "RAND")
    g_mono = modal_frac(tasks, "G", mono)[0]
    g_alt = modal_frac(tasks, "G", ("ALT",))[0]
    g_rand = modal_frac(tasks, "G", ("RAND",))[0]
    g_nm = modal_frac(tasks, "G", nm)[0]
    uc1 = bool(g_mono >= 0.80 and g_alt <= 0.15 and g_rand <= 0.30 and g_mono - g_nm >= 0.50
               and agg(tasks, "G", "mean_bit_accuracy") >= 0.65)
    l20_nm = modal_frac(tasks, "L20", nm)[0]
    l40_nm = modal_frac(tasks, "L40", nm)[0]
    l10_nm = modal_frac(tasks, "L10", nm)[0]
    uc2 = bool(l20_nm - g_nm >= 0.30 and l20_nm - l40_nm >= 0.10 and ((g_nm <= l10_nm <= l20_nm) or l10_nm >= l40_nm))
    cl_marg = agg(tasks, "CL20", "mean_bit_accuracy", nm)
    cl_exact = agg(tasks, "CL20", "exact_snapshot_fraction", nm)
    cl_modal = modal_frac(tasks, "CL20", nm)[0]
    updates = [t["protocols"]["CL20"]["feedback_updates"] for t in tasks]
    uc3 = bool(cl_marg >= 0.95 and cl_exact >= 0.65 and cl_modal >= 0.80 and float(np.median(updates)) <= 3
               and cl_exact - agg(tasks, "L20", "exact_snapshot_fraction", nm) >= 0.15)
    g_gap = agg(tasks, "G", "marginal_product", nm) - agg(tasks, "G", "exact_snapshot_fraction", nm)
    cl_gap = abs(agg(tasks, "CL20", "exact_snapshot_fraction", nm) - agg(tasks, "CL20", "marginal_product", nm))
    uc4 = bool(g_gap >= 0.10 and cl_gap <= 0.10)
    cal, unc = boot["calibrated_60nm"]["by_register_size"], boot["uncalibrated_100nm_T4K"]["by_register_size"]
    uc5 = bool(0.85 <= cal["10"]["one_shot_exact_mean"] <= 0.95 and 0.30 <= cal["100"]["one_shot_exact_mean"] <= 0.50
               and cal["100"]["retry3"]["exact_completion_mean"] >= 0.90
               and unc["10"]["one_shot_exact_mean"] <= 0.25 and unc["100"]["one_shot_exact_mean"] < 0.001
               and unc["100"]["retry3"]["exact_completion_mean"] < 0.30)
    close = [t for t in tasks if t["min_target_distance_nm"] < 20]
    far = [t for t in tasks if t["min_target_distance_nm"] >= 30]
    if len(close) >= 3 and len(far) >= 3:
        gap_l = agg(far, "L20", "mean_bit_accuracy") - agg(close, "L20", "mean_bit_accuracy")
        gap_cl = abs(agg(far, "CL20", "mean_bit_accuracy") - agg(close, "CL20", "mean_bit_accuracy"))
        uc6 = bool(gap_l >= 0.10 and gap_cl <= 0.05)
        uc6_txt = "pass" if uc6 else "fail"
    else:
        gap_l = gap_cl = None
        uc6_txt = "not testable"
    return dict(UC1="pass" if uc1 else "fail", UC2="pass" if uc2 else "fail", UC3="pass" if uc3 else "fail",
                UC4="pass" if uc4 else "fail", UC5="pass" if uc5 else "fail", UC6=uc6_txt,
                diagnostics=dict(G_modal_MONO=g_mono, G_modal_ALT=g_alt, G_modal_RAND=g_rand, G_modal_nonmono=g_nm,
                                 L10_modal_nonmono=l10_nm, L20_modal_nonmono=l20_nm, L40_modal_nonmono=l40_nm,
                                 CL20_marginal_nonmono=cl_marg, CL20_exact_nonmono=cl_exact, CL20_modal_nonmono=cl_modal,
                                 CL20_median_updates=float(np.median(updates)), G_product_minus_exact=g_gap,
                                 CL20_abs_product_gap=cl_gap, UC6_n_close=len(close), UC6_n_far=len(far),
                                 UC6_L20_gap=gap_l, UC6_CL20_gap=gap_cl))


def aggregate(tasks):
    out = {}
    for protocol in PROTOCOLS:
        mf, k, n = modal_frac(tasks, protocol)
        rec = dict(mean_bit_accuracy=agg(tasks, protocol, "mean_bit_accuracy"),
                   mean_exact_snapshot_fraction=agg(tasks, protocol, "exact_snapshot_fraction"),
                   mean_marginal_product=agg(tasks, protocol, "marginal_product"),
                   modal_exact_fraction=mf, modal_exact_count=f"{k}/{n}",
                   mean_nontarget_modal_flip_fraction=agg(tasks, protocol, "nontarget_modal_flip_fraction_vs_G"),
                   by_pattern={})
        for pat in ("MONO", "ALT", "RAND"):
            mf2, k2, n2 = modal_frac(tasks, protocol, (pat,))
            rec["by_pattern"][pat] = dict(mean_bit_accuracy=agg(tasks, protocol, "mean_bit_accuracy", (pat,)),
                                          mean_exact_snapshot_fraction=agg(tasks, protocol, "exact_snapshot_fraction", (pat,)),
                                          modal_exact_fraction=mf2, modal_exact_count=f"{k2}/{n2}")
        out[protocol] = rec
    return out


def scorecard():
    small = load("results/exp24/small_bodies.json")["summary"]["60 nm"]
    ctrl_times = load("results/exp25/control_all.json")["times"]
    survey = load("results/exp26/vmrev_survey.json")
    isotopes = load("results/exp25/register_isotopes.json")
    return dict(
        stack=["0 census: build threshold/coupling/spectrum map", "1 calibrate in equilibrium at operating point",
               "2 select qualified subset by post-crosstalk margin", "3 compile target into global cut + local offsets/frequency slots",
               "4 write", "5 read and verify", "6 correct local offsets/pulse schedule within bounded retries",
               "7 freeze/hold (remove drive, D0 or deep-margin hold)", "8 couple only after isolation test", "9 refresh/retire defective elements"],
        full_control_gates={
            "address": "post-crosstalk nearest-address margin >= 8 noise units for charge bits (Exp24 corrected criterion) or distinct nonzero quantum couplings with control time << coherence",
            "set": "marginal target probability >=0.95 per element before declaring a register element controlled",
            "hold": "total write/verify latency < retention or coherence budget/10",
            "read": "read error at least 10x below register error budget (not yet measured for charge atoms)",
            "couple": "conditional operation demonstrated without destroying address margins",
            "scale": "whole-register exact success, not mean bit accuracy alone, meets the declared budget"},
        layers=[
            dict(layer="charge atoms in 60 nm cell body (Exp24/35)", stored_single_atom_mean_P=small["mean_P"],
                 qualified_in_window_per_body=small["qualified_in_window"],
                 verdict="best atomic charge layer; multi-bit result is Exp35 Part A; readout is predicted (~7% current step), not simulated as a measured detector"),
            dict(layer="ionic O-vacancy filament (Exp24)", stored="n=1 in 60/60 trials at 0.81-0.90 V and n=2 in 40/40 at 0.95-1.00 V (KMC model)",
                 verdict="atomic-width conductance count is current-programmable in model; atom positions are stochastic, not individually placed"),
            dict(layer="silicon filament cells (Exp18-20)", stored="seed placement +30 mV: 10/10 in seed; step count exact k=2 in 20/20 ring fields at 1.05-1.1 V; winner identity only 65-79% by best predictor",
                 verdict="count and birth site are programmable; arbitrary winner identity needs sequential feedback/repair, not one global step"),
            dict(layer="spin register, purified Si (Exp25/26 + literature)",
                 stored=dict(control_times=ctrl_times, isotopes_at_10kHz={k: v for k, v in isotopes.items() if k.endswith("Dmin=10kHz")}),
                 verdict="only stored route that passes addressability + control-time + storage together at small N (model anchored to measured T2/literature)"),
            dict(layer="spin register, natural Si (Exp25/26)",
                 stored=dict(register_fraction_leak_le_0p1=survey["register_fraction_leak_le_0p1"],
                             dimer_member_fraction=survey["class_fractions"]["member"]),
                 verdict="fails full control by storage leakage from nearest-neighbour 29Si dimers, despite addressability and Lie-algebra controllability in principle"),
        ],
        thermodynamic_margin_law="For independent logistic bits, whole-register success S=prod_i p_i; a common per-bit margin m (logit units) needs approximately m >= ln(N/(1-S)) for large m. Correlated/global drives can fall far below this product (tested in UC4).",
    )


def margin_requirements():
    rows = {}
    for N in (10, 100, 1000):
        rows[str(N)] = {}
        for S in (0.90, 0.99):
            p = S ** (1 / N)
            rows[str(N)][f"S={S:g}"] = dict(per_bit_p_required=float(p), logit_margin=float(np.log(p / (1 - p))))
    return rows


def main():
    t0 = time.time()
    small = load("results/exp24/small_bodies.json")
    eligible = [run for run in small["runs"] if run["side"] == 60.0
                and sum(bool(x["in_cell_window"]) for x in run["rows"]) >= 3]
    print(f"eligible 60 nm islands: {[r['seed'] for r in eligible]}", flush=True)
    with ProcessPoolExecutor(max_workers=6) as ex:
        islands = list(ex.map(island_job, eligible))
    tasks = [t for island in islands for t in island["tasks"]]
    boot = bootstrap_scaling()
    ver = verdicts(tasks, boot)
    out = dict(
        experiment="Exp35 universal cell-atom control stack",
        run_config=dict(FINAL_SWEEPS=FINAL_SWEEPS, FEEDBACK_SWEEPS=FEEDBACK_SWEEPS, A0_eV=A0_EV,
                        GAIN_eV=GAIN_EV, AMP_CLIP_eV=AMP_CLIP_EV, temperature_K=4.0,
                        note="Flip+sparse-hop Metropolis sampler for Exp24 grand-potential model; local offsets are idealised field modes"),
        islands=islands, aggregates=aggregate(tasks), verdicts=ver, bootstrap_scaling=boot,
        margin_requirements=margin_requirements(), universal_control_scorecard=scorecard(),
        consistency=dict(exp24_60nm_mean_P=small["summary"]["60 nm"]["mean_P"],
                         exp24_60nm_n_pairs=small["summary"]["60 nm"]["n_pairs"],
                         exp24_60nm_qualified_in_window=small["summary"]["60 nm"]["qualified_in_window"]),
        wall_s=time.time() - t0,
    )
    print("VERDICTS", {k: v for k, v in ver.items() if k != "diagnostics"}, flush=True)
    print("DIAGNOSTICS", json.dumps(ver["diagnostics"], ensure_ascii=False), flush=True)
    for p in PROTOCOLS:
        a = out["aggregates"][p]
        print(f"AGG {p}: bit={a['mean_bit_accuracy']:.4f} exact_snap={a['mean_exact_snapshot_fraction']:.4f} "
              f"modal={a['modal_exact_fraction']:.3f} ({a['modal_exact_count']}) nontarget_flip={a['mean_nontarget_modal_flip_fraction']:.4f}", flush=True)
    txt = json.dumps(out, indent=1, ensure_ascii=False, default=float)
    dst = ROOT / "results" / "exp35" / "universal_control35.json"
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(txt, encoding="utf-8")
        print("WROTE", dst, flush=True)
    except Exception as ex:
        alt = Path.home() / "_exp35_out"
        alt.mkdir(parents=True, exist_ok=True)
        (alt / "universal_control35.json").write_text(txt, encoding="utf-8")
        print(f"WRITE_FAIL {dst} :: {type(ex).__name__}: {ex} -> WROTE {alt / 'universal_control35.json'}", flush=True)
    print(f"ELAPSED_S={time.time()-t0:.1f}", flush=True)
    return out


if __name__ == "__main__":
    main()

