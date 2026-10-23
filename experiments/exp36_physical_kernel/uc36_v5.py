"""Exp36 v5: surface-segment resolution for background-aware, noisy single-shot atom control.
Registered in results/exp36/predictions_exp36_v5.json after v1-v4 diagnoses and before any v5 outcome.
The v2 controller is unchanged; the physical kernel changes from the archived 4x4 reference to 6x6
(10 nm pitch/sigma) and 8x8 (7.5 nm pitch/sigma) surface-segment Gaussian/depth-decay bases.
ANA uses analog mean occupancy only as an ablation; SS1/SS5 decisions use one noisy terminal snapshot
(p=0.01) or a 3-read majority (p=0.05). Fresh verification metrics never feed the controller.
Exp35 and Exp36 v1-v4 result files are read only for reference and are never overwritten.
Output: results/exp36/universal_control36_v5.json (fallback ~/_exp36_out)."""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent


def _root():
    for p in [SCRIPT_DIR, *SCRIPT_DIR.parents]:
        if (p / "results" / "exp35").exists():
            return p
    return Path(__file__).resolve().parents[2]


ROOT = _root()
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(ROOT / "experiments" / "exp35_universal_control"))
try:
    import uc35 as uc
except ModuleNotFoundError:
    import _uc35_run as uc

ROOT = uc.ROOT
sys.path.insert(0, str(ROOT / "experiments" / "exp24_atom_program"))
from ap24_atoms import EA0, build, relax  # noqa: E402

WINDOW = uc.WINDOW
FINAL_SWEEPS = uc.FINAL_SWEEPS
FEEDBACK_SWEEPS = uc.FEEDBACK_SWEEPS
GUARD_EV = 0.010
STEP_EV = 0.030
D_CLIP_EV = 0.220
SEG_CLIP_V = 0.80
MU_BG = 0.005
CONFIRMATIONS = 5
RIDGE = 1e-3
LZ_NM = 50.0
KERNELS = {
    "SEG6": dict(n_segments_per_side=6, pitch_nm=10.0, sigma_nm=10.0),
    "SEG8": dict(n_segments_per_side=8, pitch_nm=7.5, sigma_nm=7.5),
}
MODE_KINDS = (("ANA", "ANA"), ("SS1", "SS1"), ("SS5", "SS5"), ("SS1U", "SS1-U"))
PROTOCOLS = ("G",) + tuple(f"{kernel}-{suffix}" for kernel in KERNELS for _, suffix in MODE_KINDS)


def physical_basis(r, n_segments_per_side, sigma_nm):
    """Uniform mode + square surface-segment modes at z=0 (site-energy-equivalent volts)."""
    pitch = 60.0 / int(n_segments_per_side)
    centers = (np.arange(int(n_segments_per_side)) + 0.5) * pitch
    cols = [np.ones(len(r))]
    for x0 in centers:
        for y0 in centers:
            cols.append(np.exp(-((r[:, 0] - x0) ** 2 + (r[:, 1] - y0) ** 2) / (2 * float(sigma_nm)**2))
                        * np.exp(-r[:, 2] / LZ_NM))
    return np.column_stack(cols)


def solve_actuator(B, targets, desired, mu_bg, vb_base):
    targets = np.array(targets, dtype=int)
    mask = np.ones(len(B), bool)
    mask[targets] = False
    Bt, Bbg = B[targets], B[mask]
    M = Bt.T @ Bt + float(mu_bg) * (Bbg.T @ Bbg) + RIDGE * np.eye(B.shape[1])
    a = np.linalg.solve(M, Bt.T @ np.array(desired, dtype=float))
    a[0] = float(np.clip(a[0], WINDOW[0] - vb_base, WINDOW[1] - vb_base))
    a[1:] = np.clip(a[1:], -SEG_CLIP_V, SEG_CLIP_V)
    return a


def initial_desired(addresses, bits, vb):
    raw = np.array(addresses, dtype=float) - float(vb)
    return np.where(np.array(bits) == 1, np.maximum(0.0, raw + GUARD_EV), np.minimum(0.0, raw - GUARD_EV))


def sample_metrics(state, vb_base, delta, W, D, targets, bits, sweeps, seed, track_all=False):
    sample = uc.mc_sample(state, EA0 + vb_base + D + delta, W, targets, sweeps,
                          np.random.default_rng(seed), track_all=track_all)
    out = uc.metrics(sample, bits)
    out["state"] = sample["state"]
    out["all_mean"] = sample["all_mean"]
    out["snapshots"] = sample["snapshots"]
    return out


def noisy_read(sample, bits, p_read, n_reads, rng):
    """Estimate bits from terminal snapshot(s); sensor flips occupancy bits independently."""
    snaps = sample["snapshots"]
    chosen = snaps[-min(n_reads, len(snaps)):]
    flips = (rng.random(chosen.shape) < p_read).astype(int)
    observed_n = chosen ^ flips
    est_n = (np.sum(observed_n, axis=0) * 2 > len(chosen)).astype(int)
    est_bits = 1 - est_n
    return [int(x) for x in est_bits], bool(np.array_equal(est_bits, np.array(bits)))


def run_controller(state0, vb_base, B, W, D, targets, bits, addresses, seed_base, kind):
    target_n = 1 - np.array(bits, dtype=int)
    mu = 0.0 if kind == "SS1U" else MU_BG
    desired = initial_desired(addresses, bits, vb_base)
    a = solve_actuator(B, targets, desired, mu, vb_base)
    state = np.array(state0, dtype=int, copy=True)
    rounds, updates, consecutive_exact, stopped = [], 0, 0, False
    max_rounds = 12 if kind == "ANA" else 20
    read_rng = np.random.default_rng(seed_base + 777)
    for rnd in range(max_rounds):
        delta = B @ a
        ev = sample_metrics(state, vb_base, delta, W, D, targets, bits, FEEDBACK_SWEEPS, seed_base + rnd)
        state = ev["state"]
        if kind == "ANA":
            est_bits, observed_exact = ev["modal_bits"], bool(ev["modal_exact"])
        elif kind == "SS5":
            est_bits, observed_exact = noisy_read(ev, bits, 0.05, 3, read_rng)
        else:
            est_bits, observed_exact = noisy_read(ev, bits, 0.01, 1, read_rng)
        rounds.append(dict(round=rnd, estimated_bits=est_bits, observed_exact=observed_exact,
                           mean_n_diagnostic=ev["mean_n"], marginals_diagnostic=ev["marginals"],
                           desired_corrections_eV=[float(x) for x in desired],
                           actuator_V=[float(x) for x in a], effective_global_Vb=float(vb_base + a[0]),
                           target_delta_rms_residual_eV=float(np.sqrt(np.mean((delta[targets] - desired) ** 2))),
                           background_delta_rms_eV=float(np.sqrt(np.mean(delta[np.array([i for i in range(len(B)) if i not in targets])] ** 2)))))
        if kind == "ANA":
            if ev["modal_exact"] and min(ev["marginals"]) >= 0.95:
                stopped = True
                break
            if rnd == max_rounds - 1:
                break
            desired = np.clip(desired + uc.GAIN_EV * np.clip(np.array(ev["mean_n"]) - target_n, -1, 1),
                              -D_CLIP_EV, D_CLIP_EV)
        else:
            consecutive_exact = consecutive_exact + 1 if observed_exact else 0
            if consecutive_exact >= CONFIRMATIONS:
                stopped = True
                break
            if rnd == max_rounds - 1:
                break
            if observed_exact:
                continue  # confirmation round under the same actuator; not a feedback update
            wrong = np.array(est_bits) != np.array(bits)
            step = np.where(np.array(bits) == 1, STEP_EV, -STEP_EV)
            desired = np.clip(desired + wrong * step, -D_CLIP_EV, D_CLIP_EV)
        a = solve_actuator(B, targets, desired, mu, vb_base)
        updates += 1
    return a, desired, state, rounds, updates, stopped


def flip_stats(base_all_mean, ev_all_mean, targets):
    base_modal = (np.array(base_all_mean) >= 0.5).astype(int)
    modal = (np.array(ev_all_mean) >= 0.5).astype(int)
    mask = np.ones(len(base_modal), bool)
    mask[np.array(targets, dtype=int)] = False
    flips = int(np.sum(modal[mask] != base_modal[mask]))
    denom = int(np.sum(mask))
    return dict(off_flips=flips, off_denom=denom, off_flip_frac=float(flips / denom) if denom else 0.0)


def make_record(protocol, island_seed, pat, bits, targets, addresses, depths, base_ev, ev, extra):
    oc = flip_stats(base_ev["all_mean"], ev["all_mean"], targets)
    rec = dict(protocol=protocol, island_seed=int(island_seed), pattern=pat,
               bits=[int(x) for x in bits], targets=[int(x) for x in targets],
               addresses=[round(float(x), 4) for x in addresses],
               target_depths_nm=[round(float(x), 2) for x in depths],
               marginal=float(ev["mean_bit_accuracy"]), min_marginal=float(min(ev["marginals"])),
               marginals=ev["marginals"], exact_snap=float(ev["exact_snapshot_fraction"]),
               modal_exact=bool(ev["modal_exact"]), modal_bits=ev["modal_bits"], mean_n=ev["mean_n"],
               off_flips=oc["off_flips"], off_denom=oc["off_denom"], off_flip_frac=oc["off_flip_frac"],
               n_atoms=int(len(base_ev["all_mean"])))
    rec.update(extra)
    return rec


def island_job(run):
    seed = int(run["seed"])
    rows = sorted([x for x in run["rows"] if x["in_cell_window"]], key=lambda x: x["address"])[:4]
    targets = [int(x["site"]) for x in rows]
    addresses = [float(x["address"]) for x in rows]
    r, W, D = build(60.0, False, seed)
    bases = {name: physical_basis(r, cfg["n_segments_per_side"], cfg["sigma_nm"])
             for name, cfg in KERNELS.items()}
    depths = r[targets, 2]
    patterns = uc.make_patterns(len(targets), seed)
    relaxed_cache = {}

    def relaxed(vb):
        key = round(float(vb), 10)
        if key not in relaxed_cache:
            n0, h0 = np.zeros(len(r), int), np.zeros(len(r))
            relax(n0, h0, EA0 + float(vb) + D, W)
            relaxed_cache[key] = n0
        return relaxed_cache[key]

    records, task_vb = [], {}
    for pi, (pat, bits0) in enumerate(patterns.items()):
        bits = [int(x) for x in bits0]
        vb_g = float(uc.best_global(addresses, bits))
        task_vb[pat] = vb_g
        state_g = relaxed(vb_g)
        seed_pattern = 360_000_000 + seed * 1_000_000 + pi * 100_000
        g_ev = uc.evaluate(state_g, EA0 + vb_g + D, W, targets, bits, FINAL_SWEEPS,
                           seed_pattern + 900, track_all=True)
        records.append(make_record("G", seed, pat, bits, targets, addresses, depths, g_ev, g_ev,
                                   dict(kernel=None, vb=vb_g, effective_global_Vb=vb_g, feedback_updates=0,
                                        controller_stopped=None, false_stop=None, rounds=[], actuator_V=None)))
        for ki, (kernel, cfg) in enumerate(KERNELS.items()):
            B = bases[kernel]
            for mi, (kind, suffix) in enumerate(MODE_KINDS):
                protocol = f"{kernel}-{suffix}"
                seed_base = seed_pattern + ki * 10_000 + mi * 1_000
                a, desired, state, rounds, updates, stopped = run_controller(
                    state_g, vb_g, B, W, D, targets, bits, addresses, seed_base, kind)
                ev = sample_metrics(state, vb_g, B @ a, W, D, targets, bits, FINAL_SWEEPS,
                                    seed_base + 950, track_all=True)
                extra = dict(kernel=kernel, segment_pitch_nm=float(cfg["pitch_nm"]),
                             segment_sigma_nm=float(cfg["sigma_nm"]), vb=vb_g,
                             effective_global_Vb=float(vb_g + a[0]), feedback_updates=int(updates),
                             controller_stopped=bool(stopped), false_stop=bool(stopped and not ev["modal_exact"]),
                             rounds=rounds, actuator_V=[float(x) for x in a],
                             desired_corrections_eV=[float(x) for x in desired],
                             achieved_target_delta_eV=[float(x) for x in (B @ a)[targets]],
                             background_penalty_mu=(0.0 if kind == "SS1U" else MU_BG))
                records.append(make_record(protocol, seed, pat, bits, targets, addresses, depths, g_ev, ev, extra))
    island = dict(island_seed=seed, side_nm=60.0, n_atoms=int(len(r)), targets=[int(x) for x in targets],
                  target_xy_nm=[[float(x), float(y)] for x, y in r[targets, :2]],
                  addresses=[float(x) for x in addresses], target_depths_nm=[float(x) for x in depths],
                  task_global_Vb=task_vb, vb_med=float(np.median(addresses)),
                  kernels={name: dict(**cfg, basis_columns=int(bases[name].shape[1]),
                                      basis_condition_number=float(np.linalg.cond(bases[name])),
                                      target_segment_max_response=[float(x) for x in np.max(np.abs(bases[name][targets, 1:]), axis=1)])
                           for name, cfg in KERNELS.items()})
    print(f"Exp36 island seed {seed}: {len(records)} protocol records done", flush=True)
    return island, records


def aggregate_local(records, protocol):
    def one(sel):
        if not sel:
            return dict(n_tasks=0, marginal=None, exact_snap=None, modal=None, modal_count="0/0", off_flip_frac=None)
        modal_k = int(sum(bool(x["modal_exact"]) for x in sel))
        return dict(n_tasks=int(len(sel)), marginal=float(np.mean([x["marginal"] for x in sel])),
                    exact_snap=float(np.mean([x["exact_snap"] for x in sel])),
                    modal=float(modal_k / len(sel)), modal_count=f"{modal_k}/{len(sel)}",
                    off_flip_frac=float(np.mean([x["off_flip_frac"] for x in sel])))
    sel = [x for x in records if x["protocol"] == protocol]
    return dict(all=one(sel), nonmono=one([x for x in sel if x["pattern"] != "MONO"]),
                mono=one([x for x in sel if x["pattern"] == "MONO"]))


def extra_stats(records, protocol):
    sel = [x for x in records if x["protocol"] == protocol]
    non = [x for x in sel if x["pattern"] != "MONO"]
    return dict(median_updates_nonmono=float(np.median([x["feedback_updates"] for x in non])) if non else None,
                false_stops=int(sum(bool(x.get("false_stop")) for x in sel)),
                stopped=int(sum(bool(x.get("controller_stopped")) for x in sel)))


def depth_stats(records, ana_protocol="SEG6-ANA", ss1_protocol="SEG6-SS1"):
    ana = {(x["island_seed"], x["pattern"]): x for x in records if x["protocol"] == ana_protocol}
    ss1 = {(x["island_seed"], x["pattern"]): x for x in records if x["protocol"] == ss1_protocol}
    rows = []
    for key, rec in ana.items():
        first = rec["rounds"][0]["marginals_diagnostic"] if rec["rounds"] else []
        for j, z in enumerate(rec["target_depths_nm"]):
            rows.append(dict(island_seed=rec["island_seed"], pattern=rec["pattern"], target_index=j, depth_nm=float(z),
                             first_round_marginal=float(first[j]) if j < len(first) else None,
                             final_ss1_marginal=float(ss1[key]["marginals"][j])))
    out = dict(rows=rows)
    for name, pred in (("shallow_lt75", lambda z: z < 75), ("deep_ge75", lambda z: z >= 75)):
        sel = [x for x in rows if pred(x["depth_nm"])]
        out[name] = dict(n_target_task_bits=len(sel),
                         first_round_mean=float(np.mean([x["first_round_marginal"] for x in sel])) if sel else None,
                         final_ss1_mean=float(np.mean([x["final_ss1_marginal"] for x in sel])) if sel else None)
    out["testable"] = out["shallow_lt75"]["n_target_task_bits"] >= 8 and out["deep_ge75"]["n_target_task_bits"] >= 8
    if out["testable"]:
        out["first_round_gap_shallow_minus_deep"] = out["shallow_lt75"]["first_round_mean"] - out["deep_ge75"]["first_round_mean"]
        out["final_ss1_gap_shallow_minus_deep"] = out["shallow_lt75"]["final_ss1_mean"] - out["deep_ge75"]["final_ss1_mean"]
    return out


def main():
    t0 = time.time()
    islands, records = [], []
    small = uc.load("results/exp24/small_bodies.json")
    eligible = [run for run in small["runs"] if run["side"] == 60.0
                and sum(bool(x["in_cell_window"]) for x in run["rows"]) >= 3]
    print(f"eligible 60 nm islands: {[r['seed'] for r in eligible]}", flush=True)
    with ProcessPoolExecutor(max_workers=7) as ex:
        for island, recs in ex.map(island_job, eligible):
            islands.append(island)
            records.extend(recs)
    aggs = {p: {**aggregate_local(records, p), **extra_stats(records, p)} for p in PROTOCOLS}
    ref_path = Path(ROOT) / "results" / "exp35" / "universal_control35_v3.json"
    cl20m_bg = None
    if ref_path.exists():
        cl20m_bg = json.loads(ref_path.read_text(encoding="utf-8"))["aggregates"]["CL20M"]["mean_nontarget_modal_flip_fraction"]
    ds = depth_stats(records)
    def budget_for(protocol):
        segs = np.array([c for x in records if x["protocol"] == protocol for c in x["actuator_V"][1:]])
        return dict(n_segment_coefficients=int(len(segs)),
                    frac_within_0p50V=float(np.mean(np.abs(segs) <= 0.50)),
                    frac_at_clip=float(np.mean(np.abs(segs) >= SEG_CLIP_V - 1e-12)),
                    max_abs_segment_V=float(np.max(np.abs(segs))),
                    all_effective_global_in_window=bool(all(WINDOW[0] <= x["effective_global_Vb"] <= WINDOW[1]
                                                           for x in records if x["protocol"] == protocol)))
    actuator_budget = {kernel: budget_for(f"{kernel}-SS1") for kernel in KERNELS}
    def passes(p, marg, exact, modal):
        a = aggs[p]["nonmono"]
        return bool(a["marginal"] >= marg and a["exact_snap"] >= exact and a["modal"] >= modal
                    and aggs[p]["all"]["off_flip_frac"] <= 0.05)
    full_by_kernel = {kernel: bool(passes(f"{kernel}-SS1", 0.95, 0.65, 0.80)
                                   and passes(f"{kernel}-SS5", 0.95, 0.65, 0.80))
                      for kernel in KERNELS}
    coarsest_full = next((kernel for kernel in KERNELS if full_by_kernel[kernel]), None)
    modal_gap6 = 100 * (aggs["SEG6-SS1"]["nonmono"]["modal"] - aggs["SEG6-SS5"]["nonmono"]["modal"])
    seg8_resolution_ok = bool(full_by_kernel["SEG8"]
                              and aggs["SEG8-SS1"]["nonmono"]["modal"] >= aggs["SEG6-SS1"]["nonmono"]["modal"]
                              and aggs["SEG8-SS1"]["all"]["off_flip_frac"] <= aggs["SEG6-SS1"]["all"]["off_flip_frac"])
    actuator_ok = bool(all(b["frac_within_0p50V"] >= .90 and b["frac_at_clip"] < .10
                           and b["all_effective_global_in_window"] for b in actuator_budget.values()))
    verdicts = {
        "RY1": "PASS" if passes("SEG6-ANA", .95, .65, .80) and (cl20m_bg is None or aggs["SEG6-ANA"]["all"]["off_flip_frac"] <= cl20m_bg - .015) else "FAIL",
        "RY2": "PASS" if passes("SEG6-SS1", .95, .65, .80) and aggs["SEG6-SS1"]["median_updates_nonmono"] <= 10 and aggs["SEG6-SS1"]["false_stops"] <= 1 else "FAIL",
        "RY3": "PASS" if passes("SEG6-SS5", .90, .55, .70) and modal_gap6 <= 15 else "FAIL",
        "RY4": "PASS" if seg8_resolution_ok else "FAIL",
        "RY5": ("PASS" if ds["first_round_gap_shallow_minus_deep"] >= .10 and ds["final_ss1_gap_shallow_minus_deep"] <= .05 else "FAIL") if ds["testable"] else "NOT TESTABLE",
        "RY6": "PASS" if actuator_ok else "FAIL",
        "full_control_by_kernel": full_by_kernel,
        "coarsest_full_control_kernel": coarsest_full,
        "full_control_under_noisy_readout": bool(coarsest_full),
    }
    out = dict(exp="36 v5", date="2026-10-17", task="Surface-segment resolution for background-aware, noisy single-shot universal atom control",
               registered="predictions_exp36_v5.json (after v1-v4 diagnoses, before v5 outcomes)",
               config=dict(window=WINDOW, final_sweeps=FINAL_SWEEPS, feedback_sweeps=FEEDBACK_SWEEPS, guard_eV=GUARD_EV,
                           step_eV=STEP_EV, desired_clip_eV=D_CLIP_EV, segment_clip_V=SEG_CLIP_V, mu_bg=MU_BG,
                           confirmation_rounds=CONFIRMATIONS, ridge=RIDGE, kernels=KERNELS,
                           depth_decay_nm=LZ_NM, protocols=list(PROTOCOLS)),
               sampler_gate_from_exp35="Corrected Exp24-faithful sampler passed the Exp35 Exp24-anchor gate; same module reused unchanged.",
               exp35_reference=dict(CL20M_background_overall=cl20m_bg), verdicts=verdicts, aggregates=aggs,
               depth_stats=ds, actuator_budget=actuator_budget, islands=islands, records=records,
               runtime_s=round(time.time() - t0, 1))
    outdir = Path(ROOT) / "results" / "exp36"
    try:
        outdir.mkdir(parents=True, exist_ok=True)
    except PermissionError:
        outdir = Path.home() / "_exp36_out"
        outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "universal_control36_v5.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("WROTE", outdir / "universal_control36_v5.json", "runtime", out["runtime_s"], "verdicts", verdicts, flush=True)
    for p in PROTOCOLS:
        a, aa = aggs[p]["nonmono"], aggs[p]["all"]
        print(p, "nonmono", a, "bg_all", aa["off_flip_frac"], "extra", {k: aggs[p][k] for k in ("median_updates_nonmono", "false_stops", "stopped")}, flush=True)
    print("depth_counts shallow/deep", ds["shallow_lt75"]["n_target_task_bits"], ds["deep_ge75"]["n_target_task_bits"], flush=True)


if __name__ == "__main__":
    main()

