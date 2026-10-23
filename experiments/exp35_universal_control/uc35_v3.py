"""Exp35 v3 follow-up: global-operating-point co-optimisation and anneal-assisted closed-loop control.
Registered in results/exp35/predictions_exp35_v3.json after the v2 diagnosis and before any v3 outcome.
Reuses the corrected Exp24-faithful sampler/metrics from the v2 module (uc35.py in the project copy,
_uc35_run.py in the session copy). Protocols: G baseline, L20M (median Vb + address-set local bumps),
CL20M (local + global-trim feedback), CL20A (CL20M + zero-T settle before each measurement/verification).
No physical lab control is performed; local bumps remain idealised segmented-gate/Stark field modes.
Output: results/exp35/universal_control35_v3.json (fallback ~/_exp35_out)."""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
try:
    import uc35 as uc
except ModuleNotFoundError:
    import _uc35_run as uc

ROOT = uc.ROOT
sys.path.insert(0, str(ROOT / "experiments" / "exp24_atom_program"))
from ap24_atoms import EA0, build, relax  # noqa: E402

GUARD_EV = 0.005
VB_GAIN_V = 0.015
PROTOCOLS3 = ("G", "L20M", "CL20M", "CL20A")
WINDOW = uc.WINDOW
FINAL_SWEEPS = uc.FINAL_SWEEPS
FEEDBACK_SWEEPS = uc.FEEDBACK_SWEEPS


def settle(state, eps, W):
    s = np.array(state, dtype=int, copy=True)
    h = W @ s
    relax(s, h, eps, W)
    return s


def eval3(state, vb, delta, W, D, targets, bits, sweeps, seed, track_all=False, do_settle=False):
    eps = EA0 + vb + D + delta
    st = settle(state, eps, W) if do_settle else state
    sample = uc.mc_sample(st, eps, W, targets, sweeps, np.random.default_rng(seed), track_all=track_all)
    out = uc.metrics(sample, bits)
    out["state"] = sample["state"]
    out["all_mean"] = sample["all_mean"]
    out["Vb_V"] = float(vb)
    return out


def rec_of(ev, g_modal_all, target_mask):
    rec = {k: ev[k] for k in ("mean_n", "marginals", "mean_bit_accuracy", "marginal_product",
                              "exact_snapshot_fraction", "modal_exact", "modal_bits")}
    if ev.get("all_mean") is not None and g_modal_all is not None:
        modal_all = (ev["all_mean"] >= 0.5).astype(int)
        rec["nontarget_modal_flip_fraction_vs_G"] = float(np.mean(modal_all[target_mask] != g_modal_all[target_mask]))
    rec["Vb_V"] = float(ev["Vb_V"])
    return rec


def initial_amps(addresses, bits, vb):
    base = np.array(addresses, dtype=float) - float(vb)
    signs = np.where(np.array(bits) == 1, GUARD_EV, -GUARD_EV)
    return np.clip(base + signs, -uc.AMP_CLIP_EV, uc.AMP_CLIP_EV)


def closed_loop(state0, vb0, amps0, K, W, D, targets, bits, seed_base, anneal):
    target_n = 1 - np.array(bits, dtype=int)
    vb = float(vb0)
    amps = np.array(amps0, dtype=float, copy=True)
    state = np.array(state0, dtype=int, copy=True)
    rounds, updates = [], 0
    for rnd in range(5):
        ev = eval3(state, vb, K @ amps, W, D, targets, bits, FEEDBACK_SWEEPS,
                   seed_base + rnd, track_all=False, do_settle=anneal)
        state = ev["state"]
        rounds.append(dict(round=rnd, Vb_V=vb, mean_n=ev["mean_n"], marginals=ev["marginals"],
                           modal_exact=ev["modal_exact"], amplitudes_eV=[float(x) for x in amps]))
        if ev["modal_exact"] and min(ev["marginals"]) >= 0.90:
            break
        if rnd == 4:
            break
        err = np.array(ev["mean_n"]) - target_n
        amps = np.clip(amps + uc.GAIN_EV * np.clip(err, -1, 1), -uc.AMP_CLIP_EV, uc.AMP_CLIP_EV)
        vb = float(np.clip(vb + VB_GAIN_V * float(np.mean(err)), *WINDOW))
        updates += 1
    return vb, amps, state, rounds, updates


def island_job(run):
    seed = int(run["seed"])
    rows = sorted([x for x in run["rows"] if x["in_cell_window"]], key=lambda x: x["address"])[:4]
    sites = [int(x["site"]) for x in rows]
    addresses = [float(x["address"]) for x in rows]
    r, W, D = build(60.0, False, seed)
    K20 = uc.kernel(r, sites, 20.0)
    patterns = uc.make_patterns(len(sites), seed)
    relaxed_cache = {}

    def relaxed(vb):
        key = round(float(vb), 10)
        if key not in relaxed_cache:
            n0, h0 = np.zeros(len(r), int), np.zeros(len(r))
            relax(n0, h0, EA0 + float(vb) + D, W)
            relaxed_cache[key] = n0
        return relaxed_cache[key]

    mask = np.ones(len(r), bool)
    mask[sites] = False
    tasks = []
    for p_idx, (pname, bits0) in enumerate(patterns.items()):
        bits = [int(x) for x in bits0]
        seed_base = 450000 + seed * 10000 + p_idx * 1000
        vb_g = uc.best_global(addresses, bits)
        g = uc.evaluate(relaxed(vb_g), EA0 + vb_g + D, W, sites, bits,
                        FINAL_SWEEPS, seed_base + 900, track_all=True)
        g_modal_all = (g["all_mean"] >= 0.5).astype(int)
        results = {"G": rec_of(dict(g, Vb_V=vb_g), g_modal_all, mask)}
        results["G"]["nontarget_modal_flip_fraction_vs_G"] = 0.0

        vb_med = float(np.median(addresses))
        state_med = relaxed(vb_med)
        amps0 = initial_amps(addresses, bits, vb_med)
        ev = eval3(state_med, vb_med, K20 @ amps0, W, D, sites, bits,
                   FINAL_SWEEPS, seed_base + 910, track_all=True, do_settle=False)
        rec = rec_of(ev, g_modal_all, mask)
        rec["amplitudes_eV"] = [float(x) for x in amps0]
        results["L20M"] = rec

        for name, anneal, fb_off, ver_off in (("CL20M", False, 600, 960), ("CL20A", True, 700, 970)):
            vb_f, amps_f, state_f, rounds, updates = closed_loop(
                state_med, vb_med, amps0, K20, W, D, sites, bits, seed_base + fb_off, anneal)
            ev = eval3(state_f, vb_f, K20 @ amps_f, W, D, sites, bits,
                       FINAL_SWEEPS, seed_base + ver_off, track_all=True, do_settle=anneal)
            rec = rec_of(ev, g_modal_all, mask)
            rec["amplitudes_eV"] = [float(x) for x in amps_f]
            rec["feedback_updates"] = int(updates)
            rec["feedback_rounds"] = rounds
            results[name] = rec

        tasks.append(dict(island_seed=seed, pattern=pname, bits=bits, target_sites=sites,
                          addresses_V=addresses, global_Vb_V=float(vb_g), median_Vb_V=vb_med,
                          min_target_distance_nm=float(min(np.linalg.norm(r[i] - r[j])
                                                            for a, i in enumerate(sites) for j in sites[a + 1:])),
                          protocols=results))
    print(f"v3 island seed {seed}: {len(tasks)} tasks done", flush=True)
    return dict(island_seed=seed, tasks=tasks)


def mean_key(tasks, protocol, key, patterns=None):
    vals = [t["protocols"][protocol][key] for t in tasks if patterns is None or t["pattern"] in patterns]
    return float(np.mean(vals)) if vals else None


def modal_count(tasks, protocol, patterns=None):
    vals = [bool(t["protocols"][protocol]["modal_exact"]) for t in tasks
            if patterns is None or t["pattern"] in patterns]
    return int(sum(vals)), int(len(vals)), float(np.mean(vals)) if vals else None


def aggregate3(tasks):
    out = {}
    for protocol in PROTOCOLS3:
        k, n, mf = modal_count(tasks, protocol)
        rec = dict(mean_bit_accuracy=mean_key(tasks, protocol, "mean_bit_accuracy"),
                   mean_exact_snapshot_fraction=mean_key(tasks, protocol, "exact_snapshot_fraction"),
                   mean_marginal_product=mean_key(tasks, protocol, "marginal_product"),
                   modal_exact_fraction=mf, modal_exact_count=f"{k}/{n}",
                   mean_nontarget_modal_flip_fraction=mean_key(tasks, protocol, "nontarget_modal_flip_fraction_vs_G"),
                   by_pattern={})
        for pat in ("MONO", "ALT", "RAND"):
            k2, n2, mf2 = modal_count(tasks, protocol, (pat,))
            rec["by_pattern"][pat] = dict(mean_bit_accuracy=mean_key(tasks, protocol, "mean_bit_accuracy", (pat,)),
                                          mean_exact_snapshot_fraction=mean_key(tasks, protocol, "exact_snapshot_fraction", (pat,)),
                                          modal_exact_fraction=mf2, modal_exact_count=f"{k2}/{n2}")
        out[protocol] = rec
    return out


def verdicts3(tasks):
    nm = ("ALT", "RAND")
    l_k, _, l_mf = modal_count(tasks, "L20M", nm)
    m_k, _, m_mf = modal_count(tasks, "CL20M", nm)
    a_k, _, a_mf = modal_count(tasks, "CL20A", nm)
    alt_k, _, _ = modal_count(tasks, "CL20A", ("ALT",))
    upd_m = float(np.median([t["protocols"]["CL20M"]["feedback_updates"] for t in tasks if t["pattern"] in nm]))
    amps_a = np.array([a for t in tasks for a in t["protocols"]["CL20A"]["amplitudes_eV"]])
    vbs_a = [t["protocols"]["CL20A"]["Vb_V"] for t in tasks]
    diag = dict(
        L20M_nonmono_modal_count=f"{l_k}/14", L20M_nonmono_bit=mean_key(tasks, "L20M", "mean_bit_accuracy", nm),
        CL20M_nonmono_marginal=mean_key(tasks, "CL20M", "mean_bit_accuracy", nm),
        CL20M_nonmono_exact=mean_key(tasks, "CL20M", "exact_snapshot_fraction", nm),
        CL20M_nonmono_modal_count=f"{m_k}/14", CL20M_median_updates_nonmono=upd_m,
        CL20A_nonmono_marginal=mean_key(tasks, "CL20A", "mean_bit_accuracy", nm),
        CL20A_nonmono_exact=mean_key(tasks, "CL20A", "exact_snapshot_fraction", nm),
        CL20A_nonmono_modal_count=f"{a_k}/14", CL20A_ALT_modal_count=f"{alt_k}/7",
        CL20A_minus_CL20M_modal=float(a_mf - m_mf),
        CL20A_nontarget_flip=mean_key(tasks, "CL20A", "nontarget_modal_flip_fraction_vs_G"),
        CL20M_nontarget_flip=mean_key(tasks, "CL20M", "nontarget_modal_flip_fraction_vs_G"),
        CL20A_amp_within_45mV_fraction=float(np.mean(np.abs(amps_a) <= 0.045)),
        CL20A_Vb_range=[float(min(vbs_a)), float(max(vbs_a))],
        G_consistency_note="G recomputed with the corrected sampler and same final seed offset as v2; compare with v2 aggregates in the summary")
    return dict(
        UC3_1="pass" if (l_k >= 8 and diag["L20M_nonmono_bit"] >= 0.85) else "fail",
        UC3_2="pass" if (diag["CL20M_nonmono_marginal"] >= 0.90 and diag["CL20M_nonmono_exact"] >= 0.55
                         and m_mf >= 0.65 and upd_m <= 3) else "fail",
        UC3_3="pass" if (diag["CL20A_nonmono_marginal"] >= 0.95 and diag["CL20A_nonmono_exact"] >= 0.65
                         and a_mf >= 0.80 and a_mf - m_mf >= 0.10) else "fail",
        UC3_4="pass" if alt_k >= 5 else "fail",
        UC3_5="pass" if (diag["CL20A_nontarget_flip"] <= 0.05
                         and diag["CL20A_nontarget_flip"] - diag["CL20M_nontarget_flip"] <= 0.02) else "fail",
        UC3_6="pass" if (diag["CL20A_amp_within_45mV_fraction"] >= 0.80
                         and WINDOW[0] <= diag["CL20A_Vb_range"][0] and diag["CL20A_Vb_range"][1] <= WINDOW[1]) else "fail",
        diagnostics=diag)


def main():
    t0 = time.time()
    small = uc.load("results/exp24/small_bodies.json")
    eligible = [run for run in small["runs"] if run["side"] == 60.0
                and sum(bool(x["in_cell_window"]) for x in run["rows"]) >= 3]
    with ProcessPoolExecutor(max_workers=6) as ex:
        islands = list(ex.map(island_job, eligible))
    tasks = [t for island in islands for t in island["tasks"]]
    v2 = uc.load("results/exp35/universal_control35.json")
    ver = verdicts3(tasks)
    out = dict(experiment="Exp35 v3: global-trim + anneal-assisted closed-loop control",
               registered="predictions_exp35_v3.json (after v2 diagnosis, before v3 outcomes)",
               run_config=dict(GUARD_eV=GUARD_EV, VB_GAIN_V=VB_GAIN_V, FINAL_SWEEPS=FINAL_SWEEPS,
                               FEEDBACK_SWEEPS=FEEDBACK_SWEEPS, ell_nm=20.0,
                               sampler="Exp24-faithful flip/hop sampler from v2 module"),
               islands=islands, aggregates=aggregate3(tasks), verdicts=ver,
               v2_reference=dict(verdicts=v2["verdicts"], aggregates=v2["aggregates"]),
               wall_s=time.time() - t0)
    print("V3 VERDICTS", {k: v for k, v in ver.items() if k != "diagnostics"}, flush=True)
    print("V3 DIAGNOSTICS", json.dumps(ver["diagnostics"], ensure_ascii=False), flush=True)
    for p in PROTOCOLS3:
        a = out["aggregates"][p]
        print(f"AGG3 {p}: bit={a['mean_bit_accuracy']:.4f} exact_snap={a['mean_exact_snapshot_fraction']:.4f} "
              f"modal={a['modal_exact_fraction']:.3f} ({a['modal_exact_count']}) nontarget_flip={a['mean_nontarget_modal_flip_fraction']:.4f}", flush=True)
    txt = json.dumps(out, indent=1, ensure_ascii=False, default=float)
    dst = ROOT / "results" / "exp35" / "universal_control35_v3.json"
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(txt, encoding="utf-8")
        print("WROTE", dst, flush=True)
    except Exception as ex:
        alt = Path.home() / "_exp35_out"
        alt.mkdir(parents=True, exist_ok=True)
        (alt / "universal_control35_v3.json").write_text(txt, encoding="utf-8")
        print(f"WRITE_FAIL {dst} :: {type(ex).__name__}: {ex} -> WROTE {alt / 'universal_control35_v3.json'}", flush=True)
    print(f"ELAPSED_S={time.time()-t0:.1f}", flush=True)
    return out


if __name__ == "__main__":
    main()
