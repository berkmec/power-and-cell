"""Exp25 (field): universal logic from field-coupled natural atoms in a cell body (4 K).
Every ionised acceptor acts on every other atom through its Coulomb field (W_ij = k(1/r - 1/r*), exp24 model). For an output
atom o and two input atoms i, j (forced by local gates, +-80 mV; an assumption: SKY130 cannot draw ~20 nm gates):
  o is ionised (bit 1) iff  V_b < V_o - W_oi n_i - W_oj n_j      (n = 1: ionised)
so a body voltage inside (V_o - W_oi - W_oj, V_o - max(W_oi, W_oj)) gives n_o = NAND(n_i, n_j); NAND is a universal gate.
Method: 60 nm islands (10), zero-T sweep -> isolated single-atom outputs (gaps >= 2 x 8 kT/e); the two atoms with the largest
coupling to the output are the inputs; the body voltage is scanned across the predicted window; each truth-table row is
measured with Metropolis at 4 K (background atoms free). Also: reach of one atom's field, |W_oj| > k T for all atoms within
30 nm (the programmable atom interacts with all its neighbours).
usage: python ap25_logic.py  ->  results/exp25/logic.json (predictions_logic.json written first)"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp24_atom_program"))
from ap24_atoms import build, relax, sweep, metropolis, EA0  # noqa: E402

OUT = ROOT / "results" / "exp25"
KT = 8.617333e-5 * 4.0
FORCE = 0.08


def row_prob(r, W, D, vb, i, j, o, ni, nj, rng):
    off = np.zeros(len(r))
    off[i] = -FORCE if ni else FORCE
    off[j] = -FORCE if nj else FORCE
    eps = EA0 + vb + D + off
    nocc, h = np.zeros(len(r), int), np.zeros(len(r))
    relax(nocc, h, eps, W)
    p_ion = metropolis(nocc, h, eps, W, KT, 600, rng, o)
    target = 0 if (ni and nj) else 1
    return float(p_ion if target == 1 else 1 - p_ion)


def job(seed):
    t0 = time.time()
    r, W, D = build(60.0, False, seed)
    rng = np.random.default_rng(900 + seed)
    vbs = np.round(np.arange(0.55, 1.12, 0.00025), 5)
    ev, _, _, _ = sweep(r, W, D, vbs)
    v = np.array([e["vb"] for e in ev])
    need = 2 * 8 * KT
    outs = []
    for k, e in enumerate(ev):
        if len(e["sites"]) == 1 and 0 < k < len(ev) - 1 and v[k] - v[k - 1] >= need and v[k + 1] - v[k] >= need:
            outs.append((e["sites"][0], float(v[k])))
    dist = np.sqrt(((r[:, None, :] - r[None, :, :]) ** 2).sum(-1))
    reach = []
    for o, _ in outs:
        near = (dist[o] < 30.0) & (dist[o] > 0)
        reach.append(float(np.mean(np.abs(W[o, near]) > KT)) if near.any() else None)
    tried = []
    for o, vo in outs[:8]:
        order = np.argsort(-W[o])
        i, j = [int(x) for x in order[:2]]
        wi, wj = W[o, i], W[o, j]
        lo, hi = vo - wi - wj, vo - max(wi, wj)
        best = None
        for vb in np.linspace(lo, hi, 7)[1:-1]:
            rows = [row_prob(r, W, D, vb, i, j, o, a, b, rng) for a in (0, 1) for b in (0, 1)]
            if best is None or min(rows) > min(best["rows"]):
                best = dict(vb=float(vb), rows=rows)
        tried.append(dict(output=int(o), inputs=[i, j], W_meV=[float(wi * 1e3), float(wj * 1e3)],
                          min_W_over_kT=float(min(wi, wj) / KT), V_o=vo, best=best,
                          works=bool(best and min(best["rows"]) >= 0.9)))
    return dict(seed=seed, n_atoms=len(r), n_outputs=len(outs), reach_fraction=reach, triples=tried, wall=time.time() - t0)


def thresholds(r, W, D, i, j, o, v0):
    """zero-T threshold of the output for each forced input combination (output ionised iff V_b < V_th)"""
    out = {}
    vbs = np.round(np.arange(v0 - 0.12, v0 + 0.12, 0.00025), 5)
    for ni in (0, 1):
        for nj in (0, 1):
            off = np.zeros(len(r))
            off[i] = -FORCE if ni else FORCE
            off[j] = -FORCE if nj else FORCE
            nocc, h = np.zeros(len(r), int), np.zeros(len(r))
            vth = None
            for vb in vbs:
                relax(nocc, h, EA0 + vb + D + off, W)
                if nocc[o] == 0:
                    vth = float(vb)
                    break
            out[(ni, nj)] = vth
    return out


def job_cal(seed):
    t0 = time.time()
    r, W, D = build(60.0, False, seed)
    rng = np.random.default_rng(1900 + seed)
    vbs = np.round(np.arange(0.55, 1.12, 0.00025), 5)
    ev, _, _, _ = sweep(r, W, D, vbs)
    v = np.array([e["vb"] for e in ev])
    need = 2 * 8 * KT
    outs = [(e["sites"][0], float(v[k])) for k, e in enumerate(ev)
            if len(e["sites"]) == 1 and 0 < k < len(ev) - 1 and v[k] - v[k - 1] >= need and v[k + 1] - v[k] >= need]
    tried = []
    for o, vo in outs[:8]:
        order = np.argsort(-W[o])
        i, j = [int(x) for x in order[:2]]
        th = thresholds(r, W, D, i, j, o, vo)
        if any(x is None for x in th.values()):
            tried.append(dict(output=int(o), inputs=[i, j], thresholds=None, works=False))
            continue
        lo, hi = th[(1, 1)], min(th[(0, 0)], th[(0, 1)], th[(1, 0)])
        row = dict(output=int(o), inputs=[i, j], thresholds={f"{a}{b}": th[(a, b)] for a, b in th}, window_mV=(hi - lo) * 1e3)
        if hi - lo > 0:
            vb = 0.5 * (lo + hi)
            rows = [row_prob(r, W, D, vb, i, j, o, a, b, rng) for a in (0, 1) for b in (0, 1)]
            row.update(vb=vb, rows=rows, works=bool(min(rows) >= 0.9))
        else:
            row.update(works=False)
        tried.append(row)
    return dict(seed=seed, triples=tried, wall=time.time() - t0)


def calibrated():
    pf = OUT / "predictions_logic2.json"
    if not pf.exists():
        pf.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), basis=(
            "LG1 failed (3/10 islands): in 23 of 25 failing triples only the (1,1)->0 row failed, because the window assumed "
            "neutral inputs at the reference threshold; here the output threshold is calibrated for each input combination"),
            predictions=[
                "LG4: with input-resolved threshold calibration >= 50 % of the 60 nm islands host >= 1 NAND triple whose four "
                "rows are all >= 0.9 at 4 K",
                "LG5: among triples whose calibrated NAND window is >= 2 x 8 kT/e (5.5 mV) wide, >= 70 % pass all four rows"]),
            indent=1), encoding="utf-8")
    with Pool(5) as pool:
        res = pool.map(job_cal, range(10))
    tri = [t for x in res for t in x["triples"]]
    wide = [t for t in tri if t.get("window_mV") is not None and t["window_mV"] >= 2 * 8 * KT * 1e3]
    summ = dict(islands_with_NAND=sum(1 for x in res if any(t["works"] for t in x["triples"])), islands=len(res),
                triples=len(tri), working=sum(t["works"] for t in tri), with_window=sum(1 for t in tri if (t.get("window_mV") or 0) > 0),
                wide_windows=len(wide), wide_working=sum(t["works"] for t in wide))
    print(json.dumps(summ, indent=1))
    (OUT / "logic_calibrated.json").write_text(json.dumps(dict(summary=summ, runs=res), indent=1), encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pf = OUT / "predictions_logic.json"
    if not pf.exists():
        pf.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), predictions=[
            "LG1: >= 50 % of the 60 nm islands host at least one NAND triple whose four truth-table rows are all >= 0.9 at 4 K",
            "LG2: every working triple has min(W_oi, W_oj) >= 4 kT; no triple with min W < 2 kT works",
            "LG3: the field of an output atom exceeds k T at >= 90 % of the atoms within 30 nm (it reaches all its neighbours)"]),
            indent=1), encoding="utf-8")
    with Pool(5) as pool:
        res = pool.map(job, range(10))
    tri = [t for x in res for t in x["triples"]]
    ok = [t for t in tri if t["works"]]
    islands_ok = sum(1 for x in res if any(t["works"] for t in x["triples"]))
    reach = [v for x in res for v in x["reach_fraction"] if v is not None]
    summ = dict(islands=len(res), islands_with_NAND=islands_ok, triples=len(tri), working=len(ok),
                min_W_over_kT_working=[t["min_W_over_kT"] for t in ok],
                min_W_over_kT_failing=[t["min_W_over_kT"] for t in tri if not t["works"]],
                mean_reach_fraction=float(np.mean(reach)) if reach else None)
    print(json.dumps({k: (v if not isinstance(v, list) else [round(x, 1) for x in v][:20]) for k, v in summ.items()}, indent=1))
    (OUT / "logic.json").write_text(json.dumps(dict(summary=summ, runs=res), indent=1), encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "calibrated":
        calibrated()
    else:
        main()
