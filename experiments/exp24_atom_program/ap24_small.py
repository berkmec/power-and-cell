"""Exp24 (small bodies, 4 K): test of the requirement derived from the 100 nm thermal test. Mean step <dV> ~ e x_d/(eps A)
grows as 1/A, so a smaller cell body gives larger gaps between the atoms' addresses; and the address of each atom is
calibrated in equilibrium at the operating temperature (finite-T sweep) instead of being read from a zero-temperature sweep.
Method: islands of 40 and 60 nm (10 each, same Coulomb-glass model as ap24_atoms.py), Metropolis sweep of V_b at 4 K
(0.5 mV steps, state carried over) -> <n_i>(V_b); address V_i* = V_b where <n_i> crosses 0.5. Atoms whose calibrated address
is separated from the neighbouring addresses by >= 2 x 8 kT/e are programmed at the middle of the gap below (bit 0 = B-) and
above (bit 1 = B0); the target-state probability is measured with a fresh Metropolis run from the zero-T state.
usage: python ap24_small.py run  ->  results/exp24/small_bodies.json (predictions written first: predictions_small.json)"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ap24_atoms import build, relax, metropolis, EA0, OUT  # noqa: E402

KT4 = 8.617e-5 * 4.0
WINDOW = (0.805, 1.006)          # the cell's own stable high-state body voltage (V_D from V_down to 1.2 V)


def mc_mean(nocc, h, eps, W, kT, sweeps, rng):
    """Metropolis that returns the mean occupation of every site (state carried in place)"""
    n = len(nocc)
    acc = np.zeros(n)
    cnt = 0
    for s in range(sweeps):
        for _ in range(n):
            if rng.random() < 0.5:
                i = int(rng.integers(n))
                dO = (1 - 2 * nocc[i]) * (eps[i] + h[i])
                if dO <= 0 or rng.random() < np.exp(-dO / kT):
                    sg = 1 - 2 * nocc[i]
                    nocc[i] += sg
                    h += sg * W[:, i]
            else:
                occ, emp = np.nonzero(nocc == 1)[0], np.nonzero(nocc == 0)[0]
                if len(occ) and len(emp):
                    i, j = int(rng.choice(occ)), int(rng.choice(emp))
                    dO = (eps[j] + h[j] - W[i, j]) - (eps[i] + h[i])
                    if dO <= 0 or rng.random() < np.exp(-dO / kT):
                        nocc[i], nocc[j] = 0, 1
                        h += W[:, j] - W[:, i]
        if s >= sweeps // 3:
            acc += nocc
            cnt += 1
    return acc / cnt


def job(args):
    side, seed = args
    t0 = time.time()
    r, W, D = build(side, False, seed)
    n = len(r)
    rng = np.random.default_rng(500 + seed)
    vbs = np.round(np.arange(0.55, 1.10001, 0.0005), 5)
    nocc, h = np.zeros(n, int), np.zeros(n)
    relax(nocc, h, EA0 + vbs[0] + D, W)
    prof = []
    for vb in vbs:
        prof.append(mc_mean(nocc, h, EA0 + vb + D, W, KT4, 60, rng))
    prof = np.array(prof)                                   # (V_b, site) mean ionised occupation
    addr = {}
    for k in range(n):
        x = prof[:, k]
        cross = np.nonzero((x[:-1] >= 0.5) & (x[1:] < 0.5))[0]
        if len(cross) == 1 and x[0] > 0.9 and x[-1] < 0.1:
            i = cross[0]
            addr[k] = float(vbs[i] + (x[i] - 0.5) / (x[i] - x[i + 1]) * (vbs[i + 1] - vbs[i]))
    ks = sorted(addr, key=addr.get)
    va = np.array([addr[k] for k in ks])
    need = 2 * 8 * KT4
    rows = []
    for idx, k in enumerate(ks):
        gp = va[idx] - va[idx - 1] if idx > 0 else np.inf
        gn = va[idx + 1] - va[idx] if idx + 1 < len(va) else np.inf
        if gp >= need and gn >= need and np.isfinite(gp) and np.isfinite(gn):
            row = dict(site=int(k), address=float(va[idx]), in_cell_window=bool(WINDOW[0] <= va[idx] <= WINDOW[1]))
            for label, vb in (("0", va[idx] - gp / 2), ("1", va[idx] + gn / 2)):
                nn, hh = np.zeros(n, int), np.zeros(n)
                relax(nn, hh, EA0 + vb + D, W)
                p_ion = metropolis(nn, hh, EA0 + vb + D, W, KT4, 1500, rng, k)
                row[f"P_target_{label}"] = float(p_ion if label == "0" else 1 - p_ion)
            rows.append(row)
    return dict(side=side, seed=seed, n_atoms=n, n_addressed=len(addr), n_qualified=len(rows), rows=rows,
                wall=time.time() - t0)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pf = OUT / "predictions_small.json"
    if not pf.exists():
        pf.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), predictions=[
            "SB1: with equilibrium-calibrated addresses and gaps >= 2 x 8 kT/e, the mean target-state probability at 4 K is "
            ">= 0.9 and >= 80 % of the (atom, state) pairs are >= 0.9 (40 and 60 nm islands together)",
            "SB2: qualified atoms per island inside the cell's own window (0.805-1.006 V) >= 2 for 40 nm islands (mean)",
            "SB3: qualified atoms per island grow when the body shrinks from 60 to 40 nm relative to the number of atoms "
            "(fraction of atoms qualified larger at 40 nm)"]), indent=1), encoding="utf-8")
    jobs = [(s, k) for s in (40.0, 60.0) for k in range(10)]
    with Pool(8) as pool:
        res = pool.map(job, jobs)
    summ = {}
    for side in (40.0, 60.0):
        rr = [x for x in res if x["side"] == side]
        p = [row[f"P_target_{b}"] for x in rr for row in x["rows"] for b in ("0", "1")]
        win = [sum(row["in_cell_window"] for row in x["rows"]) for x in rr]
        summ[f"{side:g} nm"] = dict(atoms=float(np.mean([x["n_atoms"] for x in rr])),
                                    addressed=float(np.mean([x["n_addressed"] for x in rr])),
                                    qualified=float(np.mean([x["n_qualified"] for x in rr])),
                                    qualified_in_window=float(np.mean(win)), mean_P=float(np.mean(p)) if p else None,
                                    frac_P_ge_0p9=float(np.mean(np.array(p) >= 0.9)) if p else None, n_pairs=len(p))
        print(f"{side:g} nm: " + ", ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}" for k, v in summ[f'{side:g} nm'].items()))
    (OUT / "small_bodies.json").write_text(json.dumps(dict(summary=summ, runs=res), indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
