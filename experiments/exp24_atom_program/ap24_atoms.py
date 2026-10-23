"""Exp24 (atoms): natural atoms trapped in a cell. The floating body of a 4 K cell is a box of frozen boron acceptors; each
atom is either neutral (bound hole, B0, bit 0) or ionized (B-, bit 1). The cell sets the hole quasi-Fermi level of the box
(body voltage V_b); the inversion layer above is a grounded conducting plane. Atom-resolved Coulomb-glass model:
  grand potential  Omega(n) = sum_i n_i (E_a0 + eV_b + D_i) + 1/2 sum_{i != j} n_i n_j W_ij,   n_i = 1: ionized
  E_a0 = -E_g + E_A = -1.17 + 0.045 eV (acceptor level relative to the conduction band at the inversion layer)
  W_ij = k (1/|r_i - r_j| - 1/|r_i - r_j*|)  (Coulomb + image in the conducting plane), D_i = -k / (4 z_i) (self-image),
  k = e^2 / (4 pi eps_Si) = 0.1231 eV nm; lateral periodic replicas for the 'device' box, none for an isolated island.
Zero-temperature relaxation with single flips and single-electron hops; V_b swept up and down in 0.25 mV steps. Output:
switching voltage of every atom (up and down sweep), events with exactly one atom, spacing statistics, hysteresis.
Continuum check: x_d(V_b) = sqrt(2 eps (1.125 - V_b) / (q N_A)), mean step e x_d / (eps A).
usage: python ap24_atoms.py predict | run [nproc] | analyze   ->  results/exp24/"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp24"
K = 1.4399645 / 11.7                     # eV nm
NA = 1.7e17 * 1e-21                      # nm^-3
EA0 = -1.17 + 0.045                      # eV
EPS = 11.7 * 8.8541878128e-12            # F/m
Q = 1.602176634e-19
DEPTH, ZMIN = 150.0, 2.0                 # nm
BOXES = [("island100", 100.0, False, range(5)), ("island200", 200.0, False, range(5)), ("device300", 300.0, True, range(2))]


def xd_nm(vb):
    return np.sqrt(np.maximum(2 * EPS * (1.125 - vb) / (Q * NA * 1e27), 0.0)) * 1e9


def build(side, periodic, seed):
    rng = np.random.default_rng(1000 + seed)
    n = rng.poisson(NA * side * side * (DEPTH - ZMIN))
    r = np.column_stack([rng.uniform(0, side, n), rng.uniform(0, side, n), rng.uniform(ZMIN, DEPTH, n)])
    W = np.zeros((n, n))
    D = -K / (4 * r[:, 2])
    shifts = [(a, b) for a in range(-2, 3) for b in range(-2, 3)] if periodic else [(0, 0)]
    for a, b in shifts:
        s = np.array([a * side, b * side, 0.0])
        d = r[:, None, :] - (r[None, :, :] + s)
        dist = np.sqrt((d**2).sum(-1))
        dimg = d.copy()
        dimg[..., 2] = r[:, None, 2] + r[None, :, 2]
        im = np.sqrt((dimg**2).sum(-1))
        if a == 0 and b == 0:
            np.fill_diagonal(dist, np.inf)
            W += K * (1 / dist - 1 / im) * (1 - np.eye(n)) - np.diag(np.zeros(n))
        else:
            term = K * (1 / dist - 1 / im)
            D += 0.5 * np.diag(term)            # own replicas (identical configuration in every cell of the lattice)
            np.fill_diagonal(term, 0.0)
            W += term
    return r, W, D


def relax(nocc, h, eps_site, W, max_iter=100000):
    """single flips and single-electron hops until no move lowers Omega; returns the list of flipped sites"""
    flipped = []
    for _ in range(max_iter):
        e = eps_site + h
        dflip = np.where(nocc == 1, -e, e)
        i = int(np.argmin(dflip))
        if dflip[i] < -1e-12:
            s = 1 if nocc[i] == 0 else -1
            nocc[i] += s
            h += s * W[:, i]
            flipped.append(i)
            continue
        occ, emp = np.nonzero(nocc == 1)[0], np.nonzero(nocc == 0)[0]
        if len(occ) == 0 or len(emp) == 0:
            break
        hop = (e[emp][None, :] - W[np.ix_(occ, emp)]) - e[occ][:, None]
        a, b = np.unravel_index(int(np.argmin(hop)), hop.shape)
        if hop[a, b] < -1e-12:
            i, j = occ[a], emp[b]
            nocc[i], nocc[j] = 0, 1
            h += W[:, j] - W[:, i]
            flipped += [i, j]
            continue
        break
    return flipped


def sweep(r, W, D, vbs, nocc=None, h=None):
    """V_b sweep; continues from (nocc, h) when given (down sweep after the up sweep -> hysteresis)"""
    n = len(r)
    if nocc is None:
        nocc, h = np.zeros(n, int), np.zeros(n)
    events = []
    relax(nocc, h, EA0 + vbs[0] + D, W)
    n0 = int(nocc.sum())
    for vb in vbs[1:]:
        f = relax(nocc, h, EA0 + vb + D, W)
        if f:
            events.append(dict(vb=float(vb), sites=[int(x) for x in f], n_ion=int(nocc.sum())))
    return events, nocc, h, n0


def run(job):
    t0 = time.time()
    name, side, periodic, seed = job
    r, W, D = build(side, periodic, seed)
    up = np.round(np.arange(0.0, 1.0501, 0.00025), 5)
    ev_up, nocc, h, n0 = sweep(r, W, D, up)
    ev_dn, _, _, n_top = sweep(r, W, D, up[::-1], nocc, h)
    res = dict(name=f"{name}_s{seed}", side=side, periodic=periodic, n_atoms=len(r), z=r[:, 2].tolist(),
               xy=r[:, :2].tolist(), n_ion_at_0=n0, n_ion_at_top=n_top, up=ev_up, down=ev_dn, wall=time.time() - t0)
    (OUT / "atom_runs").mkdir(parents=True, exist_ok=True)
    (OUT / "atom_runs" / f"{res['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    return f"{res['name']}: {len(r)} atoms, {len(ev_up)} up events, {len(ev_dn)} down events ({res['wall']:.0f} s)"


def predict():
    f = OUT / "predictions_atoms.json"
    if f.exists():
        print("already registered:", f)
        return
    rows = {}
    for name, side, periodic, _ in BOXES:
        xd = xd_nm(0.97)
        step = Q * xd * 1e-9 / (EPS * (side * 1e-9) ** 2)
        n_sw = NA * side * side * (xd_nm(0.0) - xd)
        rows[name] = dict(xd_097_nm=float(xd), mean_step_mV=step * 1e3, atoms_switched_0_to_097=float(n_sw),
                          poisson_isolated_fraction_4K=float(np.exp(-2 * 5 * 0.000345 / step)),
                          poisson_isolated_fraction_1K=float(np.exp(-2 * 5 * 0.0000862 / step)))
        print(name, {k: round(v, 4) for k, v in rows[name].items()})
    text = ["A1: atoms change state in depth order: Spearman rank correlation between switching voltage (up sweep) and "
            "depth >= 0.8 in every box",
            "A2: mean spacing of single-atom switching voltages for V_b in 0.90-1.00 V = e x_d/(eps A) within 30 % (islands)",
            "A3: number of atoms switched between V_b = 0 and 0.97 V = N_A A (x_d(0) - x_d(0.97)) within 15 %",
            "A4: fraction of isolated single-atom events (gap to both neighbours >= 5 kT/e) >= the Poisson estimate "
            "exp(-2 delta / step) at 4 K and at 1 K (Coulomb repulsion makes the staircase more regular)",
            "A5 (exploration): atoms with different up and down switching voltages (hysteresis = non-volatile 0/1)"]
    f.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), rows=rows, predictions=text), indent=1),
                 encoding="utf-8")


def analyze():
    from scipy import stats
    pr = json.loads((OUT / "predictions_atoms.json").read_text(encoding="utf-8"))
    out = {}
    for name, side, periodic, seeds in BOXES:
        agg = dict(rho=[], step=[], n_sw=[], iso4=[], iso1=[], hyst=[], single=[], n_atoms=[])
        for s in seeds:
            f = OUT / "atom_runs" / f"{name}_s{s}.json"
            if not f.exists():
                continue
            d = json.loads(f.read_text(encoding="utf-8"))
            z = np.array(d["z"])
            first_up = {}
            for e in d["up"]:
                for site in e["sites"]:
                    first_up.setdefault(site, e["vb"])
            last_dn = {}
            for e in d["down"]:
                for site in e["sites"]:
                    last_dn.setdefault(site, e["vb"])
            sites = sorted(first_up)
            if len(sites) > 3:
                agg["rho"].append(float(stats.spearmanr([first_up[k] for k in sites], z[sites])[0]))
            single = [e for e in d["up"] if len(e["sites"]) == 1]
            agg["single"].append(len(single) / max(len(d["up"]), 1))
            v = np.array([e["vb"] for e in d["up"]])
            sel = v[(v >= 0.90) & (v <= 1.0)]
            if len(sel) > 2:
                agg["step"].append(float(np.mean(np.diff(sel))))
            ion_hi = [e["n_ion"] for e in d["up"] if e["vb"] <= 0.97]
            agg["n_sw"].append(float(d["n_ion_at_0"] - (ion_hi[-1] if ion_hi else d["n_ion_at_0"])))
            vs = np.array([e["vb"] for e in d["up"]])
            for T, key in ((4.0, "iso4"), (1.0, "iso1")):
                dlt = 5 * 8.617e-5 * T
                ok = 0
                for k, e in enumerate(d["up"]):
                    if len(e["sites"]) != 1:
                        continue
                    gp = vs[k] - vs[k - 1] if k > 0 else np.inf
                    gn = vs[k + 1] - vs[k] if k + 1 < len(vs) else np.inf
                    ok += (gp >= dlt) and (gn >= dlt)
                agg[key].append(ok / max(len(d["up"]), 1))
            hy = [first_up[k] - last_dn[k] for k in sites if k in last_dn]
            agg["hyst"].append(float(np.mean(np.array(hy) > 0.00025)) if hy else 0.0)
            agg["n_atoms"].append(d["n_atoms"])
        p = pr["rows"][name]
        out[name] = {k: (float(np.mean(v)) if v else None) for k, v in agg.items()}
        out[name].update(pred=p, A1=bool(agg["rho"]) and min(agg["rho"]) >= 0.8,
                         A2=(abs(out[name]["step"] / (p["mean_step_mV"] / 1e3) - 1) <= 0.3) if (out[name]["step"] and not periodic) else None,
                         A3=abs(out[name]["n_sw"] / p["atoms_switched_0_to_097"] - 1) <= 0.15 if out[name]["n_sw"] else None,
                         A4=(out[name]["iso4"] >= p["poisson_isolated_fraction_4K"] and out[name]["iso1"] >= p["poisson_isolated_fraction_1K"]))
        print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[name].items() if k != "pred"})
    (OUT / "analysis_atoms.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def metropolis(nocc, h, eps_site, W, kT, sweeps, rng, k):
    """grand-canonical single flips + canonical hops at temperature kT (eV); returns the mean occupation of site k"""
    n = len(nocc)
    acc_k, cnt = 0.0, 0
    for s in range(sweeps):
        for _ in range(n):
            if rng.random() < 0.5:
                i = int(rng.integers(n))
                dO = (1 - 2 * nocc[i]) * (eps_site[i] + h[i])
                if dO <= 0 or rng.random() < np.exp(-dO / kT):
                    sgn = 1 - 2 * nocc[i]
                    nocc[i] += sgn
                    h += sgn * W[:, i]
            else:
                occ, emp = np.nonzero(nocc == 1)[0], np.nonzero(nocc == 0)[0]
                if len(occ) == 0 or len(emp) == 0:
                    continue
                i, j = int(rng.choice(occ)), int(rng.choice(emp))
                dO = (eps_site[j] + h[j] - W[i, j]) - (eps_site[i] + h[i])
                if dO <= 0 or rng.random() < np.exp(-dO / kT):
                    nocc[i], nocc[j] = 0, 1
                    h += W[:, j] - W[:, i]
        if s >= sweeps // 4:
            acc_k += nocc[k]
            cnt += 1
    return acc_k / cnt


def thermal_job(args):
    seed, T = args
    r, W, D = build(100.0, False, seed)
    d = json.loads((OUT / "atom_runs" / f"island100_s{seed}.json").read_text(encoding="utf-8"))
    vs = np.array([e["vb"] for e in d["up"]])
    dl = 5 * 8.617e-5 * T
    picks = []
    for i, e in enumerate(d["up"]):
        if len(e["sites"]) != 1:
            continue
        gp = vs[i] - vs[i - 1] if i else np.inf
        gn = vs[i + 1] - vs[i] if i + 1 < len(vs) else np.inf
        if gp >= dl and gn >= dl and np.isfinite(gp) and np.isfinite(gn):
            picks.append((e["sites"][0], vs[i], gp, gn))
    rng = np.random.default_rng(seed)
    sel = [picks[i] for i in sorted(rng.choice(len(picks), min(8, len(picks)), replace=False))] if picks else []
    out = []
    for k, vk, gp, gn in sel:
        row = dict(site=int(k), V_k=float(vk))
        for label, vb in (("0", vk - gp / 2), ("1", vk + gn / 2)):    # 0: atom ionized (B-), 1: neutral (B0)
            nocc, h = np.zeros(len(r), int), np.zeros(len(r))
            relax(nocc, h, EA0 + vb + D, W)
            p_ion = metropolis(nocc, h, EA0 + vb + D, W, 8.617e-5 * T, 1500, rng, k)
            row[f"P_target_{label}"] = float(p_ion if label == "0" else 1 - p_ion)
        out.append(row)
    return dict(seed=seed, T=T, n_candidates=len(picks), rows=out)


def thermal():
    pf = OUT / "predictions_thermal.json"
    if not pf.exists():
        pf.write_text(json.dumps(dict(
            registered=time.strftime("%Y-%m-%d %H:%M"),
            predictions=["TH1: for atoms classified addressable at 4 K (single-atom event isolated by >= 5 kT/e), setting V_b to "
                         "the middle of the gap below / above the threshold gives the target state (B- = 0 / B0 = 1) with mean "
                         "probability >= 0.9 in equilibrium at 4 K (Metropolis, 100 nm islands, 5 boxes, up to 8 atoms each)",
                         "TH2: the same at 1 K with the 1 K criterion gives mean probability >= 0.9"]), indent=1), encoding="utf-8")
    jobs = [(s, T) for T in (4.0, 1.0) for s in range(5)]
    with Pool(5) as pool:
        res = pool.map(thermal_job, jobs)
    summ = {}
    for T in (4.0, 1.0):
        rows = [r_ for x in res if x["T"] == T for r_ in x["rows"]]
        p = [r_[f"P_target_{b}"] for r_ in rows for b in ("0", "1")]
        summ[f"T={T:g} K"] = dict(atoms=len(rows), mean_P=float(np.mean(p)) if p else None, min_P=float(np.min(p)) if p else None,
                                  frac_ge_0p9=float(np.mean(np.array(p) >= 0.9)) if p else None)
        print(f"T = {T:g} K: {len(rows)} atoms, mean P(target) {summ[f'T={T:g} K']['mean_P']:.3f}, "
              f"min {summ[f'T={T:g} K']['min_P']:.3f}, fraction >= 0.9: {summ[f'T={T:g} K']['frac_ge_0p9']:.2f}")
    (OUT / "thermal.json").write_text(json.dumps(dict(summary=summ, runs=res), indent=1), encoding="utf-8")


def fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, axs = plt.subplots(2, 2, figsize=(14, 10))
    ax = axs[0, 0]
    ps = json.loads((OUT / "predictions_scale.json").read_text(encoding="utf-8"))["rows"]
    an = json.loads((OUT / "analysis_scale.json").read_text(encoding="utf-8"))
    for a_mV, col in ((30, "C0"), (150, "C3")):
        keys = [k for k in ps if ps[k]["a_mV"] == a_mV]
        lc = [ps[k]["lambda_c"] for k in keys]
        ax.plot(lc, [10 * ps[k]["P_prog"] for k in keys], "--o", c=col, mfc="none", label=f"öngörü (bağımsız birim), tohum {a_mV} mV")
        ax.plot(lc, [an[k]["births_on_seed"] for k in keys], "s", c=col, ms=9, label=f"simülasyon, tohum {a_mV} mV")
    ax.axvline(2.2, c="k", ls=":", lw=1)
    ax.text(2.3, 5.2, "gövde bağlanma\nuzunluğu ℓ_b ≈ 2,2 µm", fontsize=8)
    ax.axvline(0.083, c="0.5", ls=":", lw=1)
    ax.text(0.086, 8.0, "W_dep ≈ 83 nm", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("λ_c (µm)  ~ hücre boyu / 2")
    ax.set_ylabel("ilk hücrenin tohumda doğduğu alan (10 üzerinden)")
    ax.set_title("(a) hücre küçüldükçe programlanabilirlik: ℓ_b altında tohum gövdede yayılıp kayboluyor", fontsize=9)
    ax.legend(fontsize=7)
    d = json.loads((OUT / "atom_runs" / "island100_s0.json").read_text(encoding="utf-8"))
    ax = axs[0, 1]
    for key, col, lab in (("up", "C0", "gövde gerilimi artarken (nötrleşme)"), ("down", "C3", "azalırken (iyonlaşma)")):
        v = [e["vb"] for e in d[key]]
        n = [d["n_atoms"] - e["n_ion"] for e in d[key]]
        ax.step(v, n, where="post", c=col, lw=1.2, label=lab)
    ax.set_xlim(0.65, 1.06)
    ax.set_xlabel("V_b (V)")
    ax.set_ylabel("nötr (deliği bağlı) bor atomu sayısı")
    ax.set_title("(b) 100 nm'lik hücre gövdesi: atomlar tek tek 0 → 1; iki yönde farklı eşik = histerezis", fontsize=9)
    ax.legend(fontsize=8)
    ax = axs[1, 0]
    for name, col in (("island100_s0", "C0"), ("device300_s0", "C2")):
        p = OUT / "atom_runs" / f"{name}.json"
        if not p.exists():
            continue
        dd = json.loads(p.read_text(encoding="utf-8"))
        z = np.array(dd["z"])
        first = {}
        for e in dd["up"]:
            for s_ in e["sites"]:
                first.setdefault(s_, e["vb"])
        ks = sorted(first)
        ax.plot(z[ks], [first[k] for k in ks], ".", c=col, ms=4, label=name)
    ax.set_xlabel("atomun derinliği (nm)")
    ax.set_ylabel("nötrleştiği gövde gerilimi (V)")
    ax.set_title("(c) her atomun kendi 'adresi': eşik gerilimi derinliğe göre sıralı", fontsize=9)
    ax.legend(fontsize=8)
    ax = axs[1, 1]
    rows = [("T1 hücre gövdesinde tek akseptörün yükü", "✓", "✓", "~", "~"), ("T7 kenar tünelleme kabuğu (~1 atom)", "✓", "✓", "~", "~"),
            ("T3 tek ³¹P çekirdek spini (+RF)", "✓", "✓", "✓*", "✓*"), ("T4 tek atomluk metal temas", "✓", "✓", "✓", "✓"),
            ("T5 ReRAM O-boşluğu zinciri (~1 G₀)", "✓", "✓", "✓", "✓"), ("T8 akseptör spini (EDSR + hücre adresi)", "✓", "~", "~", "✓*"),
            ("T2 cephe alanıyla iyonlaştırma", "✓", "~", "✗", "~"), ("T6 Si–H bağı (titreşim merdiveni)", "✓", "✗", "✗", "~"),
            ("hücrenin kendisi atom boyuna", "✗", "✗", "✗", "✗"), ("örgü atomunu doğrudan oynatma", "✗", "✗", "✗", "✗")]
    colors = {"✓": "#7fc97f", "~": "#fdc086", "✗": "#f0027f", "✓*": "#beaed4"}
    for i, r_ in enumerate(rows):
        ax.text(-0.02, len(rows) - i - 0.5, r_[0], ha="right", va="center", fontsize=8)
        for j, v_ in enumerate(r_[1:]):
            ax.add_patch(plt.Rectangle((j, len(rows) - i - 1), 0.95, 0.9, color=colors[v_]))
            ax.text(j + 0.47, len(rows) - i - 0.55, v_, ha="center", va="center", fontsize=9)
    ax.set_xlim(-3.6, 4)
    ax.set_ylim(0, len(rows) + 0.6)
    for j, t in enumerate(("iki durum", "hücreyle\nkontrol", "tek atom\nadres", "okuma")):
        ax.text(j + 0.47, len(rows) + 0.25, t, ha="center", fontsize=8)
    ax.axis("off")
    ax.set_title("(d) teknikler ve dört ölçüt (✓* = ek RF/SET altyapısıyla)", fontsize=9)
    f.suptitle("Exp24: doğal atomları hücrelerle 0/1 programlamak — boyutla değil gerilim/frekans/akım adresiyle", fontsize=12)
    f.tight_layout(rect=(0, 0, 1, 0.95))
    f.savefig(OUT / "fig_exp24.png", dpi=95)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[1] == "fig":
        fig()
    elif sys.argv[1] == "thermal":
        thermal()
    elif sys.argv[1] == "predict":
        predict()
    elif sys.argv[1] == "analyze":
        analyze()
    elif sys.argv[1] == "run":
        js = [(n, s, p, k) for n, s, p, ks in BOXES for k in ks
              if not (OUT / "atom_runs" / f"{n}_s{k}.json").exists()]
        with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 3) as pool:
            for line in pool.imap_unordered(run, js):
                print(line, flush=True)
        print("TAMAM", flush=True)
