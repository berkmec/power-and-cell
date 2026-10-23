"""Exp26 (natural-silicon register): can the programmable atom reach ALL 29Si within 6 nm of NATURAL silicon?
Exp25: the hyperfine key alone fails in natural Si (~2119 nuclei crowd the spectrum, the electron T2 is too short).
Protocol tested here (each piece tied to a measured number; extrapolations are labelled):
 P1  programming in D+: the cell ionizes the donor (Exp24/25 charge control) -> no electron T1/T2 limit and no hyperfine
     crowding while the nuclei are written; back in D0 the electron couples to every nucleus again (A_k, couplings_dir).
 P2  MREV-16 homonuclear decoupling (Ladd et al., PRB 71, 014401 (2005), natural Si single crystal): T2 = 25.0 s at cycle
     t_c = 1.568 ms, RF ~120 Hz off resonance, T2 ~ t_c^-2 (1.3-4.7 ms), offsets scaled by s = 1/3, residual cross term
     ~ t_c^2 |D| w^2.  Model: T2(t_c, w) = min(cap, 25 s (1.568 ms / t_c)^2 (120 Hz / max(w, 120 Hz))^2);
     t_c < 1.3 ms and the w^2 law are extrapolations.
 P3  steerable field gradient (three on-chip currents): offset f_k = gamma G (R_k . g). A selective pulse of length
     tau = T2/10 resolves scaled offsets 4/tau -> W = 4/(s tau) = 120 Hz s / T2. Nucleus k is individually addressable if
     some direction g separates it from every other nucleus by d* = W / (gamma G) (in nm).
 P4  bus: 29Si-29Si dipolar network, recoupled ZZ, J = 4.748 Hz nm^3 (1 - 3 cos^2 theta) / r^3.
Derived laws: (L1) hyperfine key alone, Poisson spectrum: N_res = 0.63 sqrt(K/W) (one electron branch), 0.31 sqrt(K/W)
 (both branches), K = 3.2245 n_29 C_d = 1.19e5 Hz;  (L2) one direction: P(k isolated) = prod_j (1 - d*/r_kj)
 ~ exp(-d* S_k), S_k = sum_j 1/r_kj, because g.u is uniform on [-1, 1] for a random direction.
usage: python nr26.py law | run | bus | pulse | fixed | measured -> results/exp26/*.json (predictions written first);
figure: python nr26_fig.py"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1] / "exp25_universal_atom"))
import ap25_register as reg                                       # noqa: E402
from ap25_control import couplings_dir, lie_rank                  # noqa: E402

OUT = HERE.parents[2] / "results" / "exp26"
GN = reg.GN                                                       # 29Si |gamma|/2pi, Hz/T
C_D = 1e-7 * 2.0023 * 9.2740100783e-24 * GN * 1e27                # electron-29Si dipolar, Hz nm^3
D_NN = 1e-7 * GN**2 * 6.62607015e-34 * 1e27                       # 29Si-29Si secular ZZ, Hz nm^3 (4.748)
N29 = reg.ABUND * 8 / reg.A0_NM**3                                 # nm^-3
K_LAW = 8 / (3 * np.sqrt(3)) * 2 * np.pi / 3 * N29 * C_D           # Hz
S_MREV, T2_L, TC_L, W_L = 1 / 3, 25.0, 1.568, 120.0                # Ladd 2005 (s, ms, Hz)
T1E_DEV = 6.0                                                      # s, single donor at 1.5 T (Morello 2010)
R_MAX, N_CFG, N_DIR = 6.0, 20, 20000
BDIR = np.array([1.0, 2.0, 3.0]) / np.sqrt(14)
D_LINE, Z0 = 100e-9, 50.0                                          # drive line distance (m), impedance (ohm)


def register():
    OUT.mkdir(parents=True, exist_ok=True)
    pf = OUT / "predictions_nr26.json"
    if pf.exists():
        return
    pf.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), basis=(
        "Exp25: hyperfine key fails in natural Si; analytic laws L1 (N = 0.63 sqrt(K/W)) and L2 (exp(-d* S_k), S_k ~ "
        "350-530 nm^-1 in a 6 nm sphere) derived before any run of this script"), predictions=[
        "NR1: hyperfine key alone (D0, B || [1,2,3], electron state known): the number of individually addressable 29Si "
        "within 6 nm follows 0.63 sqrt(K/W) (K = 1.19e5 Hz) within +-40 % for W = 0.3, 1, 3, 10, 30 Hz (mean of 20 configs)",
        "NR2: with the measured device electron T1 (6 s at 1.5 T; pulse <= T1/20) and MREV-16 the hyperfine key reaches "
        "< 5 % of the 29Si within 6 nm",
        "NR3: D+ mode, steerable gradient: in >= 90 % of 20 natural configurations >= 99 % of the 29Si within 6 nm have a "
        "direction that isolates them by d* >= 0.010 nm; one fixed direction isolates < 5 % at d* = 0.010 nm",
        "NR4: at the measured MREV point (t_c = 1.568 ms, T2 <= 25 s) no gradient strength gives coverage >= 1 %",
        "NR5: full coverage (>= 99 % of the nuclei in >= 90 % of the configurations) needs t_c <= 0.4 ms if T2 follows the "
        "measured t_c^-2 law with no other cap, t_c <= 0.1 ms if T2 is capped at the measured 25 s, t_c <= 0.03 ms with an "
        "8 s cap (1/f noise, Ladd's polycrystal)",
        "NR6: bus: the 29Si dipolar graph with |J| >= 5/T2 (T2 = 25 s: 0.2 Hz) is connected in 20/20 configurations; "
        "4-nucleus J clusters with distinct offsets have full Lie rank 255 with selective X/Y (10/10) and with only global "
        "X/Y (>= 8/10)",
        "NR7: energy: at the largest feasible t_c (t_c^-2 law) the MREV drive at 100 nm needs >= 0.1 mW on average (Ladd's "
        "point: ~15 uW) and writing every nucleus costs >= 1 uJ per atom"]), indent=1), encoding="utf-8")


def configs(seed):
    rng = np.random.default_rng(seed)
    for _ in range(N_CFG):
        occ = rng.random(len(reg.SITES)) < reg.ABUND
        yield reg.SITES[occ], reg.A_SITES[occ]


def t2_model(tc_ms, w_hz, cap):
    return np.minimum(cap, T2_L * (TC_L / tc_ms) ** 2 * (W_L / np.maximum(w_hz, W_L)) ** 2)


def drive(tc_ms):
    """MREV-16 drive scaled from Ladd's point (pi/2 = 15 us at t_c = 1.568 ms; duty 16 pulses): average power into a 50 ohm
    line whose wire is D_LINE from the donor (B1 = mu0 I / (2 pi d))"""
    p = tc_ms * 1e-3 / (TC_L / 15e-3)
    b1_lin = 2 * (1 / (4 * p)) / GN
    i = 2 * np.pi * D_LINE * b1_lin / (4e-7 * np.pi)
    return dict(pi2_us=p * 1e6, b1_mT=b1_lin * 1e3, current_mA=i * 1e3, p_avg_W=0.5 * i**2 * Z0 * 16 * p / (tc_ms * 1e-3))


def lines_isolation(f, owner=None):
    """distance from each line to the nearest line of another nucleus (sorted spectrum)"""
    o = np.argsort(f)
    fs = f[o]
    n = len(f)
    ws = (np.arange(n) if owner is None else owner)[o]
    gap1, gap2 = np.full(n + 1, np.inf), np.full(n + 2, np.inf)
    gap1[1:n] = fs[1:] - fs[:-1]                                   # gap1[i] = fs[i] - fs[i-1]
    gap2[2:n] = fs[2:] - fs[:-2]                                   # gap2[i] = fs[i] - fs[i-2]
    same = np.zeros(n + 1, bool)
    same[1:n] = ws[1:] == ws[:-1]                                  # a nucleus' own second line is skipped (<= 2 lines each)
    left = np.where(same[:n], gap2[:n], gap1[:n])
    right = np.where(same[1:n + 1], gap2[2:n + 2], gap1[1:n + 1])
    out = np.empty(n)
    out[o] = np.minimum(left, right)
    return out


def law():
    """Part A: hyperfine key alone (D0), law L1, and the device-limited resolution"""
    register()
    W = np.logspace(-1.5, 2.5, 41)
    nu = GN * 1.0
    single, double, n_tot, nodec = [], [], [], []
    for P, Aiso in configs(2601):
        A = couplings_dir(P, Aiso, BDIR)
        n = len(A)
        iso1 = lines_isolation(np.abs(nu + A / 2))
        f2 = np.concatenate([np.abs(nu + A / 2), np.abs(nu - A / 2)])
        iso2 = lines_isolation(f2, np.concatenate([np.arange(n), np.arange(n)]))
        iso2 = np.minimum(iso2[:n], iso2[n:])
        single.append([int(np.sum(iso1 >= w)) for w in W])
        double.append([int(np.sum(iso2 >= w)) for w in W])
        n_tot.append(n)
        # no decoupling: static ZZ multiplets of unpolarized neighbours (half-width 3 sigma) + pulse resolution 4/tau
        d = P[:, None, :] - P[None, :, :]
        r = np.linalg.norm(d, axis=2)
        np.fill_diagonal(r, np.inf)
        cos = (d @ BDIR) / r
        J = D_NN * (1 - 3 * cos**2) / r**3
        h = 3 * 0.5 * np.sqrt(np.sum(J**2, axis=1))
        f = np.abs(nu + A / 2)
        tau = T1E_DEV / 20
        ok = [bool(np.all(np.abs(f[k] - np.delete(f, k)) >= h[k] + np.delete(h, k) + 4 / tau)) for k in range(n)]
        nodec.append(dict(n_addressable=int(np.sum(ok)), median_halfwidth_Hz=float(np.median(h))))
        print(f"cfg {len(n_tot)}: N = {n}, W=1 Hz single {single[-1][15]}, no-decoupling {nodec[-1]}", flush=True)
    single, double = np.array(single), np.array(double)
    w_dev = 4 / (S_MREV * min(T2_L / 10, T1E_DEV / 20))
    k_dev = int(np.argmin(np.abs(W - w_dev)))
    out = dict(K_Hz=K_LAW, W_Hz=W.tolist(), single_mean=single.mean(0).tolist(), double_mean=double.mean(0).tolist(),
               law_single=(0.627 * np.sqrt(K_LAW / W)).tolist(), law_double=(0.313 * np.sqrt(K_LAW / W)).tolist(),
               n_total_mean=float(np.mean(n_tot)), W_device_Hz=w_dev,
               device_single=float(np.interp(w_dev, W, single.mean(0))),
               device_fraction=float(np.interp(w_dev, W, single.mean(0)) / np.mean(n_tot)), nodecoupling=nodec,
               nearest_W_index=k_dev)
    for w in (0.3, 1, 3, 10, 30):
        m = float(np.interp(w, W, single.mean(0)))
        out[f"NR1_W={w}"] = dict(count=m, law=0.627 * np.sqrt(K_LAW / w), ratio=m / (0.627 * np.sqrt(K_LAW / w)))
    print(json.dumps({k: v for k, v in out.items() if k.startswith("NR1") or k.startswith("device")}, indent=1))
    (OUT / "law.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def hemisphere(m):
    i = np.arange(m) + 0.5
    z = i / m                                                       # uniform in z on [0, 1]: uniform on the hemisphere
    phi = np.pi * (1 + np.sqrt(5)) * i
    rho = np.sqrt(1 - z**2)
    return np.stack([rho * np.cos(phi), rho * np.sin(phi), z], 1)


def isolation(P, dirs, keep=False, chunk=500):
    """b_k = max over directions of the projected distance to the nearest other nucleus"""
    n = len(P)
    b, arg = np.zeros(n), np.zeros(n, int)
    full = np.empty((n, len(dirs)), np.float32) if keep else None
    for s in range(0, len(dirs), chunk):
        proj = P @ dirs[s:s + chunk].T
        o = np.argsort(proj, axis=0)
        ps = np.take_along_axis(proj, o, 0)
        g = np.diff(ps, axis=0)
        inf = np.full((1, g.shape[1]), np.inf)
        m = np.minimum(np.vstack([inf, g]), np.vstack([g, inf]))
        iso = np.empty_like(m)
        np.put_along_axis(iso, o, m, 0)
        if keep:
            full[:, s:s + chunk] = iso
        cm, ca = iso.max(1), iso.argmax(1) + s
        up = cm > b
        b[up], arg[up] = cm[up], ca[up]
    return b, arg, full


def refine(P, b, g0, levels=((0.03, 25), (0.004, 15))):
    """local search around each nucleus' best grid direction (two nested patches); returns a lower bound on b_k"""
    b, g = b.copy(), g0.copy()
    for k in range(len(P)):
        d = np.delete(P - P[k], k, axis=0)
        for rad, m in levels:
            u = np.cross(g[k], [1.0, 0, 0] if abs(g[k][0]) < 0.9 else [0, 1.0, 0])
            u /= np.linalg.norm(u)
            v = np.cross(g[k], u)
            a = np.linspace(-rad, rad, m)
            aa, bb = np.meshgrid(a, a)
            dirs = g[k] + aa.ravel()[:, None] * u + bb.ravel()[:, None] * v
            dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
            iso = np.min(np.abs(d @ dirs.T), axis=0)
            j = int(np.argmax(iso))
            if iso[j] > b[k]:
                b[k], g[k] = iso[j], dirs[j]
    return b, g


def cover(B, dstar):
    """fraction of nuclei with b_k >= d* for each configuration (rows) and each d* (columns)"""
    return np.array([1 - np.searchsorted(np.sort(b), dstar, side="left") / len(b) for b in B])


def coverage_map(B, caps, tcs, gGs):
    """for every T2 cap and t_c: the best coverage over gradient strengths (per configuration)"""
    res = {}
    for cap in caps:
        rows = []
        for tc in tcs:
            T2 = t2_model(tc, gGs * R_MAX, cap)
            dstar = 120.0 / (T2 * gGs)
            cov = cover(B, dstar)
            j = int(np.argmax(cov.mean(0)))
            rows.append(dict(tc_ms=float(tc), best_gG_Hz_per_nm=float(gGs[j]), T2_s=float(T2[j]), dstar_nm=float(dstar[j]),
                             mean_cov=float(cov[:, j].mean()), frac_cfg_ge_0p99=float(np.mean(cov[:, j] >= 0.99)),
                             frac_cfg_all=float(np.mean(cov[:, j] >= 1.0))))
        res["inf" if np.isinf(cap) else f"{cap:g}"] = rows
    return res


def run():
    """Part B: D+ mode, steerable gradient: best isolation distance of every 29Si within 6 nm; feasibility, time, energy"""
    register()
    dirs = hemisphere(N_DIR)
    g_fixed = BDIR
    B, Bgrid, R, fixed_iso, S, keep = [], [], [], [], [], None
    t0 = time.time()
    for c, (P, _) in enumerate(configs(2602)):
        b, arg, full = isolation(P, dirs, keep=(c == 0))
        b2, g2 = refine(P, b, dirs[arg])
        fi = lines_isolation(P @ g_fixed)
        dd = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
        np.fill_diagonal(dd, np.inf)
        S.append(np.sum(1 / dd, axis=1))
        B.append(b2), Bgrid.append(b), R.append(np.linalg.norm(P, axis=1)), fixed_iso.append(fi)
        if c == 0:
            keep = (P, full, b2, g2)
        print(f"cfg {c}: N = {len(P)}, b_k 1st pct grid {np.percentile(b, 1):.4f} refined {np.percentile(b2, 1):.4f} nm, "
              f"min {b2.min():.4f}, {time.time() - t0:.0f} s", flush=True)
    dgrid = np.logspace(-3.5, -0.5, 61)
    cov = cover(B, dgrid)
    cov_fixed = cover(fixed_iso, dgrid)
    law_fixed = np.array([[np.mean(np.exp(-d * s)) for d in dgrid] for s in S])
    d99 = np.array([np.percentile(b, 1) for b in B])
    i10 = int(np.argmin(np.abs(dgrid - 0.010)))
    tcs = np.logspace(np.log10(0.01), np.log10(2.0), 47)
    gGs = np.logspace(0, 4.5, 91)
    caps = [np.inf, 1000.0, 100.0, 25.0, 8.0]
    cmap = coverage_map(B, caps, tcs, gGs)
    feas = {}
    for key, rows in cmap.items():
        ok = [x for x in rows if x["mean_cov"] >= 0.99 and x["frac_cfg_ge_0p99"] >= 0.9]
        feas[key] = max(ok, key=lambda x: x["tc_ms"]) if ok else None
    ladd = coverage_map(B, [25.0], [TC_L], gGs)["25"][0]
    # time and energy at the largest feasible t_c of the uncapped law, greedy direction cover on configuration 0
    P0, full, b0, g0 = keep
    energy = None
    if feas["inf"]:
        fx = feas["inf"]
        need = b0 >= fx["dstar_nm"]
        M = full >= fx["dstar_nm"]
        unc = need.copy()
        unc[~M.any(1)] = False                                      # only refined directions isolate these
        extra = int(np.sum(need & ~M.any(1)))
        nd = 0
        while unc.any():
            cnt = M[unc].sum(0)
            j = int(np.argmax(cnt))
            if cnt[j] == 0:
                break
            unc &= ~M[:, j]
            nd += 1
        n_dir = nd + extra
        tau = fx["T2_s"] / 10
        dr = drive(fx["tc_ms"])
        energy = dict(t_c_ms=fx["tc_ms"], gG_Hz_per_nm=fx["best_gG_Hz_per_nm"], G_mT_per_nm=fx["best_gG_Hz_per_nm"] / GN * 1e3,
                      T2_s=fx["T2_s"], tau_s=tau, n_directions=n_dir, grid_directions=nd, refined_only=extra,
                      total_time_s=n_dir * tau, drive=dr, energy_J=dr["p_avg_W"] * n_dir * tau,
                      energy_per_atom_J=dr["p_avg_W"] * n_dir * tau / len(P0), kT_ln2_100mK_J=1.380649e-23 * 0.1 * np.log(2))
    out = dict(n_mean=float(np.mean([len(b) for b in B])), dgrid_nm=dgrid.tolist(), coverage_mean=cov.mean(0).tolist(),
               coverage_min=cov.min(0).tolist(), coverage_fixed_mean=cov_fixed.mean(0).tolist(),
               law_fixed_mean=law_fixed.mean(0).tolist(), d99_nm=d99.tolist(), d99_median_nm=float(np.median(d99)),
               b_grid_p1=[float(np.percentile(b, 1)) for b in Bgrid],
               NR3=dict(frac_cfg_99_at_0p010=float(np.mean(d99 >= 0.010)), fixed_mean_at_0p010=float(cov_fixed[:, i10].mean())),
               feasibility=cmap, feasible_largest_tc=feas, ladd_point=ladd, energy=energy,
               drive_ladd=drive(TC_L), S_mean=float(np.mean([s.mean() for s in S])),
               b_vs_r_cfg0=dict(r=R[0].round(3).tolist(), b=B[0].round(5).tolist()))
    print(json.dumps(dict(NR3=out["NR3"], d99_median=out["d99_median_nm"], feasible=feas, ladd=out["ladd_point"],
                          energy=energy), indent=1, default=float))
    (OUT / "run.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")


def bus():
    """Part C: dipolar bus connectivity, hop depth from the hub, electron couplings, Lie rank of 4-nucleus clusters"""
    register()
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components, shortest_path
    rng = np.random.default_rng(2603)
    rows, lie = [], []
    for c, (P, Aiso) in enumerate(configs(2602)):
        d = P[:, None, :] - P[None, :, :]
        r = np.linalg.norm(d, axis=2)
        np.fill_diagonal(r, np.inf)
        J = D_NN * (1 - 3 * ((d @ BDIR) / r) ** 2) / r**3
        A = couplings_dir(P, Aiso, BDIR)
        hub = int(np.argmax(np.abs(A)))
        row = dict(N=len(P), frac_A_ge_1Hz=float(np.mean(np.abs(A) >= 1)), frac_A_ge_10Hz=float(np.mean(np.abs(A) >= 10)))
        for T2 in (25.0, 100.0, 1000.0):
            G = csr_matrix(np.abs(J) >= 5 / T2)
            nc, lab = connected_components(G, directed=False)
            hops = shortest_path(G, unweighted=True, indices=hub)
            row[f"T2={T2:g}"] = dict(components=int(nc), largest=int(np.bincount(lab).max()),
                                     max_hops_from_hub=float(np.max(hops[np.isfinite(hops)])),
                                     mean_degree=float(G.sum() / len(P)))
        rows.append(row)
        if c < 10:
            k = int(rng.integers(len(P)))
            nb = np.argsort(-np.abs(J[k]))[:3]
            idx = [k] + nb.tolist()
            g = rng.normal(size=3)
            g /= np.linalg.norm(g)
            f = 100.0 * (P[idx] @ g)
            lie.append(lie_cluster(f, J[np.ix_(idx, idx)]))
        print(c, row["T2=25"], lie[-1] if c < 10 else "", flush=True)
    out = dict(configs=rows, lie=lie,
               NR6=dict(connected_T2_25=int(sum(x["T2=25"]["components"] == 1 for x in rows)),
                        lie_selective_full=int(sum(x["selective"] == 255 for x in lie)),
                        lie_global_full=int(sum(x["global"] == 255 for x in lie))))
    print(json.dumps(out["NR6"]))
    (OUT / "bus.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def lie_cluster(f, J):
    n = len(f)
    pau = {"x": np.array([[0, 1], [1, 0]], complex) / 2, "y": np.array([[0, -1j], [1j, 0]], complex) / 2,
           "z": np.array([[1, 0], [0, -1]], complex) / 2}

    def op(s, k):
        out = np.array([[1.0 + 0j]])
        for j in range(n):
            out = np.kron(out, pau[s] if j == k else np.eye(2))
        return out
    H0 = sum(f[k] * op("z", k) for k in range(n)) + sum(J[k, l] * op("z", k) @ op("z", l) for k in range(n)
                                                         for l in range(k + 1, n))
    H0 = H0 / np.linalg.norm(H0)
    sel = [H0] + [op(s, k) for k in range(n) for s in ("x", "y")]
    glo = [H0, sum(op("x", k) for k in range(n)), sum(op("y", k) for k in range(n))]
    return {"offsets_Hz": np.round(f, 2).tolist(), "selective": lie_rank(sel)[0], "global": lie_rank(glo)[0]}


def fixed():
    """correction after run: the fixed reference [1,2,3] is a lattice direction (sites stack in planes 0.036 nm apart), so it
    isolates almost nothing for a reason unrelated to crowding; the reference is redone with 10 generic random directions per
    configuration (same configurations) and compared with law L2; the lattice-direction values are kept under *_lattice_123"""
    Rn = json.loads((OUT / "run.json").read_text(encoding="utf-8"))
    dgrid = np.array(Rn["dgrid_nm"])
    rng = np.random.default_rng(2604)
    cov, law_g = [], []
    for P, _ in configs(2602):
        dd = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
        np.fill_diagonal(dd, np.inf)
        S = np.sum(1 / dd, axis=1)
        for _ in range(10):
            g = rng.normal(size=3)
            g /= np.linalg.norm(g)
            cov.append(cover([lines_isolation(P @ g)], dgrid)[0])
        law_g.append([float(np.mean(np.exp(-d * S))) for d in dgrid])
    i10 = int(np.argmin(np.abs(dgrid - 0.010)))
    if "coverage_fixed_lattice_123_mean" not in Rn:
        Rn["coverage_fixed_lattice_123_mean"] = Rn["coverage_fixed_mean"]
        Rn["NR3"]["fixed_lattice_123_at_0p010"] = Rn["NR3"].pop("fixed_mean_at_0p010")
    Rn["coverage_fixed_mean"] = np.mean(cov, 0).tolist()
    Rn["law_fixed_mean"] = np.mean(law_g, 0).tolist()
    Rn["NR3"]["fixed_generic_mean_at_0p010"] = float(np.mean(cov, 0)[i10])
    Rn["NR3"]["law_L2_at_0p010"] = float(np.mean(law_g, 0)[i10])
    for d in (0.001, 0.003, 0.010):
        i = int(np.argmin(np.abs(dgrid - d)))
        print(f"d* = {dgrid[i]:.4f} nm: generic fixed direction {np.mean(cov, 0)[i]:.4f}, L2 {np.mean(law_g, 0)[i]:.4f}, "
              f"[1,2,3] {Rn['coverage_fixed_lattice_123_mean'][i]:.4f}")
    (OUT / "run.json").write_text(json.dumps(Rn, indent=1, default=float), encoding="utf-8")


def bloch_flip(delta_tau, t2_tau, n=4000):
    """Gaussian pi pulse of unit length (sigma = 1/6, +-3 sigma), detuning in units of 1/tau, transverse decay T2 in units
    of tau (Bloch equations, RK4): population moved out of the initial state"""
    t = np.linspace(-0.5, 0.5, n + 1)
    env = np.exp(-t**2 / (2 * (1 / 6) ** 2))
    env *= np.pi / np.trapezoid(env, t)                               # area pi (rad)
    w = 2 * np.pi * delta_tau

    def f(m, om):
        x, y, z = m
        return np.array([-w * y - x / t2_tau, w * x - om * z - y / t2_tau, om * y])
    m, h = np.array([0.0, 0.0, 1.0]), t[1] - t[0]
    for i in range(n):
        o0, o1 = env[i], env[i + 1]
        om = 0.5 * (o0 + o1)
        k1 = f(m, o0)
        k2 = f(m + h / 2 * k1, om)
        k3 = f(m + h / 2 * k2, om)
        k4 = f(m + h * k3, o1)
        m = m + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return float((1 - m[2]) / 2)


def pulse():
    """check of the selection rule W = 4/(s tau): inversion of the target with T2 = 10 tau and excitation of a neighbour
    4/tau, 2/tau or 1/tau away"""
    register()
    out = dict(target_T2_10tau=bloch_flip(0.0, 10.0), target_no_decay=bloch_flip(0.0, 1e12),
               neighbour_4_over_tau=bloch_flip(4.0, 1e12), neighbour_2_over_tau=bloch_flip(2.0, 1e12),
               neighbour_1_over_tau=bloch_flip(1.0, 1e12))
    print(json.dumps(out, indent=1))
    (OUT / "pulse.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def measured():
    """Part E, measured parameters only (Ladd 2005: T2 = 25 s, offsets <= 120 Hz, s = 1/3; no t_c extrapolation).
    A register of radius R gets gG = 120 Hz / R so that every register offset stays in the measured range; nuclei outside
    the register are environment (collateral rotations there do not act on the register). Trade-off: pulse length
    tau = x T2 (Bloch inversion with decay) and neighbour distance kappa/tau (crosstalk) set d* = kappa R / (s x T2 120 Hz).
    Output: coverage of every register and the largest radius in which all 29Si are writable, for each (x, kappa)"""
    pf = OUT / "predictions_nr26_measured.json"
    if not pf.exists():
        raise SystemExit("register predictions_nr26_measured.json first")
    xs = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5]
    kappas = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    F = {x: bloch_flip(0.0, 1 / x) for x in xs}
    eps = {k: bloch_flip(k, 1e12) for k in kappas}
    print("inversion vs tau/T2:", {x: round(v, 4) for x, v in F.items()})
    print("crosstalk vs kappa:", {k: float(f"{v:.2e}") for k, v in eps.items()}, flush=True)
    radii = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
    dirs = hemisphere(N_DIR)
    B = {}
    t0 = time.time()
    for R in radii:
        B[R] = []
        for c, (P, _) in enumerate(configs(2602)):
            if c >= 10:
                break
            Q = P[np.linalg.norm(P, axis=1) <= R]
            b, arg, _ = isolation(Q, dirs)
            B[R].append(refine(Q, b, dirs[arg])[0])
        print(f"R = {R} nm: N = {np.mean([len(b) for b in B[R]]):.0f}, b_99 median "
              f"{np.median([np.percentile(b, 1) for b in B[R]]):.4f} nm, {time.time() - t0:.0f} s", flush=True)
    rows = []
    for x in xs:
        for k in kappas:
            per = []
            for R in radii:
                d = k * R / (S_MREV * x * T2_L * W_L)
                cv = cover(B[R], np.array([d]))[:, 0]
                per.append(dict(R_nm=R, dstar_nm=d, mean_cov=float(cv.mean()), frac_cfg_ge_0p99=float(np.mean(cv >= 0.99)),
                                frac_cfg_all=float(np.mean(cv >= 1.0)), n_mean=float(np.mean([len(b) for b in B[R]]))))
            ok = [p["R_nm"] for p in per if p["mean_cov"] >= 0.99 and p["frac_cfg_ge_0p99"] >= 0.9]
            r_all = max(ok) if ok else 0.0
            rows.append(dict(x=x, kappa=k, tau_s=x * T2_L, inversion=F[x], crosstalk=eps[k], R_all_nm=r_all,
                             n_all=float(next((p["n_mean"] for p in per if p["R_nm"] == r_all), 0.0)),
                             cov_6nm=per[-1]["mean_cov"], per_R=per))
    out = dict(T2_s=T2_L, offset_max_Hz=W_L, s=S_MREV, inversion=F, crosstalk=eps, rows=rows,
               b99_median_nm={R: float(np.median([np.percentile(b, 1) for b in B[R]])) for R in radii})
    for r in rows:
        if r["kappa"] in (2.0, 4.0):
            print(f"x = {r['x']:<4} kappa = {r['kappa']}: inversion {r['inversion']:.3f}, crosstalk {r['crosstalk']:.1e}, "
                  f"R_all {r['R_all_nm']} nm (N ~ {r['n_all']:.0f}), 6 nm coverage {r['cov_6nm']:.3f}", flush=True)
    (OUT / "measured.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")


if __name__ == "__main__":
    {"law": law, "run": run, "bus": bus, "pulse": pulse, "fixed": fixed, "measured": measured}[
        sys.argv[1] if len(sys.argv) > 1 else "run"]()
