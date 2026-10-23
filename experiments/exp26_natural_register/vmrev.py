"""Exp26 part F: exact many-spin check of the core assumption behind the measured-point result (summary_exp26.md §6):
a selective write embedded in MREV-16 decoupling acts like a single-spin selective pulse on the scaled offsets, with the
real 29Si-29Si dipolar couplings present.
 - cluster: target + 3 nuclei nearest in frequency along its best gradient direction + strongest dipolar partners (N = 7),
   natural configuration 0 of nr26; offsets f_j = gG (R_j - R_t).g with gG = 20 Hz/nm (measured offset range);
   secular dipolar H_D = sum J_kl (Iz Iz - (Ix Ix + Iy Iy)/2), J = 4.748 Hz nm^3 (1 - 3 cos^2) / r^3, B0 || [1,2,3]
 - MREV-8: tau X tau -Y 2tau Y tau -X 2tau -X tau Y 2tau -Y tau X tau (12 tau); MREV-16 = MREV-8 + MREV-8 with the x phases
   inverted (effective field along z, scaling 1/3 as in Ladd 2005); finite pi/2 pulses of 15 us at t_c = 1.568 ms
 - write: one small global pulse after each cycle (DANTE), Gaussian area pi over 17.5 s (the measured-point design)
 - storage: the same without the write (leakage of computational basis states)
usage: python vmrev.py  -> results/exp26/vmrev.json (predictions_nr26_vmrev.json written before)"""
import json
import time

import numpy as np
from scipy.linalg import expm

from nr26 import BDIR, D_NN, N_DIR, OUT, S_MREV, TC_L, bloch_flip, configs, hemisphere, isolation, refine

P_W, TC, GG, TAU_W = 15e-6, TC_L * 1e-3, 20.0, 17.5
PHASE = {"x": 0.0, "y": np.pi / 2, "-x": np.pi, "-y": 3 * np.pi / 2}
M8 = ["x", "-y", "y", "-x", "-x", "y", "-y", "x"]
W8 = [1, 1, 2, 1, 2, 1, 2, 1, 1]                                    # windows (tau) around the 8 pulses


def ops(n):
    s = {"x": np.array([[0, 1], [1, 0]], complex) / 2, "y": np.array([[0, -1j], [1j, 0]], complex) / 2,
         "z": np.array([[1, 0], [0, -1]], complex) / 2}

    def op(a, k):
        out = np.array([[1.0 + 0j]])
        for j in range(n):
            out = np.kron(out, s[a] if j == k else np.eye(2))
        return out
    return {a: [op(a, k) for k in range(n)] for a in "xyz"}


def couplings(P):
    n = len(P)
    J = np.zeros((n, n))
    for k in range(n):
        for l in range(k + 1, n):
            d = P[l] - P[k]
            r = np.linalg.norm(d)
            J[k, l] = J[l, k] = D_NN * (1 - 3 * (d @ BDIR / r) ** 2) / r**3
    return J


def cycle(f, J, O):
    """MREV-16 cycle propagator with finite pulses (time ordered)"""
    n = len(f)
    H = sum(2 * np.pi * f[k] * O["z"][k] for k in range(n))
    for k in range(n):
        for l in range(k + 1, n):
            if J[k, l]:
                H = H + 2 * np.pi * J[k, l] * (O["z"][k] @ O["z"][l] - 0.5 * (O["x"][k] @ O["x"][l] + O["y"][k] @ O["y"][l]))
    Sx, Sy = sum(O["x"]), sum(O["y"])
    om = (np.pi / 2) / P_W
    tau = TC / 24
    pulses = M8 + [{"x": "-x", "-x": "x"}.get(p, p) for p in M8]
    windows = W8[:-1] + [W8[-1] + W8[0]] + W8[1:]                  # two MREV-8 blocks share the junction window
    U = np.eye(H.shape[0], dtype=complex)
    for i, p in enumerate(pulses):
        w = windows[i] * tau - (P_W / 2 if i == 0 else P_W)
        U = expm(-1j * H * w) @ U
        ph = PHASE[p]
        U = expm(-1j * (H + om * (np.cos(ph) * Sx + np.sin(ph) * Sy)) * P_W) @ U
    U = expm(-1j * H * (windows[-1] * tau - P_W / 2)) @ U
    return U


def effective(f):
    """single-spin effective rotation per cycle: axis and scaling (angle / (2 pi f t_c))"""
    O = ops(1)
    U = cycle(np.array([f]), np.zeros((1, 1)), O)
    U = U / np.sqrt(np.linalg.det(U))
    ang = 2 * np.arccos(np.clip(np.real(U[0, 0] + U[1, 1]) / 2, -1, 1))
    v = np.array([np.imag(U[0, 1] + U[1, 0]), np.real(U[1, 0] - U[0, 1]), np.imag(U[0, 0] - U[1, 1])]) / 2
    axis = v / (np.linalg.norm(v) + 1e-30)
    return float(ang / (2 * np.pi * f * TC)) if f else 0.0, axis.round(4).tolist()


def evolve(U, n, write, states):
    """evolve computational basis states through M cycles; write = Gaussian DANTE pulses about x"""
    M = int(round(TAU_W / TC))
    t = (np.arange(M) + 0.5) * TC
    th = np.exp(-(t - TAU_W / 2) ** 2 / (2 * (TAU_W / 6) ** 2))
    th *= np.pi / th.sum()
    dim = 2**n
    Psi = np.zeros((dim, len(states)), complex)
    for c, s in enumerate(states):
        Psi[int("".join(map(str, s)), 2), c] = 1.0
    for m in range(M):
        Psi = U @ Psi
        if write:
            r = np.array([[np.cos(th[m] / 2), -1j * np.sin(th[m] / 2)], [-1j * np.sin(th[m] / 2), np.cos(th[m] / 2)]])
            Psi = Psi.reshape([2] * n + [len(states)])
            for k in range(n):
                Psi = np.moveaxis(np.tensordot(r, Psi, axes=([1], [k])), 0, k)
            Psi = Psi.reshape(dim, len(states))
    pop = np.abs(Psi) ** 2
    bits = np.array([[(i >> (n - 1 - k)) & 1 for k in range(n)] for i in range(dim)])     # 1 = spin down
    flips = []
    for c, s in enumerate(states):
        p1 = pop[:, c] @ bits                                        # probability each spin is down
        flips.append([float(p1[k] if s[k] == 0 else 1 - p1[k]) for k in range(n)])
    return np.array(flips)


def main():
    pf = OUT / "predictions_nr26_vmrev.json"
    if not pf.exists():
        raise SystemExit("register predictions_nr26_vmrev.json first")
    t0 = time.time()
    scale, axis = effective(50.0)
    print(f"MREV-16 effective field per cycle: scaling {scale:.4f} (Ladd: 1/3), axis {axis}", flush=True)
    P, _ = next(configs(2602))
    dirs = hemisphere(N_DIR)
    b, arg, _ = isolation(P, dirs)
    b2, g2 = refine(P, b, dirs[arg])
    dstar = 2.0 * 6.0 / (S_MREV * TAU_W * 120.0)
    dd = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    np.fill_diagonal(dd, np.inf)
    Jall = D_NN * (1 - 3 * (((P[:, None, :] - P[None, :, :]) @ BDIR) / dd) ** 2) / dd**3
    nn = dd.min(1) < 0.24
    ok = b2 >= dstar
    targets = {}
    c1 = np.where(ok & ~nn)[0]
    targets["no_NN"] = int(c1[np.argmin(b2[c1])])
    strong = np.where(ok & nn & (np.max(np.abs(np.where(dd < 0.24, Jall, 0)), axis=1) >= 250))[0]
    if len(strong):
        targets["NN_J_ge_250Hz"] = int(strong[np.argmin(b2[strong])])
        jmax = np.max(np.abs(np.where(dd < 0.24, Jall, 0)), axis=1)
        targets["NN_strongest"] = int(strong[np.argmax(jmax[strong])])
    out = dict(effective_scaling=scale, effective_axis=axis, dstar_nm=dstar, fraction_with_NN=float(nn.mean()), targets={})
    rng = np.random.default_rng(2606)
    for name, t in targets.items():
        proj = (P - P[t]) @ g2[t]
        near = [int(j) for j in np.argsort(np.abs(proj)) if j != t][:3]
        strongJ = [int(j) for j in np.argsort(-np.abs(Jall[t])) if j != t]
        idx = [t] + near
        for j in strongJ:
            if len(idx) >= 7:
                break
            if j not in idx:
                idx.append(j)
        Pc, f = P[idx], GG * proj[idx]
        J = couplings(Pc)
        n = len(idx)
        O = ops(n)
        states = [rng.integers(0, 2, n).tolist() for _ in range(16)] + [[0] * n]
        res = dict(nuclei=idx, offsets_Hz=f.round(4).tolist(), J_to_target_Hz=J[0].round(2).tolist(), b_target_nm=float(b2[t]),
                   max_abs_J_Hz=float(np.abs(J).max()))
        for label, useJ in (("dipolar", True), ("no_dipolar", False)):
            U = cycle(f, J if useJ else np.zeros_like(J), O)
            w = evolve(U, n, True, states)
            s = evolve(U, n, False, states)
            res[label] = dict(target_flip=float(w[:, 0].mean()), max_other_flip=float(w[:, 1:].max()),
                              mean_other_flip=float(w[:, 1:].mean()), storage_max_leak=float(s.max()),
                              storage_mean_leak=float(s.mean()),
                              storage_target_leak=float(s[:, 0].max()))
            print(f"{name} [{label}]: target flip {res[label]['target_flip']:.4f}, max other {res[label]['max_other_flip']:.2e}, "
                  f"storage leak max {res[label]['storage_max_leak']:.2e} ({time.time() - t0:.0f} s)", flush=True)
        res["bloch_crosstalk"] = [float(bloch_flip(S_MREV * abs(x) * TAU_W, 1e12)) for x in f[1:]]
        out["targets"][name] = res
    (OUT / "vmrev.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk in ("dipolar", "no_dipolar", "max_abs_J_Hz")}
                      for k, v in out["targets"].items()}, indent=1))


def survey():
    """storage leakage of 60 targets (20 NN-pair members, 20 coupled >= 5 Hz to a pair member, 20 others) at three cycle
    times (pulse width scaled with t_c); clusters = target + 6 strongest dipolar partners; class fractions of the register"""
    global TC, P_W
    pf = OUT / "predictions_nr26_survey.json"
    if not pf.exists():
        raise SystemExit("register predictions_nr26_survey.json first")
    P, _ = next(configs(2602))
    dirs = hemisphere(N_DIR)
    b, arg, _ = isolation(P, dirs)
    _, g2 = refine(P, b, dirs[arg])
    d = P[:, None, :] - P[None, :, :]
    dd = np.linalg.norm(d, axis=2)
    np.fill_diagonal(dd, np.inf)
    Jall = D_NN * (1 - 3 * ((d @ BDIR) / dd) ** 2) / dd**3
    member = dd.min(1) < 0.24
    coupled = ~member & (np.max(np.where(member[None, :], np.abs(Jall), 0), axis=1) >= 5)
    other = ~member & ~coupled
    rng = np.random.default_rng(2607)
    picks = {name: rng.choice(np.where(mask)[0], 20, replace=False) for name, mask in
             (("member", member), ("coupled", coupled), ("other", other))}
    tc0, pw0 = TC, P_W
    res = {name: {} for name in picks}
    for tcf in (1, 4, 16):
        TC, P_W = tc0 / tcf, pw0 / tcf
        M = int(round(TAU_W / TC))
        for name, ts in picks.items():
            leaks = []
            for t in ts:
                idx = [int(t)] + [int(j) for j in np.argsort(-np.abs(Jall[t]))[:6]]
                f = GG * ((P[idx] - P[t]) @ g2[t])
                n = len(idx)
                U = np.linalg.matrix_power(cycle(f, couplings(P[idx]), ops(n)), M)
                states = [rng.integers(0, 2, n).tolist() for _ in range(15)]
                bits0 = np.array([s[0] for s in states])
                cols = [int("".join(map(str, s)), 2) for s in states]
                pop = np.abs(U[:, cols]) ** 2
                down0 = pop[(np.arange(2**n) >> (n - 1)) & 1 == 1].sum(0)
                leaks.append(float(np.max(np.where(bits0 == 0, down0, 1 - down0))))
            res[name][f"t_c={TC * 1e3:.3f} ms"] = dict(leaks=leaks, frac_gt_0p1=float(np.mean(np.array(leaks) > 0.1)),
                                                      frac_gt_0p01=float(np.mean(np.array(leaks) > 0.01)),
                                                      median=float(np.median(leaks)))
            print(f"t_c = {TC * 1e3:.3f} ms {name}: median leak {np.median(leaks):.3g}, > 0.1: {np.mean(np.array(leaks) > 0.1):.2f}, "
                  f"> 0.01: {np.mean(np.array(leaks) > 0.01):.2f}", flush=True)
    TC, P_W = tc0, pw0
    frac = dict(member=float(member.mean()), coupled=float(coupled.mean()), other=float(other.mean()))
    keep = {}
    for key in res["member"]:
        keep[key] = float(sum(frac[c] * (1 - res[c][key]["frac_gt_0p1"]) for c in frac))
    out = dict(class_fractions=frac, survey=res, register_fraction_leak_le_0p1=keep)
    print(json.dumps(dict(class_fractions=frac, keep=keep), indent=1))
    (OUT / "vmrev_survey.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    import sys
    survey() if len(sys.argv) > 1 and sys.argv[1] == "survey" else main()
