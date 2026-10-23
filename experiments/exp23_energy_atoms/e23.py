"""Exp23: (1) where does the energy of a recycled cell go, and how much is lost; (2) could these cells move or program
atoms? Device: the exp18 3.3 V state (seed 1, three cells), implicit solver.
Part 1: the rightmost cell is killed (-300 mV band, 30 um margins) slowly (30 s ramp) and abruptly; the power is split
into supply = load resistor + ballasts + drain rail + transistor channels (cells) from the full profiles, before and after.
Part 2: quantities from the simulated states (largest drain voltage, current densities, power density) against the
thresholds for moving or changing atoms.
usage: python e23.py [V_DD = 3.3 | 5]  ->  results/exp23/energy_atoms_v<V_DD>.json"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp18_programmable"))
import regen18 as rg  # noqa: E402
from run18 import law, sp, prog_array, W, DY, LAM, R_D, R_L  # noqa: E402

OUT = ROOT / "results" / "exp23"
VH = float(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "fig" else 3.3
KW = dict(rg.KW)
R_V = R_D * LAM**2
N = int(round(W / DY))
Y = (np.arange(N) + 0.5) * DY
D = np.random.default_rng(1).normal(0.0, 3.356e-3 / np.sqrt(DY * 0.15), N)
Q = 1.602176634e-19


def step(V, regions, dt, store=False):
    r = sp.simulate(vdd_up=(VH, VH), dvdd=1.0, prog=prog_array(regions), V0=V, t_relax=dt, n_sub=2 if dt < 1 else 12,
                    implicit_load=True, store_at=[("up", VH)] if store else (),
                    **{k: v for k, v in KW.items() if k not in ("t_relax",)})
    return np.array(r["V_final"]), r


def partition(r, regions):
    """power split of a stored state (W); per-cell channel power"""
    Vb = np.array(r["profiles"][f"up_{VH:.2f}"])
    VD = np.array(r["profiles"][f"up_{VH:.2f}_VD"])
    h = r["hist"][0]
    p = law.BASE
    jd = sp.j_d(Vb, p["vov0"] + D + prog_array(regions), p["ss"])     # A per um
    I = float(jd.sum() * DY)
    vc = h["vds"]
    P = dict(supply=VH * I, load=(VH - vc) * I, ballast=float(np.sum((vc - VD) ** 2 * DY / R_V)),
             rail=float(np.sum(np.diff(VD) ** 2 / (R_D * DY))), channel=float(np.sum(VD * jd * DY)))
    P["balance_error"] = P["supply"] - P["load"] - P["ballast"] - P["rail"] - P["channel"]
    cells = rg.intervals(h)
    per_cell = []
    for s, e in cells:
        m = (Y >= s) & (Y < e)
        per_cell.append(dict(interval=[s, e], channel_W=float(np.sum(VD[m] * jd[m] * DY)),
                             body_energy_J=float(np.sum(0.5 * 1e-15 * Vb[m] ** 2 * DY))))
    return dict(P=P, cells=per_cell, I_uA=I * 1e6, vc=vc, VD_max=float(VD.max()), j_max_uA_um=float(jd.max() * 1e6))


def scenario(kind):
    g = json.loads((ROOT / "results" / "exp18" / "regen" / f"grow_v{VH:g}.json").read_text(encoding="utf-8"))
    V = np.array(g["V_final"])
    cells = rg.intervals(g["hist"][-1])
    t = cells[-1]
    zone = ((t[0] + t[1]) / 2, t[1] - t[0] + 60.0)
    V, r = step(V, [], 1.0, store=True)
    before = partition(r, [])
    supply_t = []
    if kind == "slow":
        for k in range(30):
            regs = [(zone[0], zone[1], -0.30 * (k + 1) / 30)]
            V, r = step(V, regs, 1.0)
            supply_t.append((k + 1.0, r["hist"][0]["I_uA"], 1.0))
    else:
        regs = [(zone[0], zone[1], -0.30)]
        for k in range(50):
            V, r = step(V, regs, 0.02)
            supply_t.append((0.02 * (k + 1), r["hist"][0]["I_uA"], 0.02))
    regs = [(zone[0], zone[1], -0.30)]
    for k in range(5):
        V, r = step(V, regs, 1.0)
        supply_t.append((supply_t[-1][0] + 1.0, r["hist"][0]["I_uA"], 1.0))
    V, r = step(V, regs, 1.0, store=True)
    after = partition(r, regs)
    # one-time cost of the transition: extra supply energy relative to the final steady state (step-end samples)
    e_extra = float(sum(VH * (i_uA - after["I_uA"]) * 1e-6 * dt for _, i_uA, dt in supply_t))
    dead = before["cells"][-1]["channel_W"]
    surv_before = sum(c["channel_W"] for c in before["cells"][:-1])
    rec = after["P"]["channel"] - surv_before
    return dict(before=before, after=after, current_trace_uA=supply_t, transition_extra_supply_energy_J=e_extra,
                transition_time_s=float(supply_t[-1][0]), dead_cell_channel_W=dead, recycled_into_cells_W=rec, recycling_efficiency=rec / dead,
                total_channel_change_W=after["P"]["channel"] - before["P"]["channel"],
                loss_change_W={k: after["P"][k] - before["P"][k] for k in ("load", "ballast", "rail")},
                body_energy_of_dead_cell_J=before["cells"][-1]["body_energy_J"],
                n_cells_before=len(before["cells"]), n_cells_after=len(after["cells"]))


def atoms(state):
    """numbers from a simulated state against atomic thresholds"""
    vdmax = state["VD_max"]
    j_lin = state["j_max_uA_um"] * 1e-6                     # A per um of width
    lam_cur = j_lin * LAM                                   # current carried along the rail ~ j * lambda (A)
    t_inv = 3e-9                                            # inversion-layer thickness (m), exp5 centroid ~2.5 nm
    j_channel = j_lin / (1e-6 * t_inv)                      # A/m^2 in the channel
    via = 0.17e-6 ** 2
    j_contact = j_lin * 10.0 / via                          # 10 um segment through one licon-sized contact
    li1 = 1.28e-6 * 0.1e-6
    j_li1 = lam_cur / li1
    p_lin = vdmax * j_lin * 1e6                             # W per m of width
    dT = {f"k={k} W/mK": p_lin / (np.pi * k) * np.log(1e3) for k in (1.0, 10.0, 200.0)}   # line heater on a half-space
    l_pinch = (10e-9, 30e-9)
    E_field = [(vdmax - 0.1) / l for l in l_pinch]
    mu0 = 4e-7 * np.pi
    i_cell = j_lin * 40.0                                   # a ~40 um cell
    extra = dict(
        oxide_field_V_per_nm=(0.5 / 4.148),                 # |V_G - V_body| ~ 0.5 V over t_ox (vertical, gate oxide)
        oersted_field_T_at_1um=mu0 * i_cell / (2 * np.pi * 1e-6),
        electron_zeeman_eV=9.274e-24 * mu0 * i_cell / (2 * np.pi * 1e-6) / Q, kT_4K_eV=1.380649e-23 * 4 / Q,
        acceptors_switching_charge_per_um=1e-15 * 0.7 / Q,  # c_b * body swing / q, collective and not addressable
        dopants_in_min_channel=1.7e17 * (0.42e-4 * 0.15e-4) * 70e-7,
        volume_per_dopant_nm3=1e21 / 1.7e17,
        hci_slowdown_vs_nominal_1p8V=float(np.exp(40.0 * (1 / vdmax - 1 / 1.8))),   # lifetime ~ exp(b/V_D), b ~ 40 V
        max_tail_energy_eV=2 * vdmax)                      # electron-electron scattering tail bound
    return dict(**extra,
        max_carrier_energy_eV=vdmax, thresholds_eV=dict(Si_impact_ionisation=1.1, SiO2_barrier=3.1,
                                                        Si_H_bond_single=3.6, B_diffusion_activation=3.5,
                                                        Si_lattice_displacement=13.0),
        max_lateral_field_V_per_nm=[e * 1e-9 for e in E_field],
        field_scale_for_moving_atoms_V_per_nm=dict(field_evaporation=10.0, STM_manipulation=1.0),
        current_density_A_cm2=dict(channel=j_channel * 1e-4, contact=j_contact * 1e-4, li1_rail=j_li1 * 1e-4),
        electromigration_onset_A_cm2=dict(Al_or_Cu_at_300K=1e5, note="thermally activated: negligible at 4 K"),
        joule_heating_dT_K=dT, power_per_um_W=p_lin * 1e-6)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = {k: scenario(k) for k in ("slow", "abrupt")}
    res["atoms"] = atoms(res["slow"]["before"])
    (OUT / f"energy_atoms_v{VH:g}.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    for k in ("slow", "abrupt"):
        s = res[k]
        print(f"{k}: transition {s['transition_time_s']:.1f} s, extra supply energy {s['transition_extra_supply_energy_J']*1e6:.2f} uJ")
        print(f"{k}: cells {s['n_cells_before']} -> {s['n_cells_after']} | dead cell {s['dead_cell_channel_W']*1e6:.1f} uW | "
              f"recycled into cells {s['recycled_into_cells_W']*1e6:.1f} uW (eta {s['recycling_efficiency']:.3f}) | "
              f"loss change (uW) { {kk: round(v*1e6, 2) for kk, v in s['loss_change_W'].items()} } | balance "
              f"{s['before']['P']['balance_error']*1e6:.3f}/{s['after']['P']['balance_error']*1e6:.3f} uW")
        print("   before (uW):", {kk: round(v * 1e6, 2) for kk, v in s["before"]["P"].items()})
        print("   after  (uW):", {kk: round(v * 1e6, 2) for kk, v in s["after"]["P"].items()})
    print("atoms:", json.dumps(res["atoms"], indent=1))


def fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, axs = plt.subplots(1, 2, figsize=(17, 5.5), gridspec_kw=dict(width_ratios=[1, 1.2]))
    ax = axs[0]
    labels, parts = [], {k: [] for k in ("channel", "load", "ballast", "rail")}
    for v in ("3.3", "5"):
        d = json.loads((OUT / f"energy_atoms_v{v}.json").read_text(encoding="utf-8"))
        for tag, P, n in ((f"{v} V önce", d["slow"]["before"]["P"], d["slow"]["n_cells_before"]),
                          (f"{v} V yavaş ölüm", d["slow"]["after"]["P"], d["slow"]["n_cells_after"]),
                          (f"{v} V ani ölüm", d["abrupt"]["after"]["P"], d["abrupt"]["n_cells_after"])):
            labels.append(f"{tag}\n({n} hücre)")
            for k in parts:
                parts[k].append(P[k] * 1e6)
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    for k, lab, col in (("channel", "hücreler (kanal)", "C2"), ("ballast", "balastlar", "C1"), ("rail", "dren rayı", "C4"),
                        ("load", "yük direnci R_L", "0.6")):
        ax.bar(x, parts[k], bottom=bottom, color=col, label=lab)
        bottom += np.array(parts[k])
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("güç (µW)")
    ax.set_title("(a) hücre öldürülünce gücü kalanlara geçiyor (%98,5–102,3);\nasıl kayıp her durumda yük direnci", fontsize=10)
    ax.legend(fontsize=8)
    ax = axs[1]
    a = json.loads((OUT / "energy_atoms_v3.3.json").read_text(encoding="utf-8"))["atoms"]
    rows = [("taşıyıcı enerjisi / Si–H bağı (3,6 eV)", a["max_carrier_energy_eV"] / 3.6),
            ("kuyruk enerjisi (2qV_D) / Si–H", a["max_tail_energy_eV"] / 3.6),
            ("taşıyıcı enerjisi / kafes atomu sökme (13 eV)", a["max_carrier_energy_eV"] / 13.0),
            ("yanal alan / STM atom itme (~1 V/nm)", max(a["max_lateral_field_V_per_nm"]) / 1.0),
            ("ısınma ΔT / ~600 K (katkı difüzyonu)", a["joule_heating_dT_K"]["k=10.0 W/mK"] / 600.0),
            ("akım alanı Zeeman / kT (4 K)", a["electron_zeeman_eV"] / a["kT_4K_eV"]),
            ("sıcak taşıyıcı hasar hızı / nominal 1,8 V", 1.0 / a["hci_slowdown_vs_nominal_1p8V"]),
            ("atom (0,5 nm) / hücre (~40 µm) boyu", 0.5e-9 / 40e-6)]
    y = np.arange(len(rows))
    vals = np.log10([r[1] for r in rows])
    ax.barh(y, vals, color=["C3" if v < 0 else "C2" for v in vals])
    ax.axvline(0, c="k", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xlabel("log10(mevcut / gereken)  — 0'ın solu: yetersiz")
    ax.set_title("(b) atom programlama eşikleri: hiçbiri aşılmıyor\n(3,3 V durumu; 5 V'ta da V_D tepe 0,82 V)", fontsize=10)
    f.suptitle("Exp23: geri dönüştürülen enerjideki kayıp ve atom kontrolü", fontsize=12)
    f.tight_layout(rect=(0, 0, 1, 0.94))
    f.savefig(OUT / "fig_exp23.png", dpi=95)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "fig":
        fig()
    else:
        main()
