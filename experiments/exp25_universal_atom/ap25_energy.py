"""Exp25 (energy manipulation).
E1 energy routing: when a cell dies its power goes to the survivors (exp23: eta ~ 1). Here the destination is chosen by a
   small gate bias on one survivor during a slow kill (exp18 3.3 V state, 3 cells, implicit solver, 30 s ramp of a -300 mV
   band on the rightmost cell). Measured: channel-power gain of each survivor, share of the chosen one, total eta.
E2 energy per atomic operation vs the Landauer bound k T ln 2 (4 K) and vs the cell's own running power.
E3 energy-level manipulation by the cell field (quadratic Stark): dA/A = eta_a E^2 (31P theory -2e-2 um^2/V^2, Sb measured
   -3.7e-3; Bradbury 2006), dg/g = eta_g E^2 (-1e-5); fields of a cell front: 0.9 V over 0.5-2.2 um.
usage: python ap25_energy.py run  ->  results/exp25/energy.json (predictions_energy.json written first)"""
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp18_programmable"))
import regen18 as rg  # noqa: E402
from run18 import law, sp, prog_array, W, DY, LAM, R_D  # noqa: E402

OUT = ROOT / "results" / "exp25"
VH = 3.3
KW = {k: v for k, v in rg.KW.items() if k != "t_relax"}
N = int(round(W / DY))
Y = (np.arange(N) + 0.5) * DY
D = np.random.default_rng(rg.SEED).normal(0.0, 3.356e-3 / np.sqrt(DY * 0.15), N)
KB_EV = 8.617333e-5


def step(V, regs, dt, store=False):
    r = sp.simulate(vdd_up=(VH, VH), dvdd=1.0, prog=prog_array(regs), V0=V, t_relax=dt, n_sub=12, implicit_load=True,
                    store_at=[("up", VH)] if store else (), **KW)
    return np.array(r["V_final"]), r


def cell_powers(r, regs):
    Vb = np.array(r["profiles"][f"up_{VH:.2f}"])
    VD = np.array(r["profiles"][f"up_{VH:.2f}_VD"])
    jd = sp.j_d(Vb, law.BASE["vov0"] + D + prog_array(regs), law.BASE["ss"])
    cells = rg.intervals(r["hist"][0])
    out = []
    for s, e in cells:
        m = (Y >= s) & (Y < e)
        out.append(dict(interval=[float(s), float(e)], P_uW=float(np.sum(VD[m] * jd[m] * DY) * 1e6)))
    return out, float(np.sum(VD * jd * DY) * 1e6)


def scenario(target, bias):
    g = json.loads((ROOT / "results" / "exp18" / "regen" / f"grow_v{VH:g}.json").read_text(encoding="utf-8"))
    V = np.array(g["V_final"])
    V, r = step(V, [], 1.0, store=True)
    before, tot0 = cell_powers(r, [])
    cells = [c["interval"] for c in before]
    dead = cells[-1]
    zone = ((dead[0] + dead[1]) / 2, dead[1] - dead[0] + 60.0)
    steer = [] if target is None else [((cells[target][0] + cells[target][1]) / 2, cells[target][1] - cells[target][0], bias)]
    for k in range(30):
        V, r = step(V, [(zone[0], zone[1], -0.30 * (k + 1) / 30)] + steer, 1.0)
    regs = [(zone[0], zone[1], -0.30)] + steer
    for _ in range(5):
        V, r = step(V, regs, 1.0)
    V, r = step(V, regs, 1.0, store=True)
    after, tot1 = cell_powers(r, regs)
    gains = []
    for c in before[:-1]:
        ov = [a for a in after if min(a["interval"][1], c["interval"][1]) - max(a["interval"][0], c["interval"][0]) > 0]
        gains.append(sum(a["P_uW"] for a in ov) - c["P_uW"])
    born = [a for a in after if not any(min(a["interval"][1], c["interval"][1]) - max(a["interval"][0], c["interval"][0]) > 0
                                        for c in before[:-1])]
    p_dead = before[-1]["P_uW"]
    return dict(target=target, bias_mV=bias * 1e3, before=before, after=after, dead_uW=p_dead, gains_uW=gains,
                share=[gi / max(sum(gains), 1e-12) for gi in gains], eta=(tot1 - (tot0 - p_dead)) / p_dead,
                n_before=len(before), n_after=len(after), births=len(born))


def tables():
    kT = KB_EV * 4.0
    landauer = kT * np.log(2)
    e2 = dict(landauer_4K_eV=landauer,
              charge_state_flip_eV=dict(min_irreversible=8 * kT, note="body voltage stepped across the address by 8 kT/e"),
              cell_running_power_uW=136.0, cell_energy_per_1ms_op_eV=136e-6 * 1e-3 / 1.602e-19,
              ionic_write_1us_eV=dict(n1=1.466 * 30e-6 * 1e-6 / 1.602e-19, per_atom_moved=1.466 * 30e-6 * 1e-6 / 1.602e-19 / 2),
              oxygen_barrier_eV=1.0, nuclear_flip_photon_eV=6.626e-34 * 17.235e6 / 1.602e-19,
              rf_pi_pulse_1mW_50us_eV=1e-3 * 50e-6 / 1.602e-19)
    e3 = {}
    A0, g_ghz_per_T = 117.53e6, 27.97e9
    for wf in (0.5, 1.0, 2.2):
        E = 0.9 / wf                                         # V/um at a cell front of width wf um
        e3[f"front {wf} um ({E:.2f} V/um)"] = dict(
            dA_P_theory_MHz=A0 * 2e-2 * E**2 / 1e6, dA_Sb_like_MHz=A0 * 3.7e-3 * E**2 / 1e6,
            dESR_at_1T_kHz=g_ghz_per_T * 1e-5 * E**2 / 1e3, below_ionisation=bool(E < 3.7))
    return dict(E2=e2, E3=e3)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pf = OUT / "predictions_energy.json"
    if not pf.exists():
        pf.write_text(json.dumps(dict(registered=time.strftime("%Y-%m-%d %H:%M"), predictions=[
            "EN1: without steering each of the two survivors receives >= 20 % of the dead cell's power",
            "EN2: a +30 mV gate bias on one survivor makes it receive >= 80 % of the dead cell's power (both choices)",
            "EN3: the total recycling efficiency stays within 0.95-1.05 and no cell is born or dies in the steered runs",
            "EN4 (+15 mV): the share of the chosen survivor lies between the unsteered and the +30 mV values"]), indent=1),
            encoding="utf-8")
    res = dict(tables=tables(), routing={})
    for target, bias in ((None, 0.0), (0, 0.030), (1, 0.030), (0, 0.015), (1, 0.015)):
        key = "none" if target is None else f"cell{target}_+{bias*1e3:g}mV"
        res["routing"][key] = s = scenario(target, bias)
        print(f"{key:16s} dead {s['dead_uW']:.1f} uW -> gains {[round(x, 1) for x in s['gains_uW']]} uW, share "
              f"{[round(x, 3) for x in s['share']]}, eta {s['eta']:.3f}, cells {s['n_before']}->{s['n_after']}, births {s['births']}",
              flush=True)
    (OUT / "energy.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res["tables"], indent=1))


if __name__ == "__main__":
    main()
