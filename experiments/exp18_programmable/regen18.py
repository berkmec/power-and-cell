"""Exp18 step 2: deliberate damage and regeneration of current-filament "cells" (same device as run18).
After a growth sweep to V_hold the rightmost cell (target) is damaged with a local V_T increase, and the cell population
is followed at constant V_hold:
  control                no damage
  recycle                target dissolved in place, damage kept ("converted to energy": the load line keeps the total
                         current, so its current must go to the healthy cells)
  recycle_heal           damage removed afterwards: do healthy cells re-occupy the region?
  recycle_replace_heal   + programmed division of the neighbour with a weak spot at its centre, then healed
  eject                  target dragged to the device end by a moving gate pattern (+40 mV ahead of the leading front,
                         -40 mV over the trailing edge) and dissolved there in a sink
  eject_replace          eject + programmed division of the neighbour
  fragment(_heal)        (v2) only the central 30 um of the target are damaged: is it cut into two surviving pieces?
Protocols: v1 (offsets -100 mV damage, -60 mV weak spot, +-5 um margins; first attempt), v2 (lethal -300 mV with
30 um margins, -150 mV weak spot, switched on at once; chosen after the v1 diagnosis) and v3 (as v2, but damage, cut band
and sink ramp from 0 to -300 mV over 30 s: slow death, no load-line voltage spike).
usage: python regen18.py grow | predict | predict_v2 | predict_v3 | run [nproc] [v1|v2|v3]  ->  results/exp18/regen*/*.json"""
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run18  # noqa: E402  (sets low priority, paths, prog_array, split_events)
from run18 import law, sp, split_events, prog_array, W, DY, LAM, R_D, R_L  # noqa: E402

BASEDIR = run18.ROOT / "results" / "exp18"
OUT = BASEDIR / "regen"                    # grow states and v1 runs
SEED = 1                                   # exp17 run H_lam60_s1: one filament -> 2 -> 3 -> 4 (known lineage)
VHOLDS = (3.3, 5.0)
VARIANTS = ("control", "recycle", "recycle_heal", "recycle_replace_heal", "eject", "eject_replace")
PROTOCOLS = {
    "v1": dict(damage=-0.10, margin=5.0, weak=-0.06, sink=-0.10, sink_margin=2.0, ramp=0, variants=VARIANTS, out="regen"),
    "v2": dict(damage=-0.30, margin=30.0, weak=-0.15, sink=-0.30, sink_margin=30.0, ramp=0,
               variants=VARIANTS + ("fragment", "fragment_heal"), out="regen_v2"),
    "v3": dict(damage=-0.30, margin=30.0, weak=-0.15, sink=-0.30, sink_margin=30.0, ramp=30,
               variants=VARIANTS + ("fragment", "fragment_heal"), out="regen_v3"),
}
DAMAGE, WEAK, PUSH = -0.10, -0.06, 0.04   # v1 values (used by predict() and the analysis of v1)
KW = dict(W=W, dy=DY, R_sq=1e12, R_L=R_L, r_D=R_D, contacts="distributed", r_v=R_D * LAM**2, down=False,
          sigma_scale=1.0, t_relax=1.0, seed=SEED, **law.BASE)


def intervals(h):
    return sorted((c - w / 2, c + w / 2) for c, w in zip(h["centres_um"], h["widths_um"]))


def track(prev, cur):
    best = max(cur, key=lambda iv: min(iv[1], prev[1]) - max(iv[0], prev[0]), default=None)
    return best if best is not None and min(best[1], prev[1]) - max(best[0], prev[0]) > 0 else None


def step(V, vdd, regions, hist, phase):
    r = sp.simulate(vdd_up=(vdd, vdd), dvdd=1.0, prog=prog_array(regions), V0=V, **KW)
    h = dict(r["hist"][0], phase=phase, regions=[list(map(float, g)) for g in regions])
    ev = sp.lineage(intervals(hist[-1]), intervals(h))          # chained lineage (simulate counts from an empty start)
    h.update(ev)
    hist.append(h)
    return np.array(r["V_final"]), h


def grow(vhold):
    r = sp.simulate(vdd_up=(0.6, vhold), dvdd=0.02, **KW)
    hist = [dict(h, phase="grow", regions=[]) for h in r["hist"]]
    out = dict(vhold=vhold, hist=hist, V_final=r["V_final"])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"grow_v{vhold:g}.json").write_text(json.dumps(out), encoding="utf-8")
    return out


def predict():
    """v1 predictions (written before the v1 runs)."""
    c = law.constants()
    res = {}
    for vh in VHOLDS:
        g = json.loads((OUT / f"grow_v{vh:g}.json").read_text(encoding="utf-8"))
        h = g["hist"][-1]
        cells = intervals(h)
        target, surv = cells[-1], cells[:-1]
        wt = target[1] - target[0]
        dw = wt / len(surv)
        widened = [(max(0.0, (s + e) / 2 - (e - s + dw) / 2), min(W, (s + e) / 2 + (e - s + dw) / 2)) for s, e in surv]
        ratio_after = [law.du_config(widened, iv, W, LAM, R_D, c) / c["dV_c"] for iv in widened]
        nb = surv[-1]
        dV_down_weak = law.constants(dict(law.BASE, vov0=law.BASE["vov0"] + WEAK))["V_down"] - c["V_down"]
        ratio_weak = (law.du_config(cells, nb, W, LAM, R_D, c) + dV_down_weak) / c["dV_c"]
        res[f"{vh:g}"] = dict(
            cells_before=cells, target=target, neighbour=nb, I_uA=h["I_uA"],
            wtot_over_I_jon=sum(e - s for s, e in cells) / (h["I_uA"] * 1e-6 / c["j_on"]),
            P1_wtot_conserved="w_tot/(I/j_on) stays within 0.97-1.05 during damage and recovery",
            survivors_equal_share=widened, survivors_ratio_exact_config=ratio_after,
            P2_spontaneous_split_predicted=[r_ >= 1.0 for r_ in ratio_after],
            P3_no_births="no new nucleation during any hold phase; every new cell comes from a division",
            weak_spot_dV_down_mV=dV_down_weak * 1e3, neighbour_ratio_with_weak_spot=ratio_weak,
            P4_weak_spot_split_predicted=ratio_weak >= 1.0,
            P5_eject="exploratory, no quantitative prediction")
    return res


PREDICTIONS_V2 = dict(
    written_before_v2_runs=True,
    basis="v1 diagnosis: -100 mV leaves V_down' (809 mV) below the drain voltage that the damaged region reaches once the "
          "target shrinks (836 mV); -60 mV weak spot is only ~4 mV past the saddle-node, critical slowing down; "
          "-150 mV quenched at once in a single-step test",
    P1_recycle="lethal damage (-300 mV, V_down' = 2.04 V): target dies within 1 step; with 30 um margins no piece escapes; "
               "its width reappears in the survivors; w_tot/(I/j_on) within 0.97-1.05 once the target is gone",
    P2_spontaneous="3.3 V: survivors do not divide; 5.0 V: the neighbour next to the vacated space divides by itself "
                   "(v1 equal-share calculation)",
    P3_weak_spot="a -150 mV, 6 um weak spot at the neighbour's centre divides it within 1-2 steps at both drives",
    P4_heal="after healing, the vacated region is re-occupied by expansion or migration of neighbouring cells, with no "
            "genuine nucleation (no new cell farther than 10 um from all previous cells)",
    P5_eject="the gate pattern carries the target to the device end (as in v1, 5-9 steps) and the lethal sink dissolves "
             "it there without the cell escaping back",
    P6_fragment="damaging only the central 30 um cuts the target into two pieces that survive outside the damaged band",
)


PREDICTIONS_V3 = dict(
    written_before_v3_runs=True,
    basis="v2 and a 20 ms diagnostic: an abrupt lethal kill releases the target's current at once, the common node jumps "
          "and new cells nucleate at the nearest free healthy places (no escape is possible through a -300 mV band); a "
          "gradual kill (-10 mV/s) produced no nucleation and the survivors widened (3.3 V: middle cell 58 -> 108 um)",
    P1_slow_recycle="slow death: no genuine nucleation (no new cell > 10 um from all earlier cells) during damage; the "
                    "target's width is taken up by widening of existing cells; 3.3 V: N 3 -> 2 and stays; 5.0 V: the "
                    "widened neighbour divides by itself (N back to 3)",
    P2_replace="3.3 V: the -150 mV weak spot divides the widened neighbour (N back to 3); after healing a cell enters the "
               "vacated region",
    P3_eject="the target is carried to the end and dies in the slowly ramped sink without any genuine nucleation",
    P4_fragment="a slowly ramped central band cuts the target into two surviving pieces, without genuine nucleation",
    P5_conservation="w_tot/(I/j_on) within 0.97-1.05 throughout",
)


def scaled(regions, f):
    return [(c, w, dv * f) for c, w, dv in regions]


def damage_steps(V, vh, regs, hist, phase, n, ramp, extra=()):
    """n steps with the regions regs ramped linearly over `ramp` steps (ramp = 0: full strength at once)."""
    for k in range(n):
        f = min(1.0, (k + 1) / ramp) if ramp else 1.0
        V, h = step(V, vh, scaled(regs, f) + list(extra), hist, phase)
    return V, h


def run_variant(job):
    t0 = time.time()
    vh, var, pr = job["vhold"], job["variant"], PROTOCOLS[job["protocol"]]
    g = json.loads((OUT / f"grow_v{vh:g}.json").read_text(encoding="utf-8"))
    hist = [g["hist"][-1]]
    V = np.array(g["V_final"])
    cells = intervals(hist[-1])
    target, nb = cells[-1], cells[-2]
    dmg = [((target[0] + target[1]) / 2, target[1] - target[0] + 2 * pr["margin"], pr["damage"])]
    log = dict(target=target, neighbour=nb, protocol=job["protocol"])

    def weak_on(nb_now):
        return [((nb_now[0] + nb_now[1]) / 2, 6.0, pr["weak"])] if nb_now else []

    if var == "control":
        for _ in range(45):
            V, h = step(V, vh, [], hist, "hold")
    elif var.startswith("recycle"):
        V, h = damage_steps(V, vh, dmg, hist, "damage", max(15, pr["ramp"]), pr["ramp"])
        if var == "recycle":
            V, h = damage_steps(V, vh, dmg, hist, "damage", 30, 0)
        elif var == "recycle_heal":
            for _ in range(30):
                V, h = step(V, vh, [], hist, "heal")
        else:
            weak = weak_on(track(nb, intervals(h)))
            log["weak_spot"] = weak
            V, h = damage_steps(V, vh, dmg, hist, "replace", 8, 0, extra=weak)
            for _ in range(22):
                V, h = step(V, vh, [], hist, "heal")
    elif var.startswith("fragment"):
        cut = [((target[0] + target[1]) / 2, 30.0, pr["damage"])]
        log["cut_band"] = cut
        V, h = damage_steps(V, vh, cut, hist, "damage", max(15, pr["ramp"]), pr["ramp"])
        if var == "fragment":
            V, h = damage_steps(V, vh, cut, hist, "damage", 30, 0)
        else:
            for _ in range(30):
                V, h = step(V, vh, [], hist, "heal")
    else:
        cur, k = target, 0
        while cur is not None and cur[1] < W - 1.0 and k < 80:
            regs = [(cur[1] + 7.5, 15.0, PUSH), (cur[0] + 6.0, 12.0, -PUSH)]
            V, h = step(V, vh, regs, hist, "eject")
            cur, k = track(cur, intervals(h)), k + 1
        log.update(eject_steps=k, eject_reached_end=bool(cur is not None and cur[1] >= W - 1.0), target_at_end=cur)
        sink = ([((cur[0] - pr["sink_margin"] + W) / 2, W - cur[0] + pr["sink_margin"], pr["sink"])] if cur else [])
        log["sink"] = sink
        V, h = damage_steps(V, vh, sink, hist, "sink", max(10, pr["ramp"]), pr["ramp"])
        if var == "eject_replace":
            weak = weak_on(track(nb, intervals(h)))
            log["weak_spot"] = weak
            for _ in range(8):
                V, h = step(V, vh, weak, hist, "replace")
        for _ in range(20 if var == "eject_replace" else 28):
            V, h = step(V, vh, [], hist, "heal")
    res = dict(job=job, log=log, hist=hist, splits=split_events(hist), V_final=V.tolist(), wall=time.time() - t0)
    out = BASEDIR / pr["out"]
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{job['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    n_seq = [x["n_fil"] for x in hist]
    return (f"{job['name']}: N {n_seq[0]} -> min {min(n_seq)} -> son {n_seq[-1]} | bölünme {len(res['splits'])} "
            f"({res['wall']:.0f} s)")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "grow":
        for vh in VHOLDS:
            g = grow(vh)
            h = g["hist"][-1]
            print(f"V_hold={vh}: N={h['n_fil']} hücreler={[(round(s), round(e)) for s, e in intervals(h)]} I={h['I_uA']:.0f} uA")
    elif mode in ("predict", "predict_v2", "predict_v3"):
        path = BASEDIR / {"predict": "predictions_regen.json", "predict_v2": "predictions_regen_v2.json",
                          "predict_v3": "predictions_regen_v3.json"}[mode]
        if path.exists():
            print("zaten var, üzerine yazılmadı:", path)
        else:
            pr = {"predict": predict, "predict_v2": lambda: PREDICTIONS_V2, "predict_v3": lambda: PREDICTIONS_V3}[mode]()
            path.write_text(json.dumps(pr, indent=1, ensure_ascii=False), encoding="utf-8")
            print("yazıldı:", path)
    else:
        proto = sys.argv[3] if len(sys.argv) > 3 else "v1"
        pr = PROTOCOLS[proto]
        jobs = [dict(name=f"regen_{v}_v{vh:g}", vhold=vh, variant=v, protocol=proto) for vh in VHOLDS for v in pr["variants"]]
        jobs = [j for j in jobs if not (BASEDIR / pr["out"] / f"{j['name']}.json").exists()]
        with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 3) as p:
            for line in p.imap_unordered(run_variant, jobs):
                print(line, flush=True)
        print("TAMAM", flush=True)
