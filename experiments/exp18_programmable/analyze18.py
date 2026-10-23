"""Exp18 analysis: outcomes against the predictions written before the runs, plus figures.
step 1 (place, wta, two): first-birth site vs. the predicted site (argmax of mismatch + programmed offset);
step 2 (regen): cell number, deaths, divisions (with the exact-configuration mechanism ratio), conservation of the total
filament width, re-occupation of the damaged region.
usage: python analyze18.py   ->  results/exp18/analysis.json, fig_exp18_step1.png, fig_exp18_regen.png"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import regen18 as rg  # noqa: E402
import run18  # noqa: E402
from run18 import law, W, LAM, R_D  # noqa: E402

RES = run18.ROOT / "results" / "exp18"


def load(pattern, sub="runs"):
    return {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted((RES / sub).glob(pattern))}


def born(r):
    b = r["birth"]
    return [(c - w / 2, c + w / 2) for c, w in zip(b["centres_um"], b["widths_um"])] if b else []


def overlaps(ivs, c, w, tol=2.0):
    return any(min(e, c + w / 2 + tol) - max(s, c - w / 2 - tol) > 0 for s, e in ivs)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    m = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (float(m - h), float(m + h))


def step1():
    out = {}
    pr = {m: json.loads((RES / f"predictions_{m}.json").read_text(encoding="utf-8")) for m in ("place", "wta", "two")}
    site_hits = site_n = 0
    for m in ("place", "wta", "two"):
        for name, r in load(f"{m}_*.json").items():
            yp = pr[m]["per_run"][name]["y_pred"]
            site_hits += any(s - 1 <= yp <= e + 1 for s, e in born(r))
            site_n += 1
    out["birth_site_prediction"] = dict(hits=site_hits, n=site_n)
    place = load("place_*.json")
    curve = {}
    for dv in sorted({r["job"]["regions"][0][2] for r in place.values()}):
        rs = [r for r in place.values() if r["job"]["regions"][0][2] == dv]
        k = sum(overlaps(born(r), 150.0, 10.0) for r in rs)
        only = sum(overlaps(born(r), 150.0, 10.0) and len(born(r)) == 1 for r in rs)
        curve[f"{dv:.3f}"] = dict(n=len(rs), in_seed=k, only_seed=only, ci95=wilson(k, len(rs)),
                                  predicted_MC=pr["place"]["P_birth_in_seed_MC"][f"{dv:.3f}"],
                                  predicted_per_run=sum(pr["place"]["per_run"][n_]["region"] == 0 for n_ in
                                                        [r["job"]["name"] for r in rs]))
    out["place"] = curve
    wta = {}
    for name, r in load("wta_*.json").items():
        regs = r["job"]["regions"]
        won = [k for k, (c, w, dv) in enumerate(regs) if overlaps(born(r), c, w)]
        strongest = pr["wta"]["strongest_input"][name]
        wta[name] = dict(inputs_mV=[round(dv * 1e3, 1) for c, w, dv in regs], strongest=strongest, won=won,
                         n_births=len(born(r)), predicted=pr["wta"]["per_run"][name]["region"],
                         strongest_won=strongest in won, prediction_hit=pr["wta"]["per_run"][name]["region"] in won)
    out["wta"] = wta
    two = {}
    for name, r in load("two_*.json").items():
        two[name] = dict(seed0=overlaps(born(r), 150.0, 10.0), seed1=overlaps(born(r), 450.0, 10.0),
                         n_births=len(born(r)), predicted_first=pr["two"]["per_run"][name]["region"])
    out["two"] = two
    return out


def ivs(h):
    return rg.intervals(h)


def regen(sub):
    c = law.constants()
    out = {}
    for name, r in load("regen_*.json", sub=sub).items():
        ws = r["log"].get("weak_spot") or []
        dvw = (law.constants(dict(law.BASE, vov0=law.BASE["vov0"] + ws[0][2]))["V_down"] - c["V_down"]) if ws else 0.0
        hist = r["hist"]
        target = tuple(r["log"]["target"])
        n_seq = [h["n_fil"] for h in hist]
        wt = [sum(e - s for s, e in ivs(h)) / (h["I_uA"] * 1e-6 / c["j_on"]) for h in hist]
        splits, deaths, births, genuine = [], [], 0, []
        cur_t, t_death = target, None
        for k in range(1, len(hist)):
            prev, now = ivs(hist[k - 1]), ivs(hist[k])
            births += hist[k].get("births", 0)
            genuine += [dict(step=k, phase=hist[k]["phase"], interval=q) for q in now
                        if all(q[1] < a - 10 or q[0] > b + 10 for a, b in prev)]
            if cur_t is not None:
                cur_t = rg.track(cur_t, now)
                if cur_t is None:
                    t_death = k
            for p in prev:
                kids = [q for q in now if min(p[1], q[1]) - max(p[0], q[0]) > 0]
                if len(kids) >= 2:
                    cut = 0.5 * (kids[0][1] + kids[1][0])
                    weak = [g for g in hist[k]["regions"] if g[1] == 6.0 and kids[0][1] - 3 <= g[0] <= kids[1][0] + 3]
                    du = law.du_config(prev, p, W, LAM, R_D, c)
                    splits.append(dict(step=k, phase=hist[k]["phase"], parent=p, cut=cut, kids_w=[q[1] - q[0] for q in kids],
                                       ratio=du / c["dV_c"], at_weak_spot=bool(weak),
                                       ratio_with_weak=(du + dvw) / c["dV_c"] if weak else None))
                if not kids and hist[k]["phase"] != "grow":
                    deaths.append(dict(step=k, phase=hist[k]["phase"], interval=p))
        final = ivs(hist[-1])
        out[name] = dict(n_start=n_seq[0], n_min=min(n_seq), n_end=n_seq[-1], births_after_grow=births,
                         genuine_nucleations=genuine,
                         wtot_ratio_range=[float(min(wt)), float(max(wt))], target_death_step=t_death,
                         splits=splits, deaths=deaths, final_cells=final,
                         target_region_reoccupied=any(min(e, target[1]) - max(s, target[0]) > 0 for s, e in final),
                         eject=dict((k_, v) for k_, v in r["log"].items() if k_.startswith("eject") or k_ == "target_at_end"))
    return out


def fig_step1(a):
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    ax = axs[0]
    dvs = sorted(float(k) for k in a["place"])
    obs = [a["place"][f"{d:.3f}"]["in_seed"] / a["place"][f"{d:.3f}"]["n"] for d in dvs]
    ci = [a["place"][f"{d:.3f}"]["ci95"] for d in dvs]
    ax.errorbar(np.array(dvs) * 1e3, obs, yerr=[[o - l for o, (l, h) in zip(obs, ci)], [h - o for o, (l, h) in zip(obs, ci)]],
                fmt="o", c="C0", ms=8, capsize=4, label="benzetim (10 rastgele alan / nokta, %95 aralık)")
    ax.plot(np.array(dvs) * 1e3, [a["place"][f"{d:.3f}"]["predicted_MC"] for d in dvs], "k--", lw=2,
            label="önceden kaydedilen öngörü (en büyük aşırı sürüşlü hücre doğar)")
    ax.set_xlabel("tohum: 10 µm'lik bölgede kapı aşırı sürüşü Δ (mV)")
    ax.set_ylabel("ilk hücrenin tohumda doğma olasılığı")
    ax.set_title("(a) Doğum yerini programlamak (600 µm aygıt, rastgele V_T farkı σ = 12 mV / 0,5 µm)", fontsize=10)
    ax.legend(fontsize=8)
    ax.set_ylim(-0.05, 1.05)
    ax = axs[1]
    rows = sorted(a["wta"].items())
    for i, (name, v) in enumerate(rows):
        for k, mv in enumerate(v["inputs_mV"]):
            ax.add_patch(plt.Rectangle((k - 0.4, i - 0.4), 0.8, 0.8, color=plt.cm.Blues(0.15 + mv / 40)))
            ax.text(k, i, f"{mv:.0f}", ha="center", va="center", fontsize=7)
        for k in v["won"]:
            ax.add_patch(plt.Rectangle((k - 0.45, i - 0.45), 0.9, 0.9, fill=False, ec="C3", lw=2.5))
        ax.text(6.0, i, "✓" if v["strongest_won"] else "✗", fontsize=11, va="center",
                color="C2" if v["strongest_won"] else "C3")
    ax.set_xlim(-0.6, 6.5)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xticks(range(6))
    ax.set_xticklabels([f"giriş {k + 1}\n{int(y)} µm" for k, y in enumerate(run18.INPUTS)], fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([n for n, _ in rows], fontsize=7)
    ax.set_title("(b) Kazanan hepsini alır: 6 giriş (mV), kırmızı çerçeve = hücrenin doğduğu giriş;\n"
                 "✓ en güçlü giriş kazandı", fontsize=10)
    fig.tight_layout()
    fig.savefig(RES / "fig_exp18_step1.png", dpi=95)


def fig_regen(sub, title):
    runs = load("regen_*.json", sub=sub)
    variants = [v for v in ("control", "recycle", "recycle_heal", "recycle_replace_heal", "eject", "eject_replace",
                            "fragment", "fragment_heal") if any(k.startswith(f"regen_{v}_v") for k in runs)]
    names = [f"regen_{v}_v{vh:g}" for vh in rg.VHOLDS for v in variants]
    fig, axs = plt.subplots(2, len(variants), figsize=(4 * len(variants), 9), sharey=True)
    colors = dict(hold="0.5", damage="C3", heal="C2", replace="C4", eject="C1", sink="C3")
    for ax, name in zip(axs.ravel(), names):
        r = runs.get(name)
        if not r:
            ax.set_visible(False)
            continue
        for k, h in enumerate(r["hist"]):
            for s, e in rg.intervals(h):
                ax.plot([k, k], [s, e], c=colors.get(h["phase"], "C0"), lw=2.2)
            for g in h["regions"]:
                ax.add_patch(plt.Rectangle((k - 0.5, g[0] - g[1] / 2), 1.0, g[1],
                                           color="C3" if g[2] < 0 else "C1", alpha=0.18, lw=0))
        for k in range(1, len(r["hist"])):
            prev, now = rg.intervals(r["hist"][k - 1]), rg.intervals(r["hist"][k])
            for p in prev:
                kids = [q for q in now if min(p[1], q[1]) - max(p[0], q[0]) > 0]
                if len(kids) >= 2:
                    ax.plot(k, 0.5 * (kids[0][1] + kids[1][0]), "k*", ms=11)
        ax.set_title(name.replace("regen_", "").replace("_v", "  V_hold="), fontsize=10)
        ax.set_xlabel("adım (her biri 1 s, sabit V_DD)")
    axs[0, 0].set_ylabel("konum y (µm)")
    axs[1, 0].set_ylabel("konum y (µm)")
    fig.suptitle(title + "  Renk evresi: gri bekleme, kırmızı hasar, turuncu dışarı atma, mor programlı bölünme, yeşil "
                 "iyileşme; kırmızı gölge = yerel V_T artışı, turuncu gölge = taşıyıcı desen; ★ bölünme", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(RES / f"fig_exp18_{sub}.png", dpi=75)


if __name__ == "__main__":
    a = {}
    if (RES / "predictions_place.json").exists() and load("place_*.json"):
        a["step1"] = step1()
        fig_step1(a["step1"])
    titles = {"regen": "Exp18 v1: hasar -100 mV (ölümcül değil), zayıf nokta -60 mV.",
              "regen_v2": "Exp18 v2: ani ölümcül hasar -300 mV, zayıf nokta -150 mV.",
              "regen_v3": "Exp18 v3: yavaş ölüm (hasar 30 s'de 0'dan -300 mV'a), zayıf nokta -150 mV."}
    for sub, title in titles.items():
        if (RES / sub).exists() and load("regen_*.json", sub=sub):
            a[sub] = regen(sub)
            fig_regen(sub, title)
    (RES / "analysis.json").write_text(json.dumps(a, indent=1, default=float), encoding="utf-8")
    s1 = a.get("step1")
    if s1:
        print("doğum yeri öngörüsü:", s1["birth_site_prediction"])
        for k, v in s1["place"].items():
            print(f"  Δ={float(k)*1e3:4.0f} mV: tohumda {v['in_seed']}/{v['n']} (yalnız tohumda {v['only_seed']}), "
                  f"öngörü MC {v['predicted_MC']:.2f}, koşu bazında {v['predicted_per_run']}/{v['n']}")
        print("  WTA: en güçlü kazandı", sum(v["strongest_won"] for v in s1["wta"].values()), "/", len(s1["wta"]),
              "| öngörü tuttu", sum(v["prediction_hit"] for v in s1["wta"].values()))
        print("  iki tohum: ikisi birden", sum(v["seed0"] and v["seed1"] for v in s1["two"].values()), "/", len(s1["two"]),
              "| en az biri", sum(v["seed0"] or v["seed1"] for v in s1["two"].values()))
    for sub in titles:
      print("=====", sub)
      for name, v in a.get(sub, {}).items():
        print(f"{name}: N {v['n_start']}->{v['n_min']}->{v['n_end']} gerçek doğum={len(v['genuine_nucleations'])} "
              f"wtot={v['wtot_ratio_range'][0]:.3f}-{v['wtot_ratio_range'][1]:.3f} ölüm adımı={v['target_death_step']} "
              f"yer dolu={v['target_region_reoccupied']} bölünmeler="
              f"{[(s['phase'], round(s['ratio'], 2), s['at_weak_spot']) for s in v['splits']]} {v['eject'] or ''}")
