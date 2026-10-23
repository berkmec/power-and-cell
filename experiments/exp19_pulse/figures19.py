"""Exp19 figure (results/exp19/fig_exp19.png) from results/exp19/analysis.json and the run logs."""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run19 as r19  # noqa: E402

a = json.loads((r19.OUT / "analysis.json").read_text(encoding="utf-8"))
fig, axs = plt.subplots(1, 4, figsize=(22, 5.2))

ax = axs[0]
for vh, col, mk, ls, ms in ((3.3, "C0", "o", "-", 11), (5.0, "C3", "s", "--", 6)):
    trs = [tr for tr in r19.TRS if f"{vh:g}|{tr:g}" in a["killspeed"]]
    x = [max(tr, 0.03) for tr in trs]
    ax.plot(x, [a["killspeed"][f"{vh:g}|{tr:g}"]["births"] for tr in trs], marker=mk, ls=ls, c=col, ms=ms,
            mfc="none" if vh == 3.3 else col, label=f"V_DD = {vh:g} V")
ax.set_xscale("log")
ax.set_xlabel("öldürme süresi T_r (s; en soldaki nokta = ani)")
ax.set_ylabel("yeni doğan hücre sayısı")
ax.set_title("(a) Ölüm hızı: ani ölüm yeni hücre doğurur,\nyavaş ölümde akım var olan hücrelere gider", fontsize=10)
ax.legend(fontsize=8)

ax = axs[1]
for vh, col in ((3.3, "C0"), (5.0, "C3")):
    for tr, ls in ((0.0, "-"), (30.0, ":")):
        r = json.loads((r19.OUT / "runs" / f"kill_v{vh:g}_Tr{tr:g}.json").read_text(encoding="utf-8"))
        t = [x["t"] for x in r["log"]]
        ax.plot(t, [x["vc"] for x in r["log"]], ls, c=col, lw=1.5, label=f"{vh:g} V, T_r = {tr:g} s")
ax.set_xscale("log")
ax.set_xlabel("zaman (s)")
ax.set_ylabel("ortak düğüm gerilimi V_c (V, adım sonunda)")
ax.set_title("(b) V_c: ani ölümde sıçrama adım içinde yeni hücrelerce\nsöndürülüyor; adım sonu değerleri", fontsize=10)
ax.legend(fontsize=7)

ax = axs[2]
ks = sorted(a["rsq"], key=float)
ax.bar(range(len(ks)), [a["rsq"][k]["persistent_new"] for k in ks], color="C2")
for i, k in enumerate(ks):
    ax.text(i, a["rsq"][k]["persistent_new"] + 0.08, f"önce {len(a['rsq'][k]['cells_before'])} hücre\n"
            f"geçici olay {a['rsq'][k]['births']}", ha="center", fontsize=8)
ax.set_xticks(range(len(ks)))
ax.set_xticklabels([f"R_□ = {k}" for k in ks], fontsize=8)
ax.set_ylim(0, 3.2)
ax.set_ylabel("ani ölümden sonra kalıcı yeni hücre")
ax.set_title("(c) Gövde bağlaşımı (R_□), ani ölüm, 3,3 V: kalıcı yeni hücre;\nbüyüme düzeni R_□ ile değiştiği için kıyas sınırlı", fontsize=10)

ax = axs[3]
vs = [v for v in r19.VSTEPS if f"{v:g}" in a["step"]]
ax.plot(vs, [a["step"][f"{v:g}"]["mean_seeds"] for v in vs], "o-", c="C0", ms=7, label="hücre olan tohum sayısı (ort., 4 tohum)")
ax.plot(vs, [a["step"][f"{v:g}"]["mean_off"] for v in vs], "s--", c="C3", ms=6, label="tohum dışında doğan hücre (ort.)")
if "slow" in a["step"]:
    ax.axhline(a["step"]["slow"]["mean_seeds"], c="0.4", ls=":", label="yavaş tarama: ilk doğumda hücreli tohum (ort.)")
ax.set_xlabel("ani V_DD adımı V_s (V; dinlenmeden)")
ax.set_ylabel("sayı")
ax.set_title("(d) Ani adımla çok hücre: 4 eşit tohum (+50 mV)\n10 rastgele alan", fontsize=10)
ax.legend(fontsize=7)
fig.suptitle("Exp19: ani akım değişimiyle çoğalan hücreler (model, 4 K SKY130 nFET filamanları)", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(r19.OUT / "fig_exp19.png", dpi=90)
print("ok", r19.OUT / "fig_exp19.png")
