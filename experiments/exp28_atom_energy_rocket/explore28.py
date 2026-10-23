"""Exp28 exploration (beyond the known thruster regimes, same physics): doors opened by the project data.
 nanogap : cell voltages (0.8-5 V) over SKY130-scale gaps (4.148 nm oxide ... 10 um): classical Child-Langmuir thrust density
           (8/9) eps0 (V/d)^2, exhaust velocity sqrt(2 q V / m), efficiency eta = eta_b K / (K + E_entry + V_M) where E_entry is the
           ionisation energy (neutral propellant) or only the Exp24 escape barrier Q0 (pre-charged ions of an ionic solid/liquid)
 impulse : atom-by-atom impulse bit m v_e (artificial-atom turnstile, Exp0-4)
 breath  : ambient atoms as propellant (the Exp24 species, oxygen, is the main gas near 200 km): net thrust = drag needs
           v_e = v (1 + C_D / (2 eps_in f_m)) and P/A = eps_in f_m rho v v_e^2 / (2 eta) -> no onboard mass ratio
 beamed  : external power converted at eta_pv and spent on atoms: F / P_beam = eta_pv 2 eta / v_e versus the sail's 2/c
 scaling : the cells turn extra power into more cells (Exp23): thrust in equal one-cell steps (Exp19: k = 1, 2, 4; Exp20: 3 -> 7)
usage: python explore28.py  -> results/exp28/explore.json"""
import json

import numpy as np

from aer28 import AMU, C, E, EPS0, OUT, PROP, Q_OX, ALPHA, V_M, P_CELL, V_D, ve

SPECIES = {"In+ (nötr yakıt, iyonlaşma 5,79 eV)": (114.818, 5.786), "Cs+ (nötr, 3,89 eV)": (132.905, 3.894),
           "EMI+ (iyonik sıvı, önceden yüklü, Q0 1 eV)": (111.17, Q_OX), "O- (iyonik katı, Exp24 bariyeri 1 eV)": (15.999, Q_OX)}


def nanogap():
    rows = []
    for d_nm in (4.148, 10.0, 100.0, 210.0, 1000.0, 10000.0):
        for V in (0.8, 3.3, 5.0):
            Ef = V / (d_nm * 1e-9)
            row = dict(gap_nm=d_nm, V=V, field_V_nm=V / d_nm, CL_thrust_density_N_m2=8 / 9 * EPS0 * Ef**2,
                       emission_ok_exp24_barrier=bool(V / d_nm >= Q_OX / ALPHA), emission_ok_In_evap=bool(V / d_nm >= 8.26))
            sp = {}
            for name, (m, e_in) in SPECIES.items():
                v = np.sqrt(2 * E * V / (m * AMU))
                eta = 0.9 * V / (V + e_in + V_M)
                sp[name] = dict(v_e_km_s=v / 1e3, Isp_s=v / 9.80665, eta=eta, F_over_P_uN_W=2 * eta / v * 1e6)
            row["species"] = sp
            rows.append(row)
    return rows


def impulse():
    m = PROP["In"]["m"] * AMU
    bit = m * ve("In", 5e3)
    return dict(In_5kV_impulse_bit_Ns=bit, ratio_to_10uNs=1e-5 / bit, O_5V_bit_Ns=16 * AMU * np.sqrt(2 * E * 5 / (16 * AMU)))


def breath():
    v = np.sqrt(3.986004418e14 / (6.371e6 + 200e3))
    nO = 2.5e15                                    # m^-3 at 200 km (Atmosphere Explorer-C, measured)
    rho_O = nO * 15.999 * AMU
    out = dict(v_orbit_m_s=v, n_O_m3=nO, rho_O_kg_m3=rho_O)
    rows = []
    for rho_tot in (1.5e-10, 2.0e-10, 3.0e-10):    # total density: model assumption (not verified here)
        for label, fm in (("bütün türler", 1.0), ("yalnız O", rho_O / rho_tot)):
            for eps in (0.4,):
                for eta in (0.6,):
                    CD = 2.2
                    v_e = v * (1 + CD / (2 * eps * fm))
                    P_A = eps * fm * rho_tot * v * v_e**2 / (2 * eta)
                    flux = eps * fm * rho_tot * v / (15.999 * AMU if fm < 1 else 24.0 * AMU)
                    P_ion = flux * 14.0 * E
                    rows.append(dict(rho_total=rho_tot, collected=label, mass_fraction=fm, intake=eps, eta=eta, v_e_km_s=v_e / 1e3,
                                     power_W_per_m2=P_A, ionisation_W_per_m2=P_ion))
    out["drag_compensation"] = rows
    out["law"] = "v_e = v (1 + C_D / (2 eps_in f_m)); independent of density when all species are collected"
    return out


def beamed():
    rows = []
    for v in (30e3, 60e3):
        for eta in (0.6,):
            f = 0.5 * 2 * eta / v
            rows.append(dict(v_e_km_s=v / 1e3, N_per_beamed_W=f, ratio_to_sail=f / (2 / C)))
    return rows


def scaling():
    I = [p / V_D for p in P_CELL]
    m = PROP["In"]["m"] * AMU / E
    rows = []
    for k in (1, 2, 3, 4, 7):
        F = k * 0.9 * I[1] * m * ve("In", 5e3)
        P = k * I[1] * (5e3 + V_M + PROP["In"]["I"])
        rows.append(dict(cells=k, thrust_uN=F * 1e6, power_W=P, F_over_P_uN_W=F / P * 1e6))
    spread = max(r["F_over_P_uN_W"] for r in rows) / min(r["F_over_P_uN_W"] for r in rows) - 1
    return dict(steps=rows, F_over_P_spread=spread, note="F/P per cell = (m/q) v_e / V_total is independent of the cell current")


if __name__ == "__main__" and "plot" not in __import__("sys").argv:
    res = dict(nanogap=nanogap(), impulse=impulse(), breath=breath(), beamed=beamed(), scaling=scaling())
    (OUT / "explore.json").write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float), encoding="utf-8")
    ng = [r for r in res["nanogap"] if r["V"] == 5.0]
    for r in ng:
        best = max(r["species"].items(), key=lambda kv: kv[1]["F_over_P_uN_W"])
        print(f"d {r['gap_nm']:>8} nm: F/A {r['CL_thrust_density_N_m2']:.3g} N/m2, field {r['field_V_nm']:.3g} V/nm, emission(Exp24) "
              f"{r['emission_ok_exp24_barrier']}, In {r['emission_ok_In_evap']}, best {best[0][:14]} {best[1]['F_over_P_uN_W']:.0f} uN/W "
              f"Isp {best[1]['Isp_s']:.0f} s")
    for name, s in ng[0]["species"].items():
        print(f"   5 V {name}: v_e {s['v_e_km_s']:.2f} km/s, eta {s['eta']:.2f}, F/P {s['F_over_P_uN_W']:.0f} uN/W")
    print("impulse", res["impulse"])
    for r in res["breath"]["drag_compensation"]:
        print("breath", r)
    print("beamed", res["beamed"])
    print("scaling spread", res["scaling"]["F_over_P_spread"])


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    X = json.loads((OUT / "explore.json").read_text(encoding="utf-8"))
    R = json.loads((OUT / "rocket.json").read_text(encoding="utf-8"))
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    a = ax[0]
    d = np.logspace(np.log10(1.0), 4, 200)
    for V in (0.8, 3.3, 5.0):
        a.loglog(d, 8 / 9 * EPS0 * (V / (d * 1e-9))**2, label=f"{V} V (hücre gerilimi)")
        a.axvline(V * ALPHA / Q_OX, ls=":", lw=0.8, color="gray")
    a.axhline(1.3, color="k", ls="--", label="NSTAR (ölçülen, 1,3 N/m²)")
    a.axvline(4.148, color="r", lw=0.8, label="SKY130 kapı oksidi 4,148 nm")
    a.set(xlabel="aralık d (nm)", ylabel="uzay yükü sınırı itki yoğunluğu (N/m²)",
          title="(g) nano aralık: (8/9) ε₀ (V/d)²; noktalı: Exp24 yayma eşiği")
    a.legend(fontsize=7)
    a = ax[1]
    for name in X["nanogap"][0]["species"]:
        pts = [(r["species"][name]["Isp_s"], r["species"][name]["F_over_P_uN_W"]) for r in X["nanogap"] if r["gap_nm"] == 100.0]
        a.loglog([p[0] for p in pts], [p[1] for p in pts], "o-", label=name)
    for m in R["measured"]:
        a.loglog(m["Isp_s"], m["F_over_P_uN_W"], "k^", ms=6)
        a.annotate(m["name"].split()[0], (m["Isp_s"], m["F_over_P_uN_W"]), fontsize=6)
    isp = np.logspace(2, 4, 50)
    a.loglog(isp, 2 * 0.6 / (isp * 9.80665) * 1e6, "k--", lw=0.8, label="2η/v_e, η = 0,6")
    a.axhline(1e6 / C, color="r", lw=0.8, label="foton sınırı (Exp27)")
    a.set(xlabel="özgül itki Isp (s)", ylabel="itki/güç (µN/W)", title="(h) düşük gerilim (0,8–5 V) ile ölçülen itkiçler")
    a.legend(fontsize=6)
    a = ax[2]
    rows = X["breath"]["drag_compensation"]
    for lab in ("bütün türler", "yalnız O"):
        r = [x for x in rows if x["collected"] == lab]
        a.plot([x["rho_total"] for x in r], [x["power_W_per_m2"] for x in r], "o-", label=f"{lab} (v_e {r[1]['v_e_km_s']:.0f} km/s @2e-10)")
    a.axhline(300, color="gray", ls=":", label="tipik güneş paneli ~300 W/m²")
    a.set(xlabel="200 km toplam yoğunluk (kg/m³, varsayım)", ylabel="sürüklenmeyi karşılayan güç (W/m² ağız)",
          title="(i) ortam atomlarıyla itki: kütle oranı sınırı yok")
    a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_exp28_explore.png", dpi=125)
    print("saved", OUT / "fig_exp28_explore.png")


if __name__ == "__main__" and "plot" in __import__("sys").argv:
    plot()
