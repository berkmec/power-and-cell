"""Exp20: ring (periodic) geometry and common-node capacitance in the exp12 1D filament solver."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(not (ROOT / "experiments" / "exp12_filament" / "spice" / "tab_vb0.0.dat").exists(),
                                reason="ngspice tables not generated")
sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))


def test_periodic_rail_is_symmetric_around_a_point_source():
    import spatial as sp
    jd = np.zeros(200)
    jd[0] = 1e-6
    V = sp.rail_voltage(jd, 1.0, 10.0, 0.5, "distributed", 10.0 * 60.0**2, periodic=True)
    assert V[1] == pytest.approx(V[-1], rel=1e-12) and V[5] == pytest.approx(V[-5], rel=1e-12)
    assert np.argmin(V) == 0


def test_ring_merges_a_cell_across_the_seam():
    import spatial as sp
    y = (np.arange(80) + 0.5) * 0.5
    prog = np.where((y < 2.0) | (y > 38.0), 0.05, 0.0)          # seed straddling the seam of a 40 um ring
    r = sp.simulate(W=40.0, dy=0.5, R_sq=1e12, R_L=5e3, vdd_up=(0.6, 1.4), dvdd=0.05, mismatch=False, down=False,
                    prog=prog, periodic=True, r_D=10.0, contacts="distributed", r_v=10.0 * 60.0**2)
    b = next(h for h in r["hist"] if h["n_fil"] > 0)
    assert b["n_fil"] == 1 and (b["centres_um"][0] < 5.0 or b["centres_um"][0] > 35.0)


def test_node_capacitance_delays_the_common_node():
    import spatial as sp
    kw = dict(W=20.0, dy=0.5, R_sq=1e12, R_L=5e3, dvdd=1.0, down=False, seed=2, V0=np.zeros(40))
    fast = sp.simulate(vdd_up=(1.0, 1.0), t_relax=1e-3, n_sub=2, **kw)
    slow = sp.simulate(vdd_up=(1.0, 1.0), t_relax=1e-3, n_sub=2, C_c=1e-6, vc0=0.6, **kw)   # R_L C = 5 ms
    assert fast["hist"][0]["vds"] > 0.99
    assert 0.6 < slow["hist"][0]["vds"] < 0.7                       # 1 ms of charging toward 1.0 V
    assert sp.simulate(vdd_up=(1.0, 1.0), t_relax=1e-3, n_sub=2, **kw)["hist"] == fast["hist"]
