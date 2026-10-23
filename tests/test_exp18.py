"""Exp18: programmable gate offset (prog) in the exp12 1D filament solver."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(not (ROOT / "experiments" / "exp12_filament" / "spice" / "tab_vb0.0.dat").exists(),
                                reason="ngspice tables not generated")
sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))


def test_prog_none_and_zero_are_identical():
    import spatial as sp
    kw = dict(W=20.0, dy=0.5, R_sq=1e12, R_L=5e3, vdd_up=(0.6, 1.0), dvdd=0.1, down=False, seed=3)
    assert sp.simulate(**kw)["hist"] == sp.simulate(prog=np.zeros(40), **kw)["hist"]


def test_continuing_from_a_converged_state_keeps_it():
    import spatial as sp
    kw = dict(W=40.0, dy=0.5, R_sq=1e12, R_L=5e3, dvdd=0.05, down=False, seed=4)
    r1 = sp.simulate(vdd_up=(0.6, 1.6), **kw)
    r2 = sp.simulate(vdd_up=(1.6, 1.6), V0=r1["V_final"], **kw)
    assert r2["hist"][0]["n_fil"] == r1["hist"][-1]["n_fil"]
    assert np.max(np.abs(np.array(r2["V_final"]) - np.array(r1["V_final"]))) < 1e-3


def test_seed_region_sets_the_birth_site():
    import spatial as sp
    y = (np.arange(80) + 0.5) * 0.5
    prog = np.where(np.abs(y - 12.0) < 2.0, 0.03, 0.0)            # uniform device, +30 mV overdrive around y = 12 um
    r = sp.simulate(W=40.0, dy=0.5, R_sq=1e12, R_L=5e3, vdd_up=(0.6, 1.8), dvdd=0.05, mismatch=False, down=False,
                    prog=prog)
    b = next(h for h in r["hist"] if h["n_fil"] > 0)
    assert any(abs(c - 12.0) < w / 2 + 0.5 for c, w in zip(b["centres_um"], b["widths_um"]))
