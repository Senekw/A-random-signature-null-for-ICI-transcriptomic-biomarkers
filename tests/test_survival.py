"""The Cox fallback must agree with lifelines to numerical tolerance.

The Methods specify lifelines for the overall-survival log-hazard-ratios.
``icinull.survival`` uses it when installed and otherwise falls back to a
Newton-Raphson fit of the Breslow partial likelihood, so that the survival
step is not a hard dependency. That fallback is only legitimate if it gives
the same answer, which is what these tests check -- on synthetic data with
a known effect, with ties, and with heavy censoring.

Run:  pytest tests/test_survival.py -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from icinull.survival import cox_loghr_per_sd  # noqa: E402

lifelines = pytest.importorskip("lifelines")


def _lifelines_fit(x, t, e):
    from lifelines import CoxPHFitter
    x = np.asarray(x, dtype=float)
    z = (x - x.mean()) / x.std(ddof=1)
    df = pd.DataFrame({"x": z, "T": np.asarray(t, float),
                       "E": np.asarray(e, int)})
    cph = CoxPHFitter().fit(df, duration_col="T", event_col="E")
    return float(cph.params_["x"]), float(cph.standard_errors_["x"])


def _synthetic(n=120, beta=-0.5, seed=0, tie_grid=None, censor=0.7):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=n)
    t = rng.exponential(np.exp(beta * x))
    if tie_grid:                            # coarsen times to force ties
        t = np.round(t / tie_grid) * tie_grid + tie_grid
    e = (rng.random(n) < censor).astype(int)
    return pd.Series(x), pd.Series(t), pd.Series(e)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_matches_lifelines_beta_and_se(seed):
    x, t, e = _synthetic(seed=seed)
    b_own, se_own, _ = cox_loghr_per_sd(x, t, e)
    b_ref, se_ref = _lifelines_fit(x, t, e)
    assert b_own == pytest.approx(b_ref, abs=1e-4)
    assert se_own == pytest.approx(se_ref, abs=1e-4)


def test_matches_lifelines_with_heavy_ties():
    """Breslow's tie handling is what the fallback implements."""
    x, t, e = _synthetic(n=150, seed=7, tie_grid=0.25)
    assert t.duplicated().sum() > 20, "fixture failed to create ties"
    b_own, se_own, _ = cox_loghr_per_sd(x, t, e)
    b_ref, se_ref = _lifelines_fit(x, t, e)
    assert b_own == pytest.approx(b_ref, abs=5e-3)
    assert se_own == pytest.approx(se_ref, abs=5e-3)


def test_matches_lifelines_under_heavy_censoring():
    x, t, e = _synthetic(n=200, seed=11, censor=0.25)
    b_own, se_own, _ = cox_loghr_per_sd(x, t, e)
    b_ref, se_ref = _lifelines_fit(x, t, e)
    assert b_own == pytest.approx(b_ref, abs=1e-3)
    assert se_own == pytest.approx(se_ref, abs=1e-3)


def test_recovers_the_planted_direction():
    x, t, e = _synthetic(n=300, beta=-0.8, seed=3)
    b, se, p = cox_loghr_per_sd(x, t, e)
    assert b > 0            # higher score -> longer survival -> positive beta
    assert p < 0.01


def test_degenerate_inputs_are_nan_not_an_exception():
    n = 40
    x = pd.Series(np.ones(n))                     # zero variance
    t = pd.Series(np.linspace(1, 5, n))
    e = pd.Series(np.ones(n, dtype=int))
    b, se, p = cox_loghr_per_sd(x, t, e)
    assert np.isnan(b)

    x2 = pd.Series(np.random.default_rng(0).normal(size=n))
    e0 = pd.Series(np.zeros(n, dtype=int))         # no events at all
    b2, _, _ = cox_loghr_per_sd(x2, t, e0)
    assert np.isnan(b2)
