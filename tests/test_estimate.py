"""ESTIMATE purity must be strictly monotone in the combined score.

The manuscript uses purity only through its *ranking* and the sign of its
correlation with the signature axis. That makes monotonicity the property
worth testing, and it is easy to lose: purity is ``cos(A + B*s)``, which
turns back upward once its argument passes pi, so an outlying sample can
silently invert the relationship for the highest-scoring samples -- exactly
the samples the confound analysis is about.

Run:  pytest tests/test_estimate.py -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from icinull.estimate import _PURITY_A, _PURITY_B, _rescale  # noqa: E402


def _purity(scores: pd.Series) -> pd.Series:
    return pd.Series(np.cos(_PURITY_A + _PURITY_B * _rescale(scores)),
                     index=scores.index)


CASES = {
    "clean_normal": pd.Series(np.random.default_rng(0).normal(size=200)),
    "single_high_outlier": None,     # filled below
    "single_low_outlier": None,
    "tiny_cohort": pd.Series(np.random.default_rng(2).normal(size=22)),
    "large_cohort": pd.Series(np.random.default_rng(3).normal(size=348)),
    "heavy_tailed": pd.Series(np.random.default_rng(4).standard_t(2, size=150)),
}
_o = pd.Series(np.random.default_rng(1).normal(size=200))
_o.iloc[0] = 50.0
CASES["single_high_outlier"] = _o
_u = pd.Series(np.random.default_rng(1).normal(size=200))
_u.iloc[0] = -50.0
CASES["single_low_outlier"] = _u


@pytest.mark.parametrize("name", sorted(CASES))
def test_purity_argument_stays_on_the_monotone_branch(name):
    arg = _PURITY_A + _PURITY_B * _rescale(CASES[name])
    assert arg.min() >= 0.0, f"{name}: argument below 0, cosine not monotone"
    assert arg.max() <= np.pi, f"{name}: argument past pi, purity wraps around"


@pytest.mark.parametrize("name", sorted(CASES))
def test_purity_is_strictly_decreasing_in_score(name):
    """Higher infiltration score must mean lower purity, with no ties."""
    s = CASES[name]
    rho = spearmanr(s, _purity(s)).statistic
    assert rho == pytest.approx(-1.0, abs=1e-9), f"{name}: rho={rho}"


@pytest.mark.parametrize("name", sorted(CASES))
def test_rescale_preserves_ranking_exactly(name):
    s = CASES[name]
    assert spearmanr(s, _rescale(s)).statistic == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("name", sorted(CASES))
def test_no_two_distinct_scores_collapse_to_one_value(name):
    """A clip would tie the tails together; an affine map must not."""
    s = CASES[name]
    assert _rescale(s).nunique() == s.nunique()


def test_an_outlier_does_not_shift_the_bulk():
    """Percentile anchors, not min/max: one extreme sample must not move
    everyone else's purity appreciably."""
    base = pd.Series(np.random.default_rng(5).normal(size=200))
    with_out = base.copy()
    with_out.iloc[0] = 40.0

    core = base.index[1:]
    shift = (_purity(with_out).loc[core] - _purity(base).loc[core]).abs().max()
    assert shift < 0.10, f"bulk purity moved by {shift:.3f} due to one outlier"


def test_degenerate_input_returns_finite_values():
    flat = pd.Series(np.full(30, 2.5))
    out = _rescale(flat)
    assert np.isfinite(out).all()
    assert np.isfinite(_purity(flat)).all()


def test_purity_stays_in_the_valid_cosine_range():
    for s in CASES.values():
        p = _purity(s)
        assert p.min() >= -1.0 and p.max() <= 1.0
