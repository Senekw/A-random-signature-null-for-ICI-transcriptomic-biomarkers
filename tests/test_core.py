"""Tests for the scoring, null-model and metric primitives.

These run without any downloaded cohort: they construct small synthetic
matrices with known answers. The point is that the pieces the paper's
claims rest on -- the AUROC, the null's calibration under the null
hypothesis, the direction convention, size matching -- are checked
against analytic expectations rather than against a previous run.

Run:  pytest tests/ -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from icinull import (  # noqa: E402
    auroc,
    auroc_se_hanley_mcneil,
    expressed_universe,
    logit,
    mean_z_score,
    random_null_test,
    ssgsea_score,
    stouffer,
)
from icinull.scoring import expit  # noqa: E402


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------

def make_expr(n_genes=400, n_samples=60, seed=0, signal_genes=(), effect=0.0):
    """Random log2-TPM-like matrix, optionally with planted signal.

    The first half of the samples are 'responders'; ``signal_genes`` are
    shifted up by ``effect`` in those samples.
    """
    rng = np.random.default_rng(seed)
    genes = [f"G{i:04d}" for i in range(n_genes)]
    samples = [f"S{i:03d}" for i in range(n_samples)]
    X = rng.normal(5.0, 1.0, size=(n_genes, n_samples))
    expr = pd.DataFrame(X, index=genes, columns=samples)

    y = pd.Series([1] * (n_samples // 2) + [0] * (n_samples - n_samples // 2),
                  index=samples)
    if effect and len(signal_genes):
        expr.loc[list(signal_genes), y[y == 1].index] += effect
    return expr, y


# --------------------------------------------------------------------------
# AUROC
# --------------------------------------------------------------------------

def test_auroc_perfect_separation():
    assert auroc([4, 3, 2, 1], [1, 1, 0, 0]) == pytest.approx(1.0)


def test_auroc_perfect_inversion():
    assert auroc([1, 2, 3, 4], [1, 1, 0, 0]) == pytest.approx(0.0)


def test_auroc_all_ties_is_chance():
    assert auroc([1, 1, 1, 1], [1, 1, 0, 0]) == pytest.approx(0.5)


def test_auroc_matches_sklearn():
    rng = np.random.default_rng(1)
    s = rng.normal(size=200)
    y = (rng.random(200) < 0.4).astype(int)
    from sklearn.metrics import roc_auc_score
    assert auroc(s, y) == pytest.approx(roc_auc_score(y, s), abs=1e-12)


def test_auroc_single_class_is_nan():
    assert np.isnan(auroc([1, 2, 3], [1, 1, 1]))


def test_auroc_drops_missing():
    a = auroc([1, 2, np.nan, 4], [0, 0, 1, 1])
    assert a == pytest.approx(auroc([1, 2, 4], [0, 0, 1]))


def test_hanley_mcneil_se_shrinks_with_n():
    small = auroc_se_hanley_mcneil(0.7, 10, 10)
    large = auroc_se_hanley_mcneil(0.7, 100, 100)
    assert 0 < large < small


def test_logit_roundtrip():
    for p in (0.05, 0.5, 0.63, 0.95):
        assert expit(logit(p)) == pytest.approx(p, abs=1e-9)


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------

def test_mean_z_is_centred():
    expr, _ = make_expr()
    s = mean_z_score(expr, expr.index[:20], 1)
    assert s.mean() == pytest.approx(0.0, abs=1e-9)


def test_mean_z_direction_flips_sign():
    expr, _ = make_expr()
    genes = expr.index[:15]
    up = mean_z_score(expr, genes, 1)
    dn = mean_z_score(expr, genes, -1)
    assert np.allclose(up.to_numpy(), -dn.to_numpy())


def test_mean_z_direction_flips_auroc_about_half():
    genes = [f"G{i:04d}" for i in range(10)]
    expr, y = make_expr(signal_genes=genes, effect=1.5, seed=7)
    a_up = auroc(mean_z_score(expr, genes, 1), y.to_numpy())
    a_dn = auroc(mean_z_score(expr, genes, -1), y.to_numpy())
    assert a_up > 0.9
    assert a_up + a_dn == pytest.approx(1.0, abs=1e-9)


def test_mean_z_ignores_absent_genes():
    expr, _ = make_expr()
    genes = list(expr.index[:10])
    s1 = mean_z_score(expr, genes, 1)
    s2 = mean_z_score(expr, genes + ["NOT_A_GENE", "ALSO_MISSING"], 1)
    assert np.allclose(s1.to_numpy(), s2.to_numpy())


def test_mean_z_all_genes_absent_is_nan():
    expr, _ = make_expr()
    s = mean_z_score(expr, ["NOPE1", "NOPE2"], 1)
    assert s.isna().all()


def test_mean_z_drops_zero_variance_genes():
    expr, _ = make_expr()
    expr.loc["G0000"] = 7.0                     # constant across samples
    s = mean_z_score(expr, ["G0000", "G0001"], 1)
    assert s.notna().all()
    assert np.allclose(s.to_numpy(), mean_z_score(expr, ["G0001"], 1).to_numpy())


def test_ssgsea_recovers_planted_signal():
    genes = [f"G{i:04d}" for i in range(12)]
    expr, y = make_expr(signal_genes=genes, effect=2.0, seed=3)
    a = auroc(ssgsea_score(expr, genes, 1), y.to_numpy())
    assert a > 0.85


def test_ssgsea_agrees_with_mean_z_in_direction():
    genes = [f"G{i:04d}" for i in range(12)]
    expr, y = make_expr(signal_genes=genes, effect=1.5, seed=11)
    a_mz = auroc(mean_z_score(expr, genes, 1), y.to_numpy())
    a_ss = auroc(ssgsea_score(expr, genes, 1), y.to_numpy())
    assert a_mz > 0.5 and a_ss > 0.5


# --------------------------------------------------------------------------
# the random-signature null -- the paper's central instrument
# --------------------------------------------------------------------------

def test_null_is_calibrated_under_the_null():
    """A random gene set should not beat the random null.

    With no planted signal, the empirical p-values of arbitrary gene sets
    must be roughly uniform, so the false-positive rate at p<0.05 stays
    near 5%. This is the property that makes "30.5% of tests beat null"
    interpretable at all.
    """
    expr, y = make_expr(n_genes=600, n_samples=80, seed=5)
    rng = np.random.default_rng(99)
    ps = []
    for i in range(60):
        genes = rng.choice(expr.index.to_numpy(), size=25, replace=False)
        r = random_null_test(expr, y, genes, n_draws=200, seed=1000 + i)
        ps.append(r.p_empirical)

    ps = np.array(ps)
    assert np.isfinite(ps).all()
    assert 0.15 < ps.mean() < 0.85           # roughly uniform, not skewed
    assert (ps < 0.05).mean() < 0.20         # FPR controlled at small n


def test_null_detects_real_signal():
    genes = [f"G{i:04d}" for i in range(20)]
    expr, y = make_expr(signal_genes=genes, effect=2.0, seed=13)
    r = random_null_test(expr, y, genes, n_draws=300, seed=42)
    assert r.observed_auroc > 0.85
    assert r.p_empirical < 0.01
    assert r.observed_auroc > r.null_q95
    assert r.z_vs_null > 2


def test_null_size_matches_realized_not_published_size():
    expr, y = make_expr()
    genes = list(expr.index[:10]) + ["MISSING1", "MISSING2", "MISSING3"]
    r = random_null_test(expr, y, genes, n_draws=50, seed=1)
    assert r.realized_size == 10


def test_null_is_reproducible_for_a_fixed_seed():
    expr, y = make_expr(seed=2)
    genes = expr.index[:15]
    a = random_null_test(expr, y, genes, n_draws=100, seed=20260705)
    b = random_null_test(expr, y, genes, n_draws=100, seed=20260705)
    assert a.p_empirical == b.p_empirical
    assert a.null_mean == pytest.approx(b.null_mean)


def test_null_differs_across_seeds():
    expr, y = make_expr(seed=2)
    genes = expr.index[:15]
    a = random_null_test(expr, y, genes, n_draws=100, seed=1)
    b = random_null_test(expr, y, genes, n_draws=100, seed=2)
    assert a.null_mean != b.null_mean


def test_null_p_never_zero():
    """(hits + 1) / (n + 1) keeps p bounded away from zero."""
    genes = [f"G{i:04d}" for i in range(20)]
    expr, y = make_expr(signal_genes=genes, effect=6.0, seed=17)
    r = random_null_test(expr, y, genes, n_draws=100, seed=1)
    assert r.p_empirical > 0
    assert r.p_empirical == pytest.approx(1 / 101, abs=1e-9)


def test_null_respects_a_restricted_universe():
    expr, y = make_expr()
    universe = list(expr.index[:50])
    r = random_null_test(expr, y, expr.index[:10], n_draws=30,
                         universe=universe, seed=1)
    assert r.n_draws == 30
    assert np.isfinite(r.null_mean)


def test_null_returns_nan_when_one_class_is_empty():
    expr, y = make_expr()
    y_one = pd.Series(1, index=y.index)
    r = random_null_test(expr, y_one, expr.index[:10], n_draws=10)
    assert np.isnan(r.observed_auroc)
    assert r.n_draws == 0


# --------------------------------------------------------------------------
# pooling
# --------------------------------------------------------------------------

def test_stouffer_combines_concordant_evidence():
    z, p = stouffer([2.0, 2.0, 2.0, 2.0])
    assert z == pytest.approx(4.0)
    assert p < 0.001


def test_stouffer_cancels_opposing_evidence():
    z, p = stouffer([2.0, -2.0])
    assert z == pytest.approx(0.0)
    assert p == pytest.approx(0.5)


def test_stouffer_ignores_non_finite():
    z1, _ = stouffer([2.0, np.nan, 2.0])
    z2, _ = stouffer([2.0, 2.0])
    assert z1 == pytest.approx(z2)


def test_stouffer_empty_is_nan():
    z, p = stouffer([np.nan, np.nan])
    assert np.isnan(z) and np.isnan(p)


# --------------------------------------------------------------------------
# expressed universe
# --------------------------------------------------------------------------

def test_expressed_universe_excludes_floor_only_genes():
    expr, _ = make_expr(n_genes=100, n_samples=20)
    floor = float(expr.to_numpy().min())
    expr.loc["G0000"] = floor                      # never detected
    expr.loc["G0001"] = [floor] * 18 + [8.0, 8.0]  # detected in 10%
    universe = expressed_universe(expr, 0.20)
    assert "G0000" not in universe
    assert "G0001" not in universe
    assert "G0050" in universe


# --------------------------------------------------------------------------
# the null's vectorized fast path must equal the naive reference
# --------------------------------------------------------------------------

def _reference_null_draws(expr, y, genes, universe, n_draws, seed,
                          direction=1):
    """The obvious, slow implementation: re-standardize per draw."""
    rng = np.random.default_rng(seed)
    pool = np.asarray([g for g in universe if g in expr.index], dtype=object)
    covered = [g for g in dict.fromkeys(genes) if g in expr.index]
    out = np.empty(n_draws, dtype=float)
    for i in range(n_draws):
        pick = rng.choice(pool, size=len(covered), replace=False)
        out[i] = auroc(mean_z_score(expr, pick, direction), y.to_numpy())
    return out


def test_null_fast_path_matches_reference():
    """random_null_test precomputes the z-matrix once instead of per draw.

    That is only a legitimate optimization if it changes nothing, so this
    pins it against the naive loop: same seed, same draws, bit-for-bit.
    """
    expr, y = make_expr(n_genes=600, n_samples=50, seed=4)
    genes = list(expr.index[:15])
    universe = list(expr.index)

    got = random_null_test(expr, y, genes, universe=universe,
                           n_draws=200, seed=99)
    ref = _reference_null_draws(expr, y, genes, universe, 200, 99)

    assert got.null_mean == pytest.approx(float(ref.mean()), abs=1e-12)
    assert got.null_sd == pytest.approx(float(ref.std(ddof=1)), abs=1e-12)
    assert got.null_q95 == pytest.approx(float(np.quantile(ref, 0.95)),
                                         abs=1e-12)


def test_null_fast_path_handles_zero_variance_genes():
    """A gene with no variance has an undefined z-score; it must not make
    the whole draw NaN."""
    expr, y = make_expr(n_genes=300, n_samples=40, seed=5)
    expr.iloc[0] = 7.0                      # constant across samples
    r = random_null_test(expr, y, list(expr.index[:10]),
                         universe=list(expr.index), n_draws=50, seed=1)
    assert np.isfinite(r.null_mean)
    assert r.n_draws == 50


def test_null_direction_flips_the_observed_auroc():
    expr, y = make_expr(n_genes=400, n_samples=60, seed=6,
                        signal_genes=[f"G{i:04d}" for i in range(10)],
                        effect=1.5)
    genes = [f"G{i:04d}" for i in range(10)]
    up = random_null_test(expr, y, genes, direction=1, n_draws=30, seed=2)
    dn = random_null_test(expr, y, genes, direction=-1, n_draws=30, seed=2)
    assert up.observed_auroc == pytest.approx(1.0 - dn.observed_auroc, abs=1e-12)


# --------------------------------------------------------------------------
# ssGSEA against the published definition
# --------------------------------------------------------------------------

def _ssgsea_from_the_paper(expr, genes, alpha=0.25):
    """Barbie et al. (Nature 2009) ssGSEA, transcribed from the definition.

    Independent of the optimized implementation in icinull.scoring: rank
    genes within each sample, walk the descending list accumulating
    |rank|**alpha for in-set genes (normalized to sum 1) against a uniform
    step for out-of-set genes, and integrate the difference of the two CDFs.
    """
    X = expr.to_numpy(dtype=float)
    in_set = np.array([g in set(genes) for g in expr.index])
    n = X.shape[0]
    out = []
    for j in range(X.shape[1]):
        rank = pd.Series(X[:, j]).rank(method="average").to_numpy()
        order = np.argsort(-rank, kind="stable")
        w = np.abs(rank[order]) ** alpha
        hit = in_set[order]
        hw = np.where(hit, w, 0.0)
        if hw.sum() == 0:
            out.append(np.nan)
            continue
        cdf_in = np.cumsum(hw) / hw.sum()
        cdf_out = np.cumsum(~hit) / (n - hit.sum())
        out.append(float(np.sum(cdf_in - cdf_out) / n))
    return pd.Series(out, index=expr.columns)


@pytest.mark.parametrize("alpha", [0.25, 0.5, 1.0])
def test_ssgsea_matches_the_published_definition(alpha):
    expr, _ = make_expr(n_genes=500, n_samples=25, seed=8)
    genes = list(expr.index[:20])
    got = ssgsea_score(expr, genes, 1, alpha=alpha)
    ref = _ssgsea_from_the_paper(expr, genes, alpha=alpha)
    assert np.abs(got.to_numpy() - ref.to_numpy()).max() < 1e-12


def test_ssgsea_alpha_is_actually_used():
    """A different exponent must give a different answer -- guards against
    the weight silently defaulting."""
    expr, _ = make_expr(n_genes=400, n_samples=20, seed=9)
    genes = list(expr.index[:15])
    a = ssgsea_score(expr, genes, 1, alpha=0.25)
    b = ssgsea_score(expr, genes, 1, alpha=1.0)
    assert not np.allclose(a.to_numpy(), b.to_numpy())


def test_ssgsea_precomputed_ranks_match_internal_ranking():
    expr, _ = make_expr(n_genes=300, n_samples=15, seed=10)
    genes = list(expr.index[:12])
    from icinull.scoring import rank_matrix
    a = ssgsea_score(expr, genes, 1)
    b = ssgsea_score(expr, genes, 1, ranks=rank_matrix(expr))
    assert np.abs(a.to_numpy() - b.to_numpy()).max() < 1e-12
