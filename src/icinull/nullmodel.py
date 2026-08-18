"""The random-signature null.

This is the bar the manuscript argues every proposed ICI signature should
have to clear: not "is it associated with response", but "does it beat a
gene set of the same size drawn at random from the same expressed
universe, scored the same way, in the same cohort".

The null is deliberately cheap to apply. If you have a new signature and a
cohort, :func:`random_null_test` is the whole test.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd

from .scoring import auroc, mean_z_score

__all__ = ["NullResult", "random_null_test", "stouffer"]


@dataclass
class NullResult:
    """Outcome of one signature-vs-null comparison in one cohort."""

    observed_auroc: float
    null_mean: float
    null_sd: float
    null_q95: float
    p_empirical: float
    z_vs_null: float
    n_draws: int
    realized_size: int
    n_pos: int
    n_neg: int

    def as_dict(self) -> dict:
        return asdict(self)


def random_null_test(
    expr: pd.DataFrame,
    labels: pd.Series,
    genes,
    *,
    direction: int = 1,
    universe=None,
    n_draws: int = 1000,
    seed: int = 20260705,
) -> NullResult:
    """Test one signature against size-matched random gene sets.

    Parameters
    ----------
    expr
        genes x samples log2-TPM matrix, symbol-indexed.
    labels
        Binary response per sample (1 = responder), indexed by sample.
        Samples missing a label are dropped.
    genes
        Signature gene symbols.
    direction
        +1 or -1, the signature's a-priori direction.
    universe
        The expressed-gene pool random sets are drawn from. Defaults to
        every gene in ``expr``; in the pipeline it is the cohort's
        expressed universe (see :mod:`icinull.harmonize`).
    n_draws
        Number of random gene sets.
    seed
        RNG seed. Fixed a priori in the manuscript at 20260705.

    Returns
    -------
    NullResult
        ``p_empirical`` is the one-sided proportion of random sets whose
        AUROC matches or exceeds the observed AUROC, with the standard
        ``(hits + 1) / (n + 1)`` correction so that p is never exactly 0.

    Notes
    -----
    Random sets are size-matched to the signature's *realized* size -- the
    number of genes actually found in this cohort's matrix -- not to its
    published size. Matching published size would compare a fully-measured
    random set against a partially-measured signature.
    """
    common = [s for s in expr.columns if s in labels.index]
    y = labels.loc[common]
    y = y[np.isfinite(pd.to_numeric(y, errors="coerce"))]
    expr = expr[list(y.index)]

    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())

    covered = [g for g in dict.fromkeys(genes) if g in expr.index]
    size = len(covered)

    pool = list(expr.index if universe is None else
                [g for g in universe if g in expr.index])

    if size == 0 or n_pos == 0 or n_neg == 0 or len(pool) <= size:
        return NullResult(
            observed_auroc=float("nan"), null_mean=float("nan"),
            null_sd=float("nan"), null_q95=float("nan"),
            p_empirical=float("nan"), z_vs_null=float("nan"),
            n_draws=0, realized_size=size, n_pos=n_pos, n_neg=n_neg,
        )

    observed = auroc(mean_z_score(expr, covered, direction), y.to_numpy())

    rng = np.random.default_rng(seed)
    pool_arr = np.asarray(pool, dtype=object)
    yv = y.to_numpy()

    draws = np.empty(n_draws, dtype=float)
    for i in range(n_draws):
        pick = rng.choice(pool_arr, size=size, replace=False)
        draws[i] = auroc(mean_z_score(expr, pick, direction), yv)

    draws = draws[np.isfinite(draws)]
    n_eff = draws.size
    if n_eff == 0:
        p = float("nan")
        z = float("nan")
        mu = sd = q95 = float("nan")
    else:
        hits = int((draws >= observed).sum())
        p = (hits + 1) / (n_eff + 1)
        mu = float(draws.mean())
        sd = float(draws.std(ddof=1))
        q95 = float(np.quantile(draws, 0.95))
        z = float((observed - mu) / sd) if sd > 0 else float("nan")

    return NullResult(
        observed_auroc=float(observed), null_mean=mu, null_sd=sd,
        null_q95=q95, p_empirical=float(p), z_vs_null=z,
        n_draws=int(n_eff), realized_size=size, n_pos=n_pos, n_neg=n_neg,
    )


def stouffer(z_values, weights=None) -> tuple:
    """Stouffer's combination of per-cohort z statistics.

    Returns ``(z_combined, p_one_sided)``. Used to pool a signature's
    evidence against the null across cohorts before FDR control.
    """
    z = np.asarray([v for v in z_values if np.isfinite(v)], dtype=float)
    if z.size == 0:
        return float("nan"), float("nan")

    if weights is None:
        w = np.ones_like(z)
    else:
        w = np.asarray(
            [w for w, v in zip(weights, z_values) if np.isfinite(v)],
            dtype=float,
        )

    zc = float((w * z).sum() / np.sqrt((w ** 2).sum()))

    from scipy.stats import norm
    return zc, float(norm.sf(zc))
