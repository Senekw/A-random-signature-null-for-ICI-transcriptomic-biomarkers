"""Signature scoring and discrimination metrics.

Two scoring methods, matching the manuscript:

* ``mean_z``  -- mean of per-gene z-scores, standardized across samples
  within a cohort. Deterministic, and the primary method.
* ``ssgsea``  -- single-sample GSEA (Barbie et al. 2009) with alpha=0.25,
  used as a sensitivity analysis.

Both return one score per sample, with the signature's a-priori direction
applied as a sign, so that a higher score always means "predicted
responder" regardless of whether the gene set is a response or a
resistance set.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = [
    "mean_z_score",
    "ssgsea_score",
    "score_signature",
    "auroc",
    "auroc_se_hanley_mcneil",
    "logit",
    "expit",
    "logit_se",
]


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------

def _covered(expr: pd.DataFrame, genes) -> list:
    """Signature genes present in the expression matrix, order preserved."""
    present = set(expr.index)
    seen, out = set(), []
    for g in genes:
        if g in present and g not in seen:
            seen.add(g)
            out.append(g)
    return out


def mean_z_score(expr: pd.DataFrame, genes, direction: int = 1) -> pd.Series:
    """Mean per-gene z-score across the covered signature genes.

    Parameters
    ----------
    expr
        genes x samples expression matrix (log2-TPM), symbol-indexed.
    genes
        Signature gene symbols.
    direction
        +1 or -1; multiplied into the returned score.

    Notes
    -----
    Genes are standardized across samples *within this matrix*, so the
    score is cohort-relative by construction. Genes with zero variance
    contribute nothing (they are dropped rather than producing NaN).
    """
    sub = expr.loc[_covered(expr, genes)]
    if sub.empty:
        return pd.Series(np.nan, index=expr.columns, dtype=float)

    sd = sub.std(axis=1, ddof=1)
    sub = sub.loc[sd > 0]
    if sub.empty:
        return pd.Series(np.nan, index=expr.columns, dtype=float)

    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1), axis=0)
    return direction * z.mean(axis=0)


def ssgsea_score(
    expr: pd.DataFrame,
    genes,
    direction: int = 1,
    alpha: float = 0.25,
    *,
    ranks: pd.DataFrame | None = None,
) -> pd.Series:
    """Single-sample GSEA enrichment score (Barbie et al. 2009).

    The weighted running-sum statistic over genes ranked within each
    sample, using rank**alpha weights, normalized by the number of genes.

    Parameters
    ----------
    ranks
        Optional precomputed within-sample ranks (see :func:`rank_matrix`).
        Ranking is the expensive step, so callers scoring many signatures
        against the same cohort should compute it once and pass it in.
    """
    hits = _covered(expr, genes)
    if not hits:
        return pd.Series(np.nan, index=expr.columns, dtype=float)

    R = rank_matrix(expr) if ranks is None else ranks
    n_genes = R.shape[0]
    idx = R.index.get_indexer(hits)
    idx = idx[idx >= 0]

    scores = np.empty(R.shape[1], dtype=float)
    rank_values = R.to_numpy()

    for j in range(rank_values.shape[1]):
        col = rank_values[:, j]
        order = np.argsort(-col, kind="stable")          # descending
        in_set = np.zeros(n_genes, dtype=bool)
        in_set[idx] = True
        in_set_sorted = in_set[order]
        weights = np.abs(col[order]) ** alpha

        hit_w = np.where(in_set_sorted, weights, 0.0)
        total_hit = hit_w.sum()
        if total_hit == 0:
            scores[j] = np.nan
            continue

        cdf_hit = np.cumsum(hit_w) / total_hit
        n_miss = n_genes - in_set_sorted.sum()
        cdf_miss = np.cumsum(~in_set_sorted) / n_miss
        scores[j] = np.sum(cdf_hit - cdf_miss) / n_genes

    return direction * pd.Series(scores, index=R.columns)


def rank_matrix(expr: pd.DataFrame) -> pd.DataFrame:
    """Within-sample ranks of every gene (ties averaged)."""
    return expr.rank(axis=0, method="average")


def score_signature(
    expr: pd.DataFrame,
    genes,
    direction: int = 1,
    method: str = "mean_z",
    **kwargs,
) -> pd.Series:
    """Dispatch to the requested scoring method."""
    if method == "mean_z":
        return mean_z_score(expr, genes, direction)
    if method == "ssgsea":
        return ssgsea_score(expr, genes, direction, **kwargs)
    raise ValueError(f"unknown scoring method: {method!r}")


# --------------------------------------------------------------------------
# Discrimination
# --------------------------------------------------------------------------

def auroc(scores, labels) -> float:
    """AUROC via the Mann-Whitney U identity, ties at 0.5.

    ``labels`` is 1 for responders, 0 for non-responders. Samples with a
    missing score or missing label are dropped. Returns NaN when either
    class is empty.
    """
    s = np.asarray(scores, dtype=float)
    y = np.asarray(labels, dtype=float)
    ok = np.isfinite(s) & np.isfinite(y)
    s, y = s[ok], y[ok]

    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")

    r = pd.Series(s).rank(method="average").to_numpy()
    u = r[y == 1].sum() - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def auroc_se_hanley_mcneil(a: float, n_pos: int, n_neg: int) -> float:
    """Hanley-McNeil standard error of an AUROC."""
    if not np.isfinite(a) or n_pos < 1 or n_neg < 1:
        return float("nan")
    q1 = a / (2.0 - a)
    q2 = 2.0 * a * a / (1.0 + a)
    var = (
        a * (1.0 - a)
        + (n_pos - 1) * (q1 - a * a)
        + (n_neg - 1) * (q2 - a * a)
    ) / (n_pos * n_neg)
    return float(np.sqrt(var)) if var > 0 else float("nan")


def logit(p: float, eps: float = 1e-6) -> float:
    p = min(max(float(p), eps), 1.0 - eps)
    return float(np.log(p / (1.0 - p)))


def expit(x: float) -> float:
    return float(1.0 / (1.0 + np.exp(-float(x))))


def logit_se(a: float, se: float, eps: float = 1e-6) -> float:
    """Delta-method transform of an AUROC standard error to logit scale."""
    a = min(max(float(a), eps), 1.0 - eps)
    return float(se / (a * (1.0 - a)))
