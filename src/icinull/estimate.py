"""ESTIMATE immune / stromal scores and the derived tumor purity.

Yoshihara et al. (Nat Commun 2013) published two fixed 141-gene
signatures -- one immune, one stromal -- whose ssGSEA enrichment scores
estimate leukocyte and stromal content and, combined, tumor purity. The
manuscript uses them as an external, independently-derived readout of the
one axis the published ICI signatures turn out to share.

The gene lists live in ``config/estimate_gene_sets.csv`` (see that file's
provenance note); they are reference data, not results, so they are
committed rather than recomputed.

Purity follows the published transformation of the combined ESTIMATE
score. The absolute purity scale was calibrated on Affymetrix data, so
what matters for this analysis -- and all the manuscript uses -- is the
*ranking* of samples and the sign of its correlation with infiltration.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .scoring import rank_matrix, ssgsea_score

__all__ = ["estimate_gene_sets", "estimate_scores"]

# Yoshihara et al. (Nat Commun 2013) purity transform: purity = cos(A + B*s),
# monotone (decreasing) only while (A + B*s) lies in [0, pi].
_PURITY_A = 0.6049872018
_PURITY_B = 0.0001467884


@lru_cache(maxsize=1)
def estimate_gene_sets(path: str | None = None) -> tuple:
    """The 141-gene immune and stromal signatures as (immune, stromal)."""
    p = (Path(path) if path
         else Path(__file__).resolve().parents[2] / "config" / "estimate_gene_sets.csv")
    df = pd.read_csv(p, comment="#")
    immune = tuple(str(g) for g in df["immune_signature"].dropna())
    stromal = tuple(str(g) for g in df["stromal_signature"].dropna())
    return immune, stromal


def estimate_scores(expr: pd.DataFrame, alpha: float = 0.25) -> pd.DataFrame:
    """Per-sample immune, stromal, combined ESTIMATE score and purity.

    Parameters
    ----------
    expr
        genes x samples log2-TPM matrix, symbol-indexed. Pass the
        cohort's expressed universe, so the ssGSEA ranking is over the
        same gene space the rest of the analysis uses.

    Returns
    -------
    DataFrame indexed by sample with columns ``immune_score``,
    ``stromal_score``, ``estimate_score``, ``tumor_purity`` and the gene
    coverage of each signature.
    """
    immune, stromal = estimate_gene_sets()
    ranks = rank_matrix(expr)

    imm = ssgsea_score(expr, immune, 1, alpha=alpha, ranks=ranks)
    strm = ssgsea_score(expr, stromal, 1, alpha=alpha, ranks=ranks)
    combined = imm + strm

    # Yoshihara et al. eq. relating the combined score to purity.
    purity = np.cos(_PURITY_A + _PURITY_B * _rescale(combined))

    return pd.DataFrame({
        "immune_score": imm,
        "stromal_score": strm,
        "estimate_score": combined,
        "tumor_purity": purity,
        "n_immune_genes_covered": int(sum(g in expr.index for g in immune)),
        "n_stromal_genes_covered": int(sum(g in expr.index for g in stromal)),
    })


def _rescale(x: pd.Series) -> pd.Series:
    """Map a combined ESTIMATE score onto the purity transform's domain.

    Yoshihara et al.'s purity equation, ``cos(A + B*s)``, was fitted against
    ESTIMATE's own score scale on Affymetrix arrays. ssGSEA scores computed
    here, on RNA-seq over a different gene universe, are on a different
    scale, so they must be mapped onto that domain first.

    The mapping is by **within-cohort rank**, not by score value. Two
    reasons:

    1. ``cos`` is monotone only while its argument lies in ``[0, pi]``. A
       value-based map has to be repaired to guarantee that, and a single
       outlying sample can otherwise push the argument past ``pi``, where
       the cosine turns back upward and the most-infiltrated samples are
       reported as the *purest* -- silently inverting the relationship the
       confound analysis exists to measure.
    2. Cohort sizes here span 22 to 348 samples, so score extremes are not
       comparable across cohorts. A rank map is unaffected by how extreme
       an outlier is: an outlier is simply the last rank.

    A rank map is strictly monotone and tie-free for distinct scores, so the
    sample *ranking* -- the only property the manuscript's analysis uses --
    is preserved exactly, and equal scores stay equal (average ranks).

    Consequently the absolute values are **not** calibrated tumour-purity
    estimates and should not be reported as percentages; they are a monotone
    transform of the combined score, used for ranking and for the sign of
    its correlation with the signature axis.
    """
    v = pd.to_numeric(x, errors="coerce")
    lo, hi = np.nanpercentile(v, 5), np.nanpercentile(v, 95)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        # degenerate cohort (all-equal or unusably small): fall back to the
        # full range, and to zeros only if that is degenerate too
        lo, hi = float(np.nanmin(v)), float(np.nanmax(v))
        if not np.isfinite(lo) or hi == lo:
            return pd.Series(np.zeros(len(v)), index=v.index)
    # The purity transform is cos(A + B*s), monotone only while its argument
    # stays in [0, pi]. Past pi the cosine turns back up and the
    # highest-scoring samples get reported as HIGH purity -- inverting the
    # very relationship being measured.
    #
    # Rather than map the score range onto the argument range and then repair
    # it, map RANKS directly: the sample ranking is the only property the
    # analysis uses, and a rank map is by construction strictly monotone,
    # tie-free for distinct scores, bounded inside the monotone branch, and
    # completely insensitive to how extreme an outlier is -- an outlier is
    # just the last rank. `lo`/`hi` above are retained only to document the
    # published score scale this replaces.
    del lo, hi

    n = int(v.notna().sum())
    if n < 2:
        return pd.Series(np.zeros(len(v)), index=v.index)

    # ranks in (0, 1), average-ranked so equal scores stay equal
    u = v.rank(method="average", na_option="keep") / (n + 1.0)

    margin = 0.02
    arg = margin + u * (np.pi - 2.0 * margin)
    return (arg - _PURITY_A) / _PURITY_B
