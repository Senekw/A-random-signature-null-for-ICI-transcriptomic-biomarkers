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
    purity = np.cos(0.6049872018 + 0.0001467884 * _rescale(combined))

    return pd.DataFrame({
        "immune_score": imm,
        "stromal_score": strm,
        "estimate_score": combined,
        "tumor_purity": purity,
        "n_immune_genes_covered": int(sum(g in expr.index for g in immune)),
        "n_stromal_genes_covered": int(sum(g in expr.index for g in stromal)),
    })


def _rescale(x: pd.Series) -> pd.Series:
    """Map an ssGSEA-scale combined score onto the published score range.

    The purity equation was fitted against ESTIMATE's own score scale.
    ssGSEA scores computed on a different platform and gene universe are
    on a different scale, so they are linearly mapped onto the published
    range before the transform. This preserves sample ranking exactly --
    the only property the analysis relies on -- while keeping purity in a
    plausible interval.
    """
    lo, hi = float(np.nanmin(x)), float(np.nanmax(x))
    if not np.isfinite(lo) or hi == lo:
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - lo) / (hi - lo) * 6000.0 - 3000.0
