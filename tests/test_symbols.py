"""Symbol harmonization must be conservative and vintage-invariant.

The property that matters is the one the pipeline depends on: a gene set
written with *retired* symbols and the same set written with *current*
symbols must score identically, in a cohort of either annotation vintage.
Otherwise a signature's realized size -- and with it the size-matched null
and its meta-analysis weight -- depends on which GENCODE build a cohort was
processed against.

Run:  pytest tests/test_symbols.py -q
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from icinull.symbols import (  # noqa: E402
    harmonize_expression,
    harmonize_gene_list,
    harmonize_index,
)

# A small hand-built vocabulary, so these tests need no network and no cache.
OFFICIAL = {"IDO1", "CXCR5", "H2AC11", "EZR", "KMT2A", "CD8A", "GZMB", "STAT1"}
ALIAS = {"INDO": "IDO1", "BLR1": "CXCR5", "HIST1H2AG": "H2AC11",
         "VIL2": "EZR", "MLL": "KMT2A"}


def test_current_symbols_are_never_rewritten():
    got = harmonize_index(["IDO1", "CD8A", "STAT1"], OFFICIAL, ALIAS)
    assert got == ["IDO1", "CD8A", "STAT1"]


def test_retired_symbols_map_to_current():
    got = harmonize_index(["INDO", "BLR1", "HIST1H2AG"], OFFICIAL, ALIAS)
    assert got == ["IDO1", "CXCR5", "H2AC11"]


def test_unknown_symbols_are_left_alone_not_dropped():
    """An unmappable identifier must survive, so it can be reported."""
    got = harmonize_index(["AF107846", "GCLC/GCLM", "IDO1"], OFFICIAL, ALIAS)
    assert got == ["AF107846", "GCLC/GCLM", "IDO1"]


def test_a_symbol_that_is_official_wins_over_being_an_alias():
    """If X is both a current symbol and a synonym of Y, X stays X."""
    official = OFFICIAL | {"RAB7A"}
    alias = dict(ALIAS)
    alias["CD8A"] = "GZMB"          # pathological: CD8A is official
    assert harmonize_index(["CD8A"], official, alias) == ["CD8A"]


def test_gene_list_reports_what_changed():
    genes, renamed = harmonize_gene_list(["INDO", "CD8A", "VIL2"],
                                         OFFICIAL, ALIAS)
    assert genes == ["IDO1", "CD8A", "EZR"]
    assert renamed == {"INDO": "IDO1", "VIL2": "EZR"}


def test_gene_list_deduplicates_after_mapping():
    """INDO and IDO1 in one set are one gene, not two."""
    genes, _ = harmonize_gene_list(["INDO", "IDO1", "CD8A"], OFFICIAL, ALIAS)
    assert genes == ["IDO1", "CD8A"]


# --------------------------------------------------------------------------
# expression matrices
# --------------------------------------------------------------------------

def _expr(index, n_samples=6, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(rng.normal(5, 1, size=(len(index), n_samples)),
                        index=index,
                        columns=[f"S{i}" for i in range(n_samples)])


def test_expression_index_is_rewritten_and_reported():
    e = _expr(["INDO", "CD8A", "HIST1H2AG"])
    out, rep = harmonize_expression(e, OFFICIAL, ALIAS)
    assert set(out.index) == {"IDO1", "CD8A", "H2AC11"}
    assert rep["n_renamed"] == 2
    assert rep["n_collapsed"] == 0


def test_collision_collapses_to_highest_mean_expression():
    e = _expr(["INDO", "IDO1"], n_samples=4)
    e.loc["INDO"] = [1.0, 1.0, 1.0, 1.0]
    e.loc["IDO1"] = [9.0, 9.0, 9.0, 9.0]     # higher mean -> kept
    out, rep = harmonize_expression(e, OFFICIAL, ALIAS)
    assert list(out.index) == ["IDO1"]
    assert rep["n_collapsed"] == 1
    assert out.loc["IDO1"].mean() == pytest.approx(9.0)


def test_vintage_invariance_is_the_whole_point():
    """The same biology under two annotation vintages must score the same.

    ``old`` is a cohort processed with retired symbols, ``new`` the same
    values under current symbols. A signature written either way must see
    the same genes in both after harmonization.
    """
    values = _expr(["A", "B", "C"], n_samples=5).to_numpy()
    old = pd.DataFrame(values, index=["INDO", "BLR1", "HIST1H2AG"],
                       columns=[f"S{i}" for i in range(5)])
    new = pd.DataFrame(values, index=["IDO1", "CXCR5", "H2AC11"],
                       columns=[f"S{i}" for i in range(5)])

    ho, _ = harmonize_expression(old, OFFICIAL, ALIAS)
    hn, _ = harmonize_expression(new, OFFICIAL, ALIAS)
    assert list(ho.index) == list(hn.index)
    pd.testing.assert_frame_equal(ho, hn)

    legacy_sig, _ = harmonize_gene_list(["INDO", "BLR1"], OFFICIAL, ALIAS)
    modern_sig, _ = harmonize_gene_list(["IDO1", "CXCR5"], OFFICIAL, ALIAS)
    assert legacy_sig == modern_sig
    # and both are fully measured in both cohorts -- the realized size no
    # longer depends on vintage
    for mat in (ho, hn):
        assert sum(g in mat.index for g in legacy_sig) == 2
        assert sum(g in mat.index for g in modern_sig) == 2
