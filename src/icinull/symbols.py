"""Harmonize gene symbols to one vocabulary across cohorts and signatures.

Why this exists
---------------
The compendium's cohorts were processed against different GENCODE vintages,
so the *same gene* appears under different symbols in different cohorts.
Concretely, in the 15 gate-passing cohorts:

* ``ICB_Van_Allen``, ``ICB_Miao1``, ``ICB_Braun`` and ``ICB_Puch`` carry the
  pre-2020 histone symbols (``HIST1H2AG``);
* ``ICB_Kim``, ``ICB_Gide``, ``ICB_Riaz`` and others carry the current ones
  (``H2AC11``);

and roughly 300-2,400 retired symbols per cohort are still in use.

Left alone this silently corrupts every cross-cohort comparison. A gene set
written with legacy symbols matches in the four old-annotation cohorts and
vanishes in the other eleven, so its realized size — and therefore the
size-matched null, the AUROC, and its weight in the meta-analysis — depends
on which annotation vintage a cohort happened to be processed with rather
than on biology. It also silently shrinks gene sets: 80 of the 2,635
distinct signature genes matched *no* cohort at all before this module
existed, affecting 25 signatures.

What it does
------------
Both sides of the join are mapped onto one vocabulary: current NCBI official
symbols. The mapping is built from the authoritative
``Homo_sapiens.gene_info`` table and is deliberately conservative:

* a symbol that is *already* an official symbol is never rewritten;
* an alias is rewritten only when it maps to exactly **one** official
  symbol. Ambiguous aliases are left alone and reported — ``IL8RA`` is a
  documented synonym of both ``CXCR1`` and ``CXCR2``, and guessing would
  fabricate membership.

Collapsing two rows onto one symbol uses the same rule as the rest of the
pipeline: keep the row with the highest mean expression.
"""

from __future__ import annotations

import gzip
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = [
    "GENE_INFO_URL",
    "load_symbol_map",
    "harmonize_index",
    "harmonize_expression",
    "harmonize_gene_list",
]

GENE_INFO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/"
    "Homo_sapiens.gene_info.gz"
)


def _default_cache() -> Path:
    return (Path(__file__).resolve().parents[2] / "signatures" / ".cache"
            / "Homo_sapiens.gene_info.gz")


_MAP_CACHE: dict = {}


def load_symbol_map(cache_path: str | Path | None = None,
                    offline: bool = False) -> tuple:
    """Return ``(official, alias_to_official)``.

    ``official`` is the set of current symbols (so they are never
    rewritten); ``alias_to_official`` maps each unambiguous retired symbol
    to its current one.
    """
    cache = Path(cache_path) if cache_path else _default_cache()
    key = str(cache)
    if key in _MAP_CACHE:            # parsing the table takes ~1s; do it once
        return _MAP_CACHE[key]
    if not cache.exists():
        if offline:
            raise SystemExit(
                f"gene-info table not cached at {cache} and --offline was "
                f"given. Fetch it once with:\n  curl -sfL {GENE_INFO_URL} "
                f"-o {cache}"
            )
        cache.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(GENE_INFO_URL, cache)

    official: set = set()
    syn: dict = {}
    with gzip.open(cache, "rt") as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        col = {k: i for i, k in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sym = f[col["Symbol"]]
            auth = f[col["Symbol_from_nomenclature_authority"]]
            best = auth if auth not in ("-", "") else sym
            official.add(best)
            official.add(sym)
            for s in f[col["Synonyms"]].split("|"):
                if s and s != "-":
                    syn.setdefault(s, set()).add(best)

    alias = {a: next(iter(t)) for a, t in syn.items()
             if len(t) == 1 and a not in official}
    _MAP_CACHE[key] = (official, alias)
    return official, alias


def harmonize_index(symbols, official: set, alias: dict) -> list:
    """Map a list of symbols onto the current vocabulary."""
    return [s if (s in official or s not in alias) else alias[s]
            for s in symbols]


def harmonize_expression(expr: pd.DataFrame, official: set, alias: dict
                         ) -> tuple:
    """Rewrite an expression matrix's index to current symbols.

    Returns ``(expr_harmonized, report)`` where ``report`` is a dict of
    counts: how many symbols were rewritten, and how many rows were
    collapsed as a result.
    """
    old = list(expr.index)
    new = harmonize_index(old, official, alias)
    n_renamed = sum(1 for a, b in zip(old, new) if a != b)

    out = expr.copy()
    out.index = new

    n_before = len(out)
    if out.index.duplicated().any():
        # same collapse rule as harmonization: highest mean expression wins
        order = np.argsort(-out.mean(axis=1, skipna=True).to_numpy(),
                          kind="stable")
        out = out.iloc[order]
        out = out[~out.index.duplicated(keep="first")]
        out = out.sort_index()

    return out, {
        "n_symbols_in": n_before,
        "n_renamed": int(n_renamed),
        "n_collapsed": int(n_before - len(out)),
        "n_symbols_out": int(len(out)),
    }


def harmonize_gene_list(genes, official: set, alias: dict) -> tuple:
    """Map signature genes onto current symbols.

    Returns ``(harmonized, renamed)`` where ``renamed`` is the
    ``{old: new}`` subset that changed, for provenance.
    """
    out, renamed = [], {}
    for g in dict.fromkeys(genes):
        h = g if (g in official or g not in alias) else alias[g]
        if h != g:
            renamed[g] = h
        out.append(h)
    return list(dict.fromkeys(out)), renamed
