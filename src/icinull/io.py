"""Loading harmonized cohorts, the signature library, and config.

Everything downstream of ``R/01_harmonize.R`` reads through this module,
so there is exactly one definition of "the expressed universe", "the
response label", and "an evaluable signature x cohort test".
"""

from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

__all__ = [
    "Cohort",
    "load_config",
    "load_signatures",
    "list_cohorts",
    "load_cohort",
    "expressed_universe",
    "fetch_text",
    "repo_root",
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_config(path: str | Path | None = None) -> dict:
    p = Path(path) if path else repo_root() / "config" / "analysis_config.yaml"
    return yaml.safe_load(Path(p).read_text())


def load_signatures(path: str | Path | None = None) -> dict:
    p = Path(path) if path else repo_root() / "signatures" / "signatures.json"
    if not Path(p).exists():
        raise FileNotFoundError(
            f"{p} not found -- run scripts/01_build_signatures.py first."
        )
    return json.loads(Path(p).read_text())


# --------------------------------------------------------------------------
# cohorts
# --------------------------------------------------------------------------

@dataclass
class Cohort:
    """One harmonized cohort: expression, clinical, and derived labels."""

    name: str
    expr: pd.DataFrame          # genes x samples, log2-TPM
    clin: pd.DataFrame          # samples x fields
    universe: list              # expressed-gene universe (random-draw pool)

    def labels(self, definition: str = "primary", config: dict | None = None) -> pd.Series:
        """Binary response labels (1 = responder), indexed by sample.

        ``primary`` uses the compendium's R/NR label (R = CR/PR,
        NR = SD/PD). ``strict`` uses raw RECIST and discards stable
        disease (CR/PR vs PD), the manuscript's sensitivity analysis.
        """
        cfg = config or load_config()
        spec = cfg["response"][definition]

        if definition == "primary" and "response" in self.clin.columns:
            y = self.clin["response"].map({"R": 1, "NR": 0})
        else:
            if "recist" not in self.clin.columns:
                return pd.Series(dtype=float)
            mapping = {v: 1 for v in spec["responders"]}
            mapping.update({v: 0 for v in spec["non_responders"]})
            y = self.clin["recist"].map(mapping)

        return y.dropna().astype(int)

    @property
    def cancer_type(self) -> str:
        """Every distinct ``cancer_type`` value in the cohort, "|"-joined.

        Note that in some cohorts this column records the biopsy site
        rather than the primary tumour, so this is a superset of the
        cohort's disease. Use :attr:`primary_cancer_type` when counting
        cancer types.
        """
        if "cancer_type" not in self.clin.columns:
            return ""
        v = sorted({str(x) for x in self.clin["cancer_type"].dropna()})
        return "|".join(v)

    @property
    def primary_cancer_type(self) -> str:
        """The cohort's modal ``cancer_type`` — its primary tumour type."""
        if "cancer_type" not in self.clin.columns:
            return ""
        v = self.clin["cancer_type"].dropna().astype(str)
        return str(v.mode().iloc[0]) if len(v) else ""

    @property
    def treatment(self) -> str:
        if "treatment" not in self.clin.columns:
            return ""
        v = sorted({str(x) for x in self.clin["treatment"].dropna()})
        return "|".join(v)


def list_cohorts(harmonized_dir: str | Path | None = None) -> list:
    d = Path(harmonized_dir) if harmonized_dir else repo_root() / "data" / "harmonized"
    return sorted(p.name.replace("_expr.tsv.gz", "")
                  for p in Path(d).glob("*_expr.tsv.gz"))


def expressed_universe(expr: pd.DataFrame, min_detected_fraction: float = 0.20) -> list:
    """Genes detected in at least ``min_detected_fraction`` of samples.

    The compendium's log2-TPM matrices encode non-detection as the matrix
    floor (a large negative constant from the log transform of zero), so
    detection is "strictly above this cohort's floor" rather than a fixed
    numeric cutoff, which would not transfer across cohorts with different
    pseudocounts.
    """
    floor = float(np.nanmin(expr.to_numpy()))
    detected = expr.gt(floor).mean(axis=1)
    return list(expr.index[detected >= min_detected_fraction])


def load_cohort(
    name: str,
    harmonized_dir: str | Path | None = None,
    config: dict | None = None,
) -> Cohort:
    """Read one harmonized cohort from disk."""
    d = Path(harmonized_dir) if harmonized_dir else repo_root() / "data" / "harmonized"
    cfg = config or load_config()

    expr = pd.read_csv(d / f"{name}_expr.tsv.gz", sep="\t", index_col=0)
    clin = pd.read_csv(d / f"{name}_clin.tsv", sep="\t", index_col="sample_id",
                       low_memory=False)

    # Map the expression index onto current NCBI official symbols. The
    # cohorts were processed against different GENCODE vintages, so without
    # this a gene set written with legacy symbols matches in some cohorts and
    # vanishes in others -- see icinull.symbols for the full rationale.
    if cfg.get("expression", {}).get("harmonize_symbols", True):
        from .symbols import harmonize_expression, load_symbol_map
        official, alias = load_symbol_map()
        expr, _ = harmonize_expression(expr, official, alias)

    common = [s for s in expr.columns if s in clin.index]
    expr = expr[common]
    clin = clin.loc[common]

    universe = expressed_universe(
        expr, cfg["expression"]["expressed_gene_rule"]["min_detected_fraction"]
    )
    return Cohort(name=name, expr=expr, clin=clin, universe=universe)


# --------------------------------------------------------------------------
# small network helper (used only by the signature builder)
# --------------------------------------------------------------------------

def fetch_text(url: str, cache_path: str | Path, offline: bool = False) -> str:
    """GET ``url`` as text, caching to ``cache_path``."""
    cache = Path(cache_path)
    if cache.exists():
        return cache.read_text()
    if offline:
        raise SystemExit(f"--offline given but cache miss: {cache}")

    cache.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=60) as r:
        text = r.read().decode("utf-8")
    cache.write_text(text)
    return text
