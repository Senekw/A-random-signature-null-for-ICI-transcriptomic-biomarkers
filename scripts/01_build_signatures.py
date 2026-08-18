#!/usr/bin/env python3
"""Step 01 -- build the 102-signature library from upstream sources.

The library is *assembled*, never transcribed. Gene membership comes from
two maintained upstream repositories:

* the 12 canonical named ICI predictors, from the curated CSVs in
  ``bhklab/SignatureSets`` (the signature companion of the PredictioR
  compendium), plus three single-gene checkpoint markers;
* 90 immune / immuno-oncology gene sets from IOBR's
  ``signature_collection``, restricted to the roster in
  ``config/signature_roster_iobr.txt``.

Writes ``signatures/signatures.json`` (name -> genes, direction, source)
and ``signatures/signature_provenance.csv`` (one row per signature with
its upstream path, PMID and size), which are the two artifacts the
manuscript's Methods names.

Requires network access to raw.githubusercontent.com. Re-running is
idempotent; pass ``--offline`` to rebuild from the cached
``signatures/.cache`` payloads.

Usage:  python scripts/01_build_signatures.py [--offline]
"""

from __future__ import annotations

import argparse
import io
import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]

# The R interpreter used for the two steps that read IOBR's .rda objects.
# Overridable because a conda/renv R is commonly not the one on PATH:
#   RSCRIPT=/path/to/Rscript python scripts/01_build_signatures.py
RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript") or "Rscript"

sys.path.insert(0, str(ROOT / "src"))

from icinull.io import fetch_text  # noqa: E402

SIGSETS_RAW = "https://raw.githubusercontent.com/bhklab/SignatureSets/main"
IOBR_RAW = "https://raw.githubusercontent.com/IOBR/IOBR/master"
IOBR_RDA_URL = f"{IOBR_RAW}/data/signature_collection.rda"

OUT_DIR = ROOT / "signatures"
CACHE = OUT_DIR / ".cache"


# --------------------------------------------------------------------------
# canonical named predictors
# --------------------------------------------------------------------------

def build_canonical(offline: bool) -> tuple[dict, list]:
    spec = yaml.safe_load(
        (ROOT / "config" / "canonical_signatures.yaml").read_text()
    )
    sigs: dict = {}
    prov: list = []

    for entry in spec:
        name = entry["name"]

        if entry["source"] == "literal":
            genes = list(entry["genes"])
            upstream = "literal (single-gene marker)"
        else:
            url = f"{SIGSETS_RAW}/{entry['path']}"
            text = fetch_text(
                url,
                CACHE / f"sigsets__{Path(entry['path']).name}",
                offline=offline,
            )
            df = pd.read_csv(io.StringIO(text))
            col = "gene_name" if "gene_name" in df.columns else df.columns[0]
            genes = sorted({str(g) for g in df[col].dropna() if str(g).strip()})
            upstream = f"bhklab/SignatureSets:{entry['path']}"

        sigs[name] = {
            "genes": genes,
            "direction": int(entry["direction"]),
            "collection": "canonical",
            "label": entry.get("label", name),
        }
        prov.append({
            "signature": name,
            "label": entry.get("label", name),
            "collection": "canonical",
            "n_genes": len(genes),
            "direction": int(entry["direction"]),
            "pmid": entry.get("pmid", ""),
            "reference": entry.get("reference", ""),
            "upstream_source": upstream,
        })

    return sigs, prov


# --------------------------------------------------------------------------
# IOBR curated gene sets
# --------------------------------------------------------------------------

# Sets whose high scores mark resistance, stroma, EMT, TGF-beta,
# suppressive myeloid content or proliferation get direction -1. This is an
# a-priori assignment on each set's stated biology, fixed before any
# performance metric was computed.
NEGATIVE_PATTERNS = (
    "EMT", "TGF", "MDSC", "TAM", "CAF", "Treg", "T_cell_exhaustion",
    "T_cell_regulatory", "ICB_resistance", "Exhausted", "Cell_cycle",
    "CellCycle", "DNA_replication", "Histones", "WNT", "FGFR3",
    "Neutrophils", "Macrophages", "Mast_cells", "TMEscoreB",
    "Nucleotide_excision_repair", "Base_excision_repair",
    "Mismatch_Repair", "Homologous_recombination", "DDR",
    "SW480_cancer_cells", "Normal_mucosa", "Lymph_vessels",
    "Pan_F_TBRs", "GPAGs", "PPAGs",
)


def iobr_direction(name: str) -> int:
    return -1 if any(p.lower() in name.lower() for p in NEGATIVE_PATTERNS) else 1


def build_iobr(offline: bool) -> tuple[dict, list]:
    """Read the IOBR collection via a short R call, then apply the roster."""
    roster = [
        ln.strip()
        for ln in (ROOT / "config" / "signature_roster_iobr.txt")
        .read_text().splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]

    payload = CACHE / "iobr_signature_collection.json"
    if not payload.exists():
        if offline:
            raise SystemExit(
                f"--offline given but {payload} is absent; run once online."
            )
        rda = CACHE / "signature_collection.rda"
        rda.parent.mkdir(parents=True, exist_ok=True)
        if not rda.exists():
            urllib.request.urlretrieve(IOBR_RDA_URL, rda)

        r_code = (
            f'e <- new.env(); load("{rda.as_posix()}", envir = e); '
            f'sc <- e$signature_collection; '
            f'cat(jsonlite::toJSON(sc, auto_unbox = FALSE), '
            f'file = "{payload.as_posix()}")'
        )
        subprocess.run([RSCRIPT, "-e", r_code], check=True,
                       capture_output=True, text=True)

    coll = json.loads(payload.read_text())

    missing = [n for n in roster if n not in coll]
    if missing:
        raise SystemExit(
            "Roster names absent from IOBR signature_collection "
            f"({len(missing)}): {missing[:10]}"
        )

    sigs: dict = {}
    prov: list = []
    for name in roster:
        genes = sorted({str(g) for g in coll[name] if str(g).strip()})
        d = iobr_direction(name)
        sigs[name] = {
            "genes": genes,
            "direction": d,
            "collection": "iobr",
            "label": name.replace("_", " "),
        }
        prov.append({
            "signature": name,
            "label": name.replace("_", " "),
            "collection": "iobr",
            "n_genes": len(genes),
            "direction": d,
            "pmid": "",
            "reference": "Zeng D, et al. IOBR. Front Immunol 2021;12:687975.",
            "upstream_source": "IOBR/IOBR:data/signature_collection.rda",
        })

    return sigs, prov


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true",
                    help="rebuild from signatures/.cache without network")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    canon, canon_prov = build_canonical(args.offline)
    iobr, iobr_prov = build_iobr(args.offline)

    overlap = set(canon) & set(iobr)
    if overlap:
        # A canonical hand-encoded definition wins over the IOBR copy of
        # the same name; the collision is reported, not silently resolved.
        print(f"note: {len(overlap)} name(s) in both collections, "
              f"canonical kept: {sorted(overlap)}")
        for k in overlap:
            iobr.pop(k)
        iobr_prov = [r for r in iobr_prov if r["signature"] not in overlap]

    sigs = {**canon, **iobr}
    prov = pd.DataFrame(canon_prov + iobr_prov)

    (OUT_DIR / "signatures.json").write_text(json.dumps(sigs, indent=2))
    prov.to_csv(OUT_DIR / "signature_provenance.csv", index=False)

    print(f"signature library: {len(sigs)} signatures "
          f"({len(canon)} canonical + {len(iobr)} IOBR)")
    print(f"  gene-set size: median {int(prov.n_genes.median())}, "
          f"range {prov.n_genes.min()}-{prov.n_genes.max()}")
    print(f"  direction -1 (resistance/stroma): {int((prov.direction < 0).sum())}")
    print(f"  -> {OUT_DIR / 'signatures.json'}")
    print(f"  -> {OUT_DIR / 'signature_provenance.csv'}")

    if len(sigs) != 102:
        print(f"WARNING: expected 102 signatures, built {len(sigs)}. "
              "Check config/signature_roster_iobr.txt against the upstream "
              "collection -- IOBR may have renamed or removed a set.")


if __name__ == "__main__":
    main()
