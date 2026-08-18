#!/usr/bin/env python3
"""Verify the IOBR signature roster against the upstream collection.

The 90 IOBR gene sets are named in ``config/signature_roster_iobr.txt``
as the union of the ``sig_group`` categories listed in that file's header,
restricted to sets present in ``signature_collection`` and disjoint from
the 12 canonical named predictors.

This script re-derives that union from the current upstream IOBR release
and reports any drift: sets IOBR has renamed, removed, or added to those
categories since the analysis was run. It does not silently update the
roster -- a change in the roster changes the analysis, so it is a decision
for a human, recorded as a deviation from the pre-specification.

Requires R with the IOBR data objects; both are fetched from GitHub.

Usage:  python scripts/verify_roster.py
"""

from __future__ import annotations

import re
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# The R interpreter used for the two steps that read IOBR's .rda objects.
# Overridable because a conda/renv R is commonly not the one on PATH:
#   RSCRIPT=/path/to/Rscript python scripts/01_build_signatures.py
RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript") or "Rscript"

ROSTER = ROOT / "config" / "signature_roster_iobr.txt"
CANON = ROOT / "config" / "canonical_signatures.yaml"

IOBR_RAW = "https://raw.githubusercontent.com/IOBR/IOBR/master/data"

R_SCRIPT = r"""
args <- commandArgs(trailingOnly = TRUE)
sc_file <- args[1]; sg_file <- args[2]; groups <- args[-(1:2)]
e <- new.env()
load(sc_file, envir = e); load(sg_file, envir = e)
sc <- e$signature_collection; sg <- e$sig_group
have <- intersect(groups, names(sg))
missing_groups <- setdiff(groups, names(sg))
if (length(missing_groups))
  cat("MISSING_GROUP", paste(missing_groups, collapse = ","), "\n")
u <- sort(intersect(unique(unlist(sg[have])), names(sc)))
cat("DERIVED", paste(u, collapse = ","), "\n")
"""


def main() -> None:
    if not ROSTER.exists():
        raise SystemExit(f"roster not found: {ROSTER}")

    text = ROSTER.read_text()
    roster = [ln.strip() for ln in text.splitlines()
              if ln.strip() and not ln.lstrip().startswith("#")]

    # the categories are recorded in the header, one per commented line
    # after the "categories:" sentence
    groups = [ln.lstrip("# ").strip() for ln in text.splitlines()
              if re.match(r"^#\s{2,}\w", ln)]
    groups = [g for g in groups if g and " " not in g]

    print(f"roster           : {len(roster)} signatures")
    print(f"source categories: {len(groups)} ({', '.join(groups)})")

    canon_names: list = []
    if CANON.exists():
        canon_names = re.findall(r"^- name:\s*(\S+)", CANON.read_text(),
                                 flags=re.M)
    overlap = sorted(set(roster) & set(canon_names))
    print(f"canonical overlap: {len(overlap)}"
          + (f"  {overlap}" if overlap else "  (correct: disjoint)"))

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for f in ("signature_collection.rda", "sig_group.rda"):
            urllib.request.urlretrieve(f"{IOBR_RAW}/{f}", td / f)
        rs = td / "verify.R"
        rs.write_text(R_SCRIPT)
        try:
            out = subprocess.run(
                [RSCRIPT, str(rs), str(td / "signature_collection.rda"),
                 str(td / "sig_group.rda"), *groups],
                capture_output=True, text=True, check=True,
            ).stdout
        except FileNotFoundError:
            raise SystemExit("Rscript not on PATH -- needed to read the "
                             "IOBR .rda objects.")
        except subprocess.CalledProcessError as e:
            raise SystemExit(f"R failed:\n{e.stderr[-800:]}")

    for line in out.splitlines():
        if line.startswith("MISSING_GROUP"):
            print(f"WARNING: sig_group categories no longer in IOBR: "
                  f"{line.split(None, 1)[1]}")

    derived = []
    for line in out.splitlines():
        if line.startswith("DERIVED"):
            derived = [s for s in line.split(None, 1)[1].strip().split(",") if s]
    derived = sorted(set(derived) - set(canon_names))

    print(f"upstream derived : {len(derived)} signatures")

    removed = sorted(set(roster) - set(derived))
    added = sorted(set(derived) - set(roster))

    if not removed and not added:
        print("\nOK: roster matches the current upstream collection exactly.")
        return

    if removed:
        print(f"\nin roster, NOT in current upstream ({len(removed)}) -- "
              f"renamed or removed by IOBR:")
        for s in removed:
            print(f"  - {s}")
    if added:
        print(f"\nin current upstream, NOT in roster ({len(added)}) -- "
              f"added by IOBR since the analysis:")
        for s in added:
            print(f"  + {s}")

    print("\nThe roster is part of the pre-specification. Do not edit it to "
          "match upstream without recording the change as a deviation in "
          "docs/methods_prespec.md -- the committed roster is what the "
          "published results used.")
    sys.exit(1)


if __name__ == "__main__":
    main()
