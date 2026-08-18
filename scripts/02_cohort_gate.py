#!/usr/bin/env python3
"""Step 02 -- apply the pre-specified cohort inclusion gate.

Three criteria, fixed before any performance metric was computed:

  G1  >= 20 response-labeled samples
  G2  minority response class >= 5
  G3  >= 10,000 measured genes, so that published signatures and random
      gene sets are drawn from a comparable universe

Writes ``results/cohort_gate.csv`` -- every cohort with its criterion
values, pass/fail per criterion, and the failing criterion for exclusions.
This is the table the manuscript's Methods refers to as ``cohort_gate.csv``.

Usage:  python scripts/02_cohort_gate.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from icinull import list_cohorts, load_cohort, load_config  # noqa: E402


def main() -> None:
    cfg = load_config()
    gate = cfg["inclusion_gate"]
    rows = []

    cohorts = list_cohorts()
    if not cohorts:
        raise SystemExit(
            "No harmonized cohorts found -- run R/00_download_cohorts.R "
            "then R/01_harmonize.R first."
        )

    for name in cohorts:
        c = load_cohort(name, config=cfg)
        y = c.labels("primary", cfg)

        n_lab = int(y.size)
        n_r = int((y == 1).sum())
        n_nr = int((y == 0).sum())
        minority = min(n_r, n_nr)
        n_genes = int(len(c.universe))

        g1 = n_lab >= gate["min_response_labeled"]
        g2 = minority >= gate["min_minority_class"]
        g3 = n_genes >= gate["min_measured_genes"]

        failing = [g for g, ok in
                   (("G1", g1), ("G2", g2), ("G3", g3)) if not ok]

        os_col = "survival_time_os"
        n_os = (int(c.clin[os_col].notna().sum())
                if os_col in c.clin.columns else 0)

        rows.append({
            "cohort": name,
            "n_samples": int(c.expr.shape[1]),
            "n_response_labeled": n_lab,
            "n_responder": n_r,
            "n_nonresponder": n_nr,
            "responder_rate": round(n_r / n_lab, 4) if n_lab else float("nan"),
            "minority_class": minority,
            "n_expressed_genes": n_genes,
            "n_with_os": n_os,
            "cancer_type": c.cancer_type,
            "primary_cancer_type": getattr(c, "primary_cancer_type", ""),
            "treatment": c.treatment,
            "G1_min_labeled": g1,
            "G2_min_minority": g2,
            "G3_min_genes": g3,
            "passes_gate": not failing,
            "failing_criteria": ";".join(failing),
        })

    gate_df = pd.DataFrame(rows).sort_values(
        ["passes_gate", "n_response_labeled"], ascending=[False, False]
    )
    out = ROOT / "results" / "cohort_gate.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    gate_df.to_csv(out, index=False)

    passed = gate_df[gate_df.passes_gate]
    print(f"cohorts evaluated : {len(gate_df)}")
    print(f"cohorts passing   : {len(passed)}")
    print(f"samples (passing)  : {int(passed.n_samples.sum())}")
    print(f"response-labeled   : {int(passed.n_response_labeled.sum())}")
    print(f"with OS            : {int(passed.n_with_os.sum())}")
    if len(passed):
        pooled = passed.n_responder.sum() / passed.n_response_labeled.sum()
        print(f"pooled responder rate: {pooled:.3f}")
        # Count primary tumour types, not per-sample biopsy sites (see the
        # note in R/01_harmonize.R).
        col = ("primary_cancer_type"
               if "primary_cancer_type" in passed
               and passed.primary_cancer_type.notna().all()
               else "cancer_type")
        types = sorted({t for v in passed[col].astype(str)
                        for t in v.split("|") if t and t != "nan"})
        print(f"cancer types       : {len(types)} ({', '.join(types)})")

    excl = gate_df[~gate_df.passes_gate]
    if len(excl):
        print("\nexclusions:")
        for _, r in excl.iterrows():
            print(f"  {r.cohort:<16} {r.failing_criteria:<10} "
                  f"(labeled={r.n_response_labeled}, "
                  f"minority={r.minority_class}, genes={r.n_expressed_genes})")

    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
