#!/usr/bin/env python3
"""Step 03 -- score every signature in every passing cohort.

Produces, for each signature x cohort:

* the per-sample score matrix (needed by the PCA / axis analyses and by
  the leave-one-cohort-out ceiling), written per cohort to
  ``results/scores/<cohort>_<method>.tsv.gz``;
* the response AUROC with its Hanley-McNeil standard error, and the
  overall-survival log-hazard-ratio per SD of score, in
  ``results/perf_per_sig_cohort.csv``.

Both scoring methods run: mean-z (primary) and ssGSEA (sensitivity).

A signature x cohort test is evaluable only when the cohort has at least
``evaluability.min_samples`` labeled samples with a finite score and both
response classes present; non-evaluable combinations are written with a
NaN AUROC and a reason, rather than dropped silently, so the denominator
of "1,463 signature x cohort tests" is auditable.

Usage:
  python scripts/03_score_signatures.py                 # both methods
  python scripts/03_score_signatures.py --method mean_z
  python scripts/03_score_signatures.py --response strict
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from icinull import (  # noqa: E402
    auroc,
    auroc_se_hanley_mcneil,
    load_cohort,
    load_config,
    load_signatures,
    logit,
    logit_se,
)
from icinull.scoring import mean_z_score, rank_matrix, ssgsea_score  # noqa: E402
from icinull.survival import cox_loghr_per_sd  # noqa: E402


def passing_cohorts() -> list:
    gate = ROOT / "results" / "cohort_gate.csv"
    if not gate.exists():
        raise SystemExit("results/cohort_gate.csv missing -- "
                         "run scripts/02_cohort_gate.py first.")
    g = pd.read_csv(gate)
    return list(g.loc[g.passes_gate, "cohort"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["mean_z", "ssgsea", "both"],
                    default="both")
    ap.add_argument("--response", choices=["primary", "strict"],
                    default="primary")
    args = ap.parse_args()

    cfg = load_config()
    sigs = load_signatures()
    methods = (["mean_z", "ssgsea"] if args.method == "both"
               else [args.method])

    ev = cfg["evaluability"]
    min_cov = cfg["scoring"]["min_genes_covered"]
    alpha = cfg["scoring"]["ssgsea_alpha"]

    score_dir = ROOT / "results" / "scores"
    score_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for cohort in passing_cohorts():
        c = load_cohort(cohort, config=cfg)
        y = c.labels(args.response, cfg)
        univ_expr = c.expr.loc[c.universe]
        print(f"== {cohort}: {univ_expr.shape[0]} expressed genes, "
              f"{y.size} labeled ({int((y == 1).sum())} R)")

        os_time = (c.clin["survival_time_os"]
                   if "survival_time_os" in c.clin.columns else None)
        os_event = (c.clin["event_occurred_os"]
                    if "event_occurred_os" in c.clin.columns else None)

        for method in methods:
            ranks = rank_matrix(univ_expr) if method == "ssgsea" else None
            mat = {}

            for name, spec in sigs.items():
                genes = spec["genes"]
                direction = int(spec["direction"])
                covered = [g for g in genes if g in c.expr.index]

                if len(covered) < min_cov:
                    rows.append(_row(cohort, name, spec, method, args.response,
                                     len(genes), len(covered),
                                     reason="insufficient_gene_coverage"))
                    continue

                if method == "mean_z":
                    s = mean_z_score(c.expr, covered, direction)
                else:
                    s = ssgsea_score(univ_expr, covered, direction,
                                     alpha=alpha, ranks=ranks)
                mat[name] = s

                idx = [i for i in y.index if i in s.index]
                sv = s.loc[idx]
                yv = y.loc[idx]
                ok = np.isfinite(sv.to_numpy())
                sv, yv = sv[ok], yv[ok]

                n_pos = int((yv == 1).sum())
                n_neg = int((yv == 0).sum())
                evaluable = (
                    sv.size >= ev["min_samples"]
                    and (not ev["require_both_classes"] or (n_pos and n_neg))
                )
                if not evaluable:
                    rows.append(_row(cohort, name, spec, method, args.response,
                                     len(genes), len(covered),
                                     n_pos=n_pos, n_neg=n_neg,
                                     reason="not_evaluable"))
                    continue

                a = auroc(sv.to_numpy(), yv.to_numpy())
                se = auroc_se_hanley_mcneil(a, n_pos, n_neg)

                loghr = hr_se = hr_p = float("nan")
                if os_time is not None and os_event is not None:
                    loghr, hr_se, hr_p = cox_loghr_per_sd(
                        s, os_time, os_event
                    )

                rows.append(_row(
                    cohort, name, spec, method, args.response,
                    len(genes), len(covered), n_pos=n_pos, n_neg=n_neg,
                    auroc=a, auroc_se=se,
                    logit_auroc=logit(a), logit_auroc_se=logit_se(a, se),
                    os_loghr_per_sd=loghr, os_loghr_se=hr_se, os_p=hr_p,
                    reason="",
                ))

            if mat:
                sm = pd.DataFrame(mat)
                sm.index.name = "sample_id"
                sm.to_csv(score_dir / f"{cohort}_{method}.tsv.gz", sep="\t")

    perf = pd.DataFrame(rows)
    suffix = "" if args.response == "primary" else f"_{args.response}"
    out = ROOT / "results" / f"perf_per_sig_cohort{suffix}.csv"
    perf.to_csv(out, index=False)

    done = perf[perf.auroc.notna()]
    print(f"\nsignature x cohort tests: {len(perf)} "
          f"({len(done)} evaluable, {len(perf) - len(done)} not)")
    for method, grp in done.groupby("method"):
        print(f"  {method:<8} n={len(grp):<5} mean AUROC {grp.auroc.mean():.3f}  "
              f"median {grp.auroc.median():.3f}  "
              f"frac<0.55 {(grp.auroc < 0.55).mean():.3f}")
    print(f"-> {out}")
    print(f"-> {score_dir}/<cohort>_<method>.tsv.gz")


def _row(cohort, name, spec, method, response, n_genes, n_covered, **kw) -> dict:
    row = {
        "cohort": cohort,
        "signature": name,
        "collection": spec.get("collection", ""),
        "direction": int(spec.get("direction", 1)),
        "method": method,
        "response_definition": response,
        "n_genes_published": n_genes,
        "n_genes_covered": n_covered,
        "coverage": round(n_covered / n_genes, 4) if n_genes else float("nan"),
        "n_pos": kw.get("n_pos", 0),
        "n_neg": kw.get("n_neg", 0),
        "auroc": kw.get("auroc", float("nan")),
        "auroc_se": kw.get("auroc_se", float("nan")),
        "logit_auroc": kw.get("logit_auroc", float("nan")),
        "logit_auroc_se": kw.get("logit_auroc_se", float("nan")),
        "os_loghr_per_sd": kw.get("os_loghr_per_sd", float("nan")),
        "os_loghr_se": kw.get("os_loghr_se", float("nan")),
        "os_p": kw.get("os_p", float("nan")),
        "excluded_reason": kw.get("reason", ""),
    }
    return row


if __name__ == "__main__":
    main()
