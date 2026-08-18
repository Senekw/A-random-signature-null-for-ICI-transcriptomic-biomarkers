#!/usr/bin/env python3
"""Step 08 -- collect every number quoted in the manuscript.

Reads the result tables produced by steps 02-07 and writes
``results/headline_numbers.json``: one machine-readable record per claim
in the paper, each carrying the value, the table it came from, and the
column that produced it.

The point is auditability. Any number in the manuscript should be
traceable to a row of a committed CSV via this file, and re-running the
pipeline should regenerate it. Nothing here is hard-coded: if a table is
missing, the corresponding entry is null rather than a remembered value.

Usage:  python scripts/08_headline_numbers.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
sys.path.insert(0, str(ROOT / "src"))


def read(name: str):
    p = RES / name
    return pd.read_csv(p) if p.exists() else None


def num(x):
    """JSON-safe float (NaN -> None)."""
    if x is None:
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if not np.isfinite(v) else round(v, 6)


def main() -> None:
    out: dict = {"provenance": {}, "cohorts": {}, "performance": {},
                 "null_calibration": {}, "single_axis": {},
                 "confound": {}, "meta_analysis": {}, "ceiling": {}}

    # ---- cohorts -----------------------------------------------------------
    gate = read("cohort_gate.csv")
    if gate is not None:
        p = gate[gate.passes_gate]
        out["cohorts"] = {
            "n_evaluated": int(len(gate)),
            "n_passing_gate": int(len(p)),
            "n_excluded": int((~gate.passes_gate).sum()),
            "exclusions": {r.cohort: r.failing_criteria
                           for _, r in gate[~gate.passes_gate].iterrows()},
            "n_samples": int(p.n_samples.sum()),
            "n_response_labeled": int(p.n_response_labeled.sum()),
            "n_with_os": int(p.n_with_os.sum()),
            "pooled_responder_rate": num(
                p.n_responder.sum() / p.n_response_labeled.sum()
                if p.n_response_labeled.sum() else None
            ),
            "cancer_types": sorted({t for v in p.cancer_type
                                    for t in str(v).split("|") if t}),
            "treatment_classes": sorted({t for v in p.treatment
                                         for t in str(v).split("|") if t}),
        }
        out["cohorts"]["n_cancer_types"] = len(out["cohorts"]["cancer_types"])
        out["cohorts"]["n_treatment_classes"] = len(
            out["cohorts"]["treatment_classes"]
        )
        out["provenance"]["cohorts"] = "results/cohort_gate.csv"

    # ---- per-signature performance ----------------------------------------
    perf = read("perf_per_sig_cohort.csv")
    if perf is not None:
        mz = perf[(perf.method == "mean_z") & perf.auroc.notna()]
        ss = perf[(perf.method == "ssgsea") & perf.auroc.notna()]
        by_sig = mz.groupby("signature").auroc.mean().sort_values(ascending=False)

        out["performance"] = {
            "n_signatures": int(perf.signature.nunique()),
            "n_tests_mean_z": int(len(mz)),
            "mean_auroc": num(mz.auroc.mean()),
            "median_auroc": num(mz.auroc.median()),
            "frac_tests_below_0.55": num((mz.auroc < 0.55).mean()),
            "median_gene_coverage": num(mz.coverage.median()),
            "best_signatures_mean_auroc": {
                s: num(v) for s, v in by_sig.head(10).items()
            },
        }
        if len(ss):
            j = (mz.groupby("signature").auroc.mean()
                 .to_frame("mean_z")
                 .join(ss.groupby("signature").auroc.mean().rename("ssgsea"))
                 .dropna())
            out["performance"]["scoring_method_pearson_r"] = num(
                j.mean_z.corr(j.ssgsea)
            )
            out["performance"]["scoring_method_spearman_rho"] = num(
                j.mean_z.corr(j.ssgsea, method="spearman")
            )
        out["provenance"]["performance"] = "results/perf_per_sig_cohort.csv"

    # ---- the null ----------------------------------------------------------
    null = read("null_results.csv")
    pool = read("null_pooled.csv")
    if null is not None:
        canon = null[null.collection == "canonical"]
        by_cohort = null.groupby("cohort").beats_null_p05.mean()
        out["null_calibration"] = {
            "n_signature_cohort_tests": int(len(null)),
            "n_draws_per_test": int(null.n_draws.max()),
            "seed": int(null.seed.min()) if "seed" in null else None,
            "frac_beating_null_p05": num(null.beats_null_p05.mean()),
            "frac_canonical_beating_null": num(canon.beats_null_p05.mean())
            if len(canon) else None,
            "median_observed_auroc": num(null.observed_auroc.median()),
            "median_null_q95": num(null.null_q95.median()),
            "mean_observed_auroc": num(null.observed_auroc.mean()),
            "mean_null_auroc": num(null.null_mean.mean()),
            "mean_advantage_over_null": num(
                (null.observed_auroc - null.null_mean).mean()
            ),
            "null_beating_rate_by_cohort": {
                k: num(v) for k, v in by_cohort.sort_values().items()
            },
            "null_beating_range": [num(by_cohort.min()), num(by_cohort.max())],
        }
        if pool is not None:
            out["null_calibration"]["n_signatures_fdr05_positive"] = int(
                pool.beats_null_fdr05.sum()
            )
            out["null_calibration"]["n_signatures_pooled"] = int(len(pool))
        out["provenance"]["null_calibration"] = "results/null_results.csv"

    # ---- single axis -------------------------------------------------------
    for method in ("mean_z", "ssgsea"):
        pca = read(f"axis_pca_{method}.csv")
        ven = read(f"venet_alignment_{method}.csv")
        if pca is None and ven is None:
            continue
        block: dict = {}
        if pca is not None:
            block.update({
                "pc1_var_explained_mean": num(pca.pc1_var_explained.mean()),
                "pc1_var_explained_max": num(pca.pc1_var_explained.max()),
                "pc1_corr_gep_mean": num(pca.pc1_corr_gep.mean()),
            })
        if ven is not None and len(ven) > 2:
            from scipy.stats import pearsonr
            r, p = pearsonr(ven.mean_corr_with_gep, ven.mean_auroc)
            block.update({
                "venet_r": num(r),
                "venet_r2": num(r ** 2),
                "venet_p": None if not np.isfinite(p) else float(f"{p:.3g}"),
                "venet_n_signatures": int(len(ven)),
            })
        out["single_axis"][method] = block
        out["provenance"][f"single_axis.{method}"] = (
            f"results/axis_pca_{method}.csv, results/venet_alignment_{method}.csv"
        )

    # ---- infiltration / purity confound ------------------------------------
    pca = read("axis_pca_mean_z.csv")
    part = read("partialled_auroc_mean_z.csv")
    est = read("estimate_scores.csv")
    if pca is not None:
        out["confound"]["axis_corr_estimate_immune_mean"] = num(
            pca.pc1_corr_estimate_immune.mean()
        )
        out["confound"]["axis_corr_tumor_purity_mean"] = num(
            pca.pc1_corr_tumor_purity.mean()
        )
    if part is not None:
        pm = part.groupby("signature").agg(
            raw=("auroc_raw", "mean"), part=("auroc_partialled", "mean")
        ).dropna()
        out["confound"].update({
            "mean_auroc_raw": num(pm.raw.mean()),
            "mean_auroc_after_removing_infiltration": num(pm.part.mean()),
            "best10_mean_auroc": num(pm.raw.nlargest(10).mean()),
        })
    if est is not None and part is not None:
        out["provenance"]["confound"] = (
            "results/partialled_auroc_mean_z.csv, results/estimate_scores.csv"
        )

    # ---- meta-analysis -----------------------------------------------------
    meta = read("meta_results.csv")
    if meta is not None:
        ok = meta[meta.ci_low.notna()]
        if len(ok):
            b = ok.iloc[0]
            out["meta_analysis"] = {
                "n_meta_analyzable": int(len(ok)),
                "best_signature": str(b.signature),
                "best_pooled_auroc": num(b.pooled_auroc),
                "best_ci": [num(b.ci_low), num(b.ci_high)],
                "n_ci_excluding_0.5": int(ok.ci_excludes_half.sum()),
                "median_I2": num(ok.I2.median()),
                "median_tau2": num(ok.tau2.median()),
            }
            out["provenance"]["meta_analysis"] = "results/meta_results.csv"

    # ---- transferable ceiling ---------------------------------------------
    ceil = read("ceiling_loco.csv")
    opt = read("optimism.csv")
    if ceil is not None:
        out["ceiling"]["loco_auroc"] = {
            str(r.predictor): num(r.loco_auroc) for _, r in ceil.iterrows()
        }
        gep = ceil[ceil.predictor == "gep_axis"]
        if len(gep) and np.isfinite(gep.ci_low.iloc[0]):
            out["ceiling"]["gep_axis_ci"] = [num(gep.ci_low.iloc[0]),
                                             num(gep.ci_high.iloc[0])]
        out["provenance"]["ceiling"] = "results/ceiling_loco.csv"
    if opt is not None:
        d = dict(zip(opt.quantity, opt.auroc))
        out["ceiling"]["in_cohort_best_selected"] = num(
            d.get("best_signature_selected_in_cohort")
        )
        out["ceiling"]["within_cohort_cv_panel"] = num(
            d.get("within_cohort_cv_panel")
        )
        transferable = d.get("transferable_ceiling_gep_axis")
        best = d.get("best_signature_selected_in_cohort")
        if transferable is not None and best is not None:
            out["ceiling"]["optimism_gap"] = num(best - transferable)

    dest = RES / "headline_numbers.json"
    dest.write_text(json.dumps(out, indent=2, sort_keys=False))

    filled = sum(1 for k, v in out.items()
                 if k != "provenance" and isinstance(v, dict) and v)
    print(f"headline_numbers.json written ({filled}/7 result blocks "
          f"populated)")
    for block in ("cohorts", "performance", "null_calibration",
                  "single_axis", "confound", "meta_analysis", "ceiling"):
        state = "ok" if out.get(block) else "EMPTY (upstream step not run)"
        print(f"  {block:<18} {state}")
    print(f"-> {dest}")


if __name__ == "__main__":
    main()
