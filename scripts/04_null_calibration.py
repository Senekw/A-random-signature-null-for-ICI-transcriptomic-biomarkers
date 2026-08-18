#!/usr/bin/env python3
"""Step 04 -- calibrate every signature against size-matched random sets.

This is the paper's central test. For each signature in each passing
cohort, draw ``null_model.n_draws`` gene sets of the signature's realized
(coverage-limited) size from that cohort's expressed universe, score them
by mean-z, and compute the one-sided empirical p-value that a random set
matches or exceeds the published signature's AUROC.

Per-signature evidence is then pooled across cohorts by Stouffer's method
on the per-cohort z-versus-null statistics, with Benjamini-Hochberg FDR
control.

Writes ``results/null_results.csv`` (one row per signature x cohort) and
``results/null_pooled.csv`` (one row per signature).

This is the expensive step: n_signatures x n_cohorts x n_draws AUROC
evaluations. It parallelizes across cohorts.

Usage:
  python scripts/04_null_calibration.py
  python scripts/04_null_calibration.py --draws 100        # fast smoke run
  python scripts/04_null_calibration.py --seed 1           # seed stability
  python scripts/04_null_calibration.py --jobs 4
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from icinull import (  # noqa: E402
    load_cohort,
    load_config,
    load_signatures,
    random_null_test,
    stouffer,
)


def passing_cohorts() -> list:
    g = pd.read_csv(ROOT / "results" / "cohort_gate.csv")
    return list(g.loc[g.passes_gate, "cohort"])


def run_cohort(cohort: str, draws: int, seed: int, response: str) -> list:
    """All signatures against the null in one cohort."""
    cfg = load_config()
    sigs = load_signatures()
    c = load_cohort(cohort, config=cfg)
    y = c.labels(response, cfg)
    min_cov = cfg["scoring"]["min_genes_covered"]
    ev = cfg["evaluability"]

    rows = []
    for i, (name, spec) in enumerate(sorted(sigs.items())):
        covered = [g for g in spec["genes"] if g in c.expr.index]
        if len(covered) < min_cov:
            continue

        # Per-signature seed offset: every signature x cohort draw is
        # reproducible on its own, and no two share a random stream.
        sig_seed = seed + 1000 * i + abs(hash(cohort)) % 997

        res = random_null_test(
            c.expr, y, covered,
            direction=int(spec["direction"]),
            universe=c.universe,
            n_draws=draws,
            seed=sig_seed,
        )
        d = res.as_dict()
        if (d["n_pos"] + d["n_neg"]) < ev["min_samples"]:
            continue

        d.update({
            "cohort": cohort,
            "signature": name,
            "collection": spec.get("collection", ""),
            "direction": int(spec["direction"]),
            "seed": sig_seed,
            "beats_null_p05": bool(np.isfinite(d["p_empirical"])
                                   and d["p_empirical"] < 0.05),
            "above_null_q95": bool(np.isfinite(d["observed_auroc"])
                                   and np.isfinite(d["null_q95"])
                                   and d["observed_auroc"] > d["null_q95"]),
        })
        rows.append(d)

    print(f"  {cohort}: {len(rows)} tests, "
          f"{sum(r['beats_null_p05'] for r in rows)} beat null at p<0.05")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=None,
                    help="random sets per signature (default: config)")
    ap.add_argument("--seed", type=int, default=None,
                    help="RNG seed (default: config; 20260705 as published)")
    ap.add_argument("--response", choices=["primary", "strict"],
                    default="primary")
    ap.add_argument("--jobs", type=int, default=1,
                    help="cohorts to process in parallel")
    args = ap.parse_args()

    cfg = load_config()
    draws = args.draws or cfg["null_model"]["n_draws"]
    seed = args.seed if args.seed is not None else cfg["seed"]
    cohorts = passing_cohorts()

    print(f"null calibration: {len(cohorts)} cohorts x {draws} draws, "
          f"seed {seed}, response={args.response}")

    rows: list = []
    parallel_ok = args.jobs > 1
    if parallel_ok:
        try:
            with ProcessPoolExecutor(max_workers=args.jobs) as ex:
                futs = {ex.submit(run_cohort, c, draws, seed, args.response): c
                        for c in cohorts}
                for f in as_completed(futs):
                    rows.extend(f.result())
        except (OSError, PermissionError, NotImplementedError) as exc:
            # Some sandboxes and container runtimes deny the POSIX
            # semaphores multiprocessing needs. Fall back rather than fail:
            # the result is identical, only slower.
            print(f"note: parallel execution unavailable ({exc}); "
                  f"running serially")
            rows, parallel_ok = [], False

    if not parallel_ok:
        for c in cohorts:
            rows.extend(run_cohort(c, draws, seed, args.response))

    null = pd.DataFrame(rows)
    tag = ""
    if args.response != "primary":
        tag += f"_{args.response}"
    if args.seed is not None and args.seed != cfg["seed"]:
        tag += f"_seed{args.seed}"

    out = ROOT / "results" / f"null_results{tag}.csv"
    null.to_csv(out, index=False)

    # ---- pool across cohorts per signature --------------------------------
    pooled = []
    for sig, grp in null.groupby("signature"):
        z, p = stouffer(grp.z_vs_null.to_numpy())
        pooled.append({
            "signature": sig,
            "collection": grp.collection.iloc[0],
            "n_cohorts": int(grp.shape[0]),
            "mean_observed_auroc": float(grp.observed_auroc.mean()),
            "mean_null_auroc": float(grp.null_mean.mean()),
            "mean_advantage": float((grp.observed_auroc - grp.null_mean).mean()),
            "n_cohorts_beating_null": int(grp.beats_null_p05.sum()),
            "frac_cohorts_beating_null": float(grp.beats_null_p05.mean()),
            "stouffer_z": z,
            "stouffer_p": p,
        })
    pool = pd.DataFrame(pooled)

    from statsmodels.stats.multitest import multipletests
    ok = pool.stouffer_p.notna()
    pool.loc[ok, "fdr_q"] = multipletests(
        pool.loc[ok, "stouffer_p"], method="fdr_bh"
    )[1]
    pool["beats_null_fdr05"] = (
        (pool.fdr_q < 0.05) & (pool.mean_advantage > 0)
    )
    pool = pool.sort_values("mean_observed_auroc", ascending=False)
    pool.to_csv(ROOT / "results" / f"null_pooled{tag}.csv", index=False)

    # ---- headline summary --------------------------------------------------
    n_tests = len(null)
    beat = int(null.beats_null_p05.sum())
    canon = null[null.collection == "canonical"]
    med_obs = float(null.observed_auroc.median())
    med_q95 = float(null.null_q95.median())

    print(f"\nsignature x cohort tests      : {n_tests}")
    print(f"beat size-matched null p<0.05 : {beat} ({beat / n_tests:.1%})")
    if len(canon):
        print(f"canonical named predictors    : "
              f"{canon.beats_null_p05.mean():.1%} of cohorts")
    print(f"median observed AUROC         : {med_obs:.3f}")
    print(f"median null 95th percentile   : {med_q95:.3f}")
    print(f"mean observed vs null         : "
          f"{null.observed_auroc.mean():.3f} vs {null.null_mean.mean():.3f}")
    print(f"signatures at FDR<0.05 (+ eff): "
          f"{int(pool.beats_null_fdr05.sum())} / {len(pool)}")

    print("\nnull-beating rate by cohort (a cohort property, not a "
          "signature one):")
    by_cohort = (null.groupby("cohort").beats_null_p05.mean()
                 .sort_values(ascending=False))
    for coh, frac in by_cohort.items():
        print(f"  {coh:<16} {frac:.2%}")

    print(f"\n-> {out}")
    print(f"-> {ROOT / 'results' / f'null_pooled{tag}.csv'}")


if __name__ == "__main__":
    main()
