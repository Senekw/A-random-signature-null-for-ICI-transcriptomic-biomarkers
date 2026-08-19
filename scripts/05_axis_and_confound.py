#!/usr/bin/env python3
"""Step 05 -- single-axis collapse and the infiltration/purity confound.

Three analyses, all reading the score matrices written by step 03:

1. **PCA of the signature score space.** Per cohort, PCA on the
   standardized samples x signature-score matrix; record PC1 variance
   explained and PC1's correlation with the T-cell-inflamed GEP.

2. **The Venet test.** Regress each signature's mean response-AUROC on
   how strongly it aligns with the shared axis (its mean correlation with
   the T-cell-inflamed GEP across cohorts). The R-squared is the fraction
   of between-signature performance variance explained by one number.

3. **The ESTIMATE confound.** Score ESTIMATE's fixed 141-gene immune and
   stromal signatures, derive a tumor-purity estimate, correlate them with
   the signature axis, and recompute every signature's AUROC after
   linearly partialling the infiltration score out -- the test of whether
   anything specific survives removal of the confound.

Writes ``results/axis_pca.csv``, ``results/venet_alignment.csv``,
``results/estimate_scores.csv`` and ``results/partialled_auroc.csv``.

Usage:  python scripts/05_axis_and_confound.py [--method mean_z]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from icinull import auroc, load_cohort, load_config  # noqa: E402
from icinull.estimate import estimate_scores  # noqa: E402


def passing_cohorts() -> list:
    g = pd.read_csv(ROOT / "results" / "cohort_gate.csv")
    return list(g.loc[g.passes_gate, "cohort"])


def partial_out(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Residuals of y after linear regression on x (both 1-D)."""
    ok = np.isfinite(y) & np.isfinite(x)
    out = np.full_like(y, np.nan, dtype=float)
    if ok.sum() < 3 or np.nanstd(x[ok]) == 0:
        return out
    b = np.polyfit(x[ok], y[ok], 1)
    out[ok] = y[ok] - np.polyval(b, x[ok])
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["mean_z", "ssgsea"], default="mean_z")
    args = ap.parse_args()

    cfg = load_config()
    ref = cfg["axis"]["reference_signature"]
    score_dir = ROOT / "results" / "scores"
    res_dir = ROOT / "results"

    pca_rows, est_rows, part_rows = [], [], []
    align_by_sig: dict = {}
    auroc_by_sig: dict = {}

    for cohort in passing_cohorts():
        sf = score_dir / f"{cohort}_{args.method}.tsv.gz"
        if not sf.exists():
            print(f"  {cohort}: no score matrix for {args.method}, skipped")
            continue

        S = pd.read_csv(sf, sep="\t", index_col=0)
        S = S.loc[:, S.std(ddof=1) > 0].dropna(axis=1, how="any")
        if S.shape[1] < 3 or ref not in S.columns:
            print(f"  {cohort}: score matrix unusable "
                  f"({S.shape[1]} signatures, ref present={ref in S.columns})")
            continue

        c = load_cohort(cohort, config=cfg)
        y = c.labels("primary", cfg)

        # ---- 1. PCA of the signature score space -------------------------
        Z = (S - S.mean()) / S.std(ddof=1)
        U, sv, Vt = np.linalg.svd(Z.to_numpy(), full_matrices=False)
        var_exp = sv ** 2 / (sv ** 2).sum()
        pc1 = pd.Series(U[:, 0] * sv[0], index=S.index)

        gep = S[ref]
        r_pc1_gep = float(np.corrcoef(pc1, gep)[0, 1])
        # Sign of a principal component is arbitrary; orient PC1 so that
        # "more positive" means "more inflamed", then report |r|.
        if r_pc1_gep < 0:
            pc1, r_pc1_gep = -pc1, -r_pc1_gep

        # ---- 3. ESTIMATE infiltration / purity ---------------------------
        est = estimate_scores(c.expr.loc[c.universe])
        common = [s for s in S.index if s in est.index]
        immune = est.loc[common, "immune_score"]
        stromal = est.loc[common, "stromal_score"]
        purity = est.loc[common, "tumor_purity"]

        r_axis_immune = float(np.corrcoef(pc1.loc[common], immune)[0, 1])
        r_axis_purity = float(np.corrcoef(pc1.loc[common], purity)[0, 1])

        pca_rows.append({
            "cohort": cohort,
            "method": args.method,
            "n_signatures": int(S.shape[1]),
            "n_samples": int(S.shape[0]),
            "pc1_var_explained": float(var_exp[0]),
            "pc2_var_explained": float(var_exp[1]) if var_exp.size > 1 else np.nan,
            "pc1_corr_gep": r_pc1_gep,
            "pc1_corr_estimate_immune": r_axis_immune,
            "pc1_corr_tumor_purity": r_axis_purity,
        })

        est_out = est.copy()
        est_out.insert(0, "cohort", cohort)
        est_out.insert(1, "sample_id", est_out.index)
        est_rows.append(est_out)

        # ---- 2 + 3. per-signature alignment, AUROC, partialled AUROC -----
        idx = [s for s in y.index if s in S.index]
        yv = y.loc[idx].to_numpy()
        inf_v = immune.reindex(idx).to_numpy()

        for sig in S.columns:
            sv_all = S[sig]
            r_align = float(np.corrcoef(sv_all, gep)[0, 1])
            align_by_sig.setdefault(sig, []).append(r_align)

            a_raw = auroc(sv_all.loc[idx].to_numpy(), yv)
            auroc_by_sig.setdefault(sig, []).append(a_raw)

            resid = partial_out(sv_all.loc[idx].to_numpy(), inf_v)
            a_part = auroc(resid, yv)

            part_rows.append({
                "cohort": cohort,
                "signature": sig,
                "method": args.method,
                "auroc_raw": a_raw,
                "auroc_partialled": a_part,
                "delta": a_part - a_raw if np.isfinite(a_part) else np.nan,
                "corr_with_gep": r_align,
                "corr_with_infiltration": float(
                    np.corrcoef(sv_all.loc[idx], inf_v)[0, 1]
                ),
            })

        print(f"  {cohort}: PC1 {var_exp[0]:.1%} var, r(PC1,GEP)={r_pc1_gep:.2f}, "
              f"r(PC1,immune)={r_axis_immune:.2f}, "
              f"r(PC1,purity)={r_axis_purity:.2f}")

    if not pca_rows:
        raise SystemExit("No cohorts produced a usable score matrix.")

    pca = pd.DataFrame(pca_rows)
    pca.to_csv(res_dir / f"axis_pca_{args.method}.csv", index=False)
    pd.concat(est_rows).to_csv(res_dir / "estimate_scores.csv", index=False)
    part = pd.DataFrame(part_rows)
    part.to_csv(res_dir / f"partialled_auroc_{args.method}.csv", index=False)

    # ---- the Venet regression ---------------------------------------------
    venet = pd.DataFrame({
        "signature": list(align_by_sig),
        "mean_corr_with_gep": [float(np.nanmean(v)) for v in align_by_sig.values()],
        "mean_auroc": [float(np.nanmean(auroc_by_sig[s])) for s in align_by_sig],
    }).dropna()

    from scipy.stats import pearsonr

    # The reference signature correlates 1.0 with itself by construction, so
    # it contributes a fixed anchor point at the extreme of the x-range and
    # inflates the fit slightly. Report the regression both with and without
    # it; the headline is the self-excluded version, since "alignment to the
    # axis predicts performance" should not be partly carried by the axis
    # being perfectly aligned with itself.
    venet["is_reference"] = venet.signature == ref
    r, p = pearsonr(venet.mean_corr_with_gep, venet.mean_auroc)
    ex = venet[~venet.is_reference]
    r_ex, p_ex = pearsonr(ex.mean_corr_with_gep, ex.mean_auroc)
    venet.attrs["r"] = r_ex
    venet.to_csv(res_dir / f"venet_alignment_{args.method}.csv", index=False)

    print(f"\nPC1 variance explained : mean {pca.pc1_var_explained.mean():.1%}, "
          f"max {pca.pc1_var_explained.max():.1%}")
    print(f"PC1 vs T-cell-inflamed : mean r = {pca.pc1_corr_gep.mean():.2f}")
    print(f"PC1 vs ESTIMATE immune : mean r = "
          f"{pca.pc1_corr_estimate_immune.mean():.2f}")
    print(f"PC1 vs tumor purity    : mean r = "
          f"{pca.pc1_corr_tumor_purity.mean():.2f}")
    print(f"\nVenet test, reference excluded (n={len(ex)}): "
          f"r = {r_ex:.3f}, R^2 = {r_ex ** 2:.3f}, p = {p_ex:.2e}   [headline]")
    print(f"Venet test, reference included (n={len(venet)}): "
          f"r = {r:.3f}, R^2 = {r ** 2:.3f}, p = {p:.2e}")

    pm = part.groupby("signature").agg(
        raw=("auroc_raw", "mean"), partialled=("auroc_partialled", "mean")
    ).dropna()
    print(f"mean AUROC raw {pm.raw.mean():.3f} -> "
          f"partialled {pm.partialled.mean():.3f} "
          f"(infiltration removed)")
    print(f"\n-> {res_dir}/axis_pca_{args.method}.csv")
    print(f"-> {res_dir}/venet_alignment_{args.method}.csv")
    print(f"-> {res_dir}/estimate_scores.csv")
    print(f"-> {res_dir}/partialled_auroc_{args.method}.csv")


if __name__ == "__main__":
    main()
