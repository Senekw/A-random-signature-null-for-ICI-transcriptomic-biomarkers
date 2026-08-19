#!/usr/bin/env python3
"""Step 06 -- the transferable ceiling under leave-one-cohort-out validation.

A null result is only half the story. This step measures what bulk tumor
RNA can actually deliver when evaluated the way a deployable biomarker
would be: train on all cohorts but one, predict the held-out cohort, pool
the held-out predictions, and compute a single AUROC.

Predictors compared:

  gep_axis                       the single T-cell-inflamed GEP score
  estimate_immune                a single ESTIMATE immune score
  full_panel                     L2-logistic over all 102 signatures
  full_panel_plus_infiltration   the panel plus the infiltration score
  random_panel                   a size-matched panel of random gene sets

with a cohort-level (cluster) bootstrap CI on the axis ceiling, and the
in-cohort optimism comparison: the best signature *selected within* each
cohort, and a within-cohort cross-validated panel, versus what transfers.

Writes ``results/ceiling_loco.csv`` and ``results/optimism.csv``.

Usage:  python scripts/06_transferable_ceiling.py [--bootstrap 2000]
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


def passing_cohorts() -> list:
    g = pd.read_csv(ROOT / "results" / "cohort_gate.csv")
    return list(g.loc[g.passes_gate, "cohort"])


def load_panel(method: str = "mean_z") -> tuple:
    """Signature scores + labels + infiltration for all passing cohorts."""
    cfg = load_config()
    est = pd.read_csv(ROOT / "results" / "estimate_scores.csv")
    est = est.set_index(["cohort", "sample_id"])

    blocks, labels, infil = [], [], []
    for cohort in passing_cohorts():
        sf = ROOT / "results" / "scores" / f"{cohort}_{method}.tsv.gz"
        if not sf.exists():
            continue
        S = pd.read_csv(sf, sep="\t", index_col=0)
        c = load_cohort(cohort, config=cfg)
        y = c.labels("primary", cfg)
        idx = [s for s in y.index if s in S.index]
        if len(idx) < 10:
            continue

        S = S.loc[idx]
        S.index = pd.MultiIndex.from_product([[cohort], idx],
                                             names=["cohort", "sample_id"])
        blocks.append(S)
        labels.append(pd.Series(y.loc[idx].to_numpy(), index=S.index))

        if (cohort, idx[0]) in est.index:
            infil.append(est.loc[cohort].loc[idx, "immune_score"]
                         .set_axis(S.index))
        else:
            infil.append(pd.Series(np.nan, index=S.index))

    if not blocks:
        raise SystemExit("No score matrices found -- run step 03 first.")

    # The panel needs a fixed feature set across folds, so it is the
    # intersection of signatures scored in every cohort. A cohort whose
    # assay covers very few signature genes (a targeted panel rather than
    # whole transcriptome) would shrink that intersection toward zero and
    # silently redefine the predictor, so such cohorts are reported and
    # dropped from the panel analysis -- they still count in the
    # per-signature and null analyses, where each test stands alone.
    sizes = {b.index.get_level_values("cohort")[0]: b.shape[1] for b in blocks}
    med = float(np.median(list(sizes.values())))
    keep = [b for b in blocks
            if b.shape[1] >= max(0.5 * med, 10)]
    dropped = {c: n for c, n in sizes.items()
               if n < max(0.5 * med, 10)}
    if dropped:
        print("dropped from panel analysis (too few signatures scored -- "
              "sparse or targeted assay):")
        for c, n in dropped.items():
            print(f"  {c}: {n} signatures scored (median across cohorts "
                  f"{med:.0f})")
    if not keep:
        raise SystemExit("No cohort has enough scored signatures for the "
                         "panel analysis.")
    blocks = keep

    common = set(blocks[0].columns)
    for b in blocks[1:]:
        common &= set(b.columns)
    cols = sorted(common)
    if not cols:
        raise SystemExit("No signature is scored in every retained cohort.")
    print(f"panel features: {len(cols)} signatures common to "
          f"{len(blocks)} cohorts")

    keep_idx = pd.concat(blocks).index
    labels = [s.loc[s.index.isin(keep_idx)] for s in labels]
    infil = [s.loc[s.index.isin(keep_idx)] for s in infil]

    X = pd.concat([b[cols] for b in blocks])
    y = pd.concat(labels)
    inf = pd.concat(infil)
    return X, y, inf, cfg


def loco_auroc(X: pd.DataFrame, y: pd.Series, C: float, seed: int = 0) -> tuple:
    """Pooled held-out AUROC under leave-one-cohort-out validation."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    cohorts = list(dict.fromkeys(X.index.get_level_values("cohort")))
    preds, truth, keys = [], [], []

    for held in cohorts:
        tr = X.index.get_level_values("cohort") != held
        te = ~tr
        Xtr, Xte = X[tr], X[te]
        ytr, yte = y[tr], y[te]
        if ytr.nunique() < 2 or yte.nunique() < 2:
            continue

        sc = StandardScaler().fit(Xtr)
        m = LogisticRegression(C=C, max_iter=5000)
        m.fit(sc.transform(Xtr), ytr)
        p = m.predict_proba(sc.transform(Xte))[:, 1]

        preds.append(p)
        truth.append(yte.to_numpy())
        keys.extend([held] * len(p))

    if not preds:
        return float("nan"), pd.DataFrame()

    p = np.concatenate(preds)
    t = np.concatenate(truth)
    pooled = pd.DataFrame({"cohort": keys, "pred": p, "y": t})
    return auroc(p, t), pooled


def single_feature_loco(X: pd.DataFrame, y: pd.Series, col: str,
                        standardize: str = "within") -> tuple:
    """A single score used directly as the predictor (no model is fitted).

    A single score needs no training, so leave-one-cohort-out is degenerate
    for it and the only real choice is how to put cohorts on a common scale
    before pooling predictions. Raw scores cannot be pooled: their location
    and spread differ per cohort, so pooling them would measure cohort
    offsets rather than discrimination.

    ``standardize``:

    * ``"within"``  -- z-score inside each cohort. Rank-preserving within a
      cohort, so every per-cohort AUROC is untouched; but it does use the
      held-out cohort's own mean and SD, which is information a deployed
      single-gene test would not have for a lone patient.
    * ``"train"``   -- z-score the held-out cohort using the mean and SD
      pooled over the other cohorts, matching the standardization rule the
      Methods state for the trained panel.

    Both are reported. They differ only in how cohorts are aligned, so a
    large gap between them is itself a finding about cross-cohort
    comparability rather than about the signature.
    """
    if standardize not in ("within", "train"):
        raise ValueError(f"standardize must be 'within' or 'train', "
                         f"got {standardize!r}")

    parts, truth, keys = [], [], []
    for cohort, grp in X.groupby(level="cohort"):
        s = grp[col]
        if standardize == "within":
            mu, sd = s.mean(), s.std(ddof=1)
        else:
            other = X.loc[X.index.get_level_values("cohort") != cohort, col]
            mu, sd = other.mean(), other.std(ddof=1)
        if not np.isfinite(sd) or sd == 0:
            continue
        parts.append(((s - mu) / sd).to_numpy())
        truth.append(y.loc[grp.index].to_numpy())
        keys.extend([cohort] * len(grp))

    if not parts:
        return float("nan"), pd.DataFrame(columns=["cohort", "pred", "y"])

    p = np.concatenate(parts)
    t = np.concatenate(truth)
    pooled = pd.DataFrame({"cohort": keys, "pred": p, "y": t})
    return auroc(p, t), pooled


def cluster_bootstrap(pooled: pd.DataFrame, n: int, seed: int) -> tuple:
    """Cohort-level bootstrap CI of a pooled AUROC."""
    rng = np.random.default_rng(seed)
    cohorts = pooled.cohort.unique()
    draws = []
    for _ in range(n):
        pick = rng.choice(cohorts, size=len(cohorts), replace=True)
        sub = pd.concat([pooled[pooled.cohort == c] for c in pick])
        a = auroc(sub.pred.to_numpy(), sub.y.to_numpy())
        if np.isfinite(a):
            draws.append(a)
    if not draws:
        return float("nan"), float("nan")
    return (float(np.quantile(draws, 0.025)),
            float(np.quantile(draws, 0.975)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["mean_z", "ssgsea"], default="mean_z")
    ap.add_argument("--bootstrap", type=int, default=None)
    args = ap.parse_args()

    X, y, inf, cfg = load_panel(args.method)
    C = cfg["ceiling"]["C"]
    seed = cfg["seed"]
    n_boot = args.bootstrap or cfg["ceiling"]["cluster_bootstrap_n"]
    ref = cfg["axis"]["reference_signature"]

    n_cohorts = X.index.get_level_values("cohort").nunique()
    print(f"panel: {X.shape[1]} signatures x {X.shape[0]} labeled samples "
          f"across {n_cohorts} cohorts")
    if n_cohorts < 3:
        print(f"WARNING: only {n_cohorts} cohorts in the panel. "
              f"Leave-one-cohort-out still runs, but each fold trains on "
              f"{n_cohorts - 1} cohort(s), so the trained-panel rows "
              f"(full_panel, full_panel_plus_infiltration, random_panel) are "
              f"dominated by fold noise: they can fall below chance and can "
              f"order arbitrarily with respect to each other. Do NOT read a "
              f"directional result off them. Only the full compendium "
              f"(>=3, published: 15 cohorts) supports the ceiling comparison.")

    rows = []

    a_gep, pooled_gep = single_feature_loco(X, y, ref)
    lo, hi = cluster_bootstrap(pooled_gep, n_boot, seed)
    rows.append({"predictor": "gep_axis", "n_features": 1,
                 "loco_auroc": a_gep, "ci_low": lo, "ci_high": hi,
                 "standardization": "within_cohort"})

    # Same axis, standardized on the training cohorts instead -- the rule the
    # Methods state for the panel. Reported alongside rather than instead, so
    # the ceiling does not silently depend on which convention was chosen.
    a_gep_tr, pooled_gep_tr = single_feature_loco(X, y, ref,
                                                 standardize="train")
    lo_tr, hi_tr = cluster_bootstrap(pooled_gep_tr, n_boot, seed)
    rows.append({"predictor": "gep_axis_train_standardized", "n_features": 1,
                 "loco_auroc": a_gep_tr, "ci_low": lo_tr, "ci_high": hi_tr,
                 "standardization": "training_cohorts"})

    Xi = X.copy()
    Xi["__infiltration"] = inf.reindex(X.index)
    if Xi["__infiltration"].notna().all():
        a_est, pooled_est = single_feature_loco(Xi, y, "__infiltration")
        rows.append({"predictor": "estimate_immune", "n_features": 1,
                     "loco_auroc": a_est, "ci_low": np.nan, "ci_high": np.nan,
                     "standardization": "within_cohort"})

    a_panel, _ = loco_auroc(X, y, C)
    rows.append({"predictor": "full_panel", "n_features": X.shape[1],
                 "loco_auroc": a_panel, "ci_low": np.nan, "ci_high": np.nan})

    if Xi["__infiltration"].notna().all():
        a_pi, _ = loco_auroc(Xi, y, C)
        rows.append({"predictor": "full_panel_plus_infiltration",
                     "n_features": Xi.shape[1], "loco_auroc": a_pi,
                     "ci_low": np.nan, "ci_high": np.nan})

    # random panel: same number of features, drawn as random gene sets
    # The Methods specify "a size-matched random 102-gene-set panel", i.e.
    # matched to the LIBRARY size, not to however many signatures survived
    # the cross-cohort intersection. Using X.shape[1] silently shrinks the
    # random panel (98 here) and makes it a weaker comparator than the one
    # the paper describes.
    n_library = len(pd.read_csv(ROOT / "signatures" / "signature_provenance.csv"))
    rand = _random_panel(n_library, args.method, seed, cfg)
    if rand is not None:
        a_rand, _ = loco_auroc(rand.loc[X.index], y, C)
        rows.append({"predictor": "random_panel", "n_features": rand.shape[1],
                     "loco_auroc": a_rand, "ci_low": np.nan, "ci_high": np.nan})

    ceiling = pd.DataFrame(rows).sort_values("loco_auroc", ascending=False)
    ceiling.to_csv(ROOT / "results" / "ceiling_loco.csv", index=False)

    # ---- in-cohort optimism ------------------------------------------------
    best_in, panel_in = [], []
    for cohort, grp in X.groupby(level="cohort"):
        yv = y.loc[grp.index]
        if yv.nunique() < 2:
            continue
        aur = [auroc(grp[c].to_numpy(), yv.to_numpy()) for c in grp.columns]
        aur = [a for a in aur if np.isfinite(a)]
        if aur:
            best_in.append(max(max(aur), 1 - min(aur)))
        panel_in.append(_within_cohort_cv(grp, yv, C, seed))

    opt = pd.DataFrame([
        {"quantity": "best_signature_selected_in_cohort",
         "auroc": float(np.nanmean(best_in))},
        {"quantity": "within_cohort_cv_panel",
         "auroc": float(np.nanmean([v for v in panel_in if np.isfinite(v)]))},
        {"quantity": "transferable_ceiling_gep_axis", "auroc": a_gep},
        {"quantity": "transferable_full_panel", "auroc": a_panel},
    ])
    opt["optimism_gap_vs_transferable"] = opt.auroc - a_gep
    opt.to_csv(ROOT / "results" / "optimism.csv", index=False)

    print("\nleave-one-cohort-out held-out AUROC:")
    for _, r in ceiling.iterrows():
        ci = (f"  [{r.ci_low:.3f}, {r.ci_high:.3f}]"
              if np.isfinite(r.ci_low) else "")
        print(f"  {r.predictor:<30} {r.loco_auroc:.3f}{ci}")

    print("\nin-cohort optimism:")
    for _, r in opt.iterrows():
        print(f"  {r.quantity:<36} {r.auroc:.3f}")
    print(f"\n-> {ROOT / 'results' / 'ceiling_loco.csv'}")
    print(f"-> {ROOT / 'results' / 'optimism.csv'}")


def _within_cohort_cv(X: pd.DataFrame, y: pd.Series, C: float, seed: int) -> float:
    """Stratified 5-fold CV AUROC of the trained panel inside one cohort."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler

    if y.nunique() < 2 or int(y.sum()) < 3 or int((y == 0).sum()) < 3:
        return float("nan")

    k = min(5, int(y.sum()), int((y == 0).sum()))
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    p = np.full(len(y), np.nan)
    Xv, yv = X.to_numpy(), y.to_numpy()

    for tr, te in skf.split(Xv, yv):
        sc = StandardScaler().fit(Xv[tr])
        m = LogisticRegression(C=C, max_iter=5000)
        m.fit(sc.transform(Xv[tr]), yv[tr])
        p[te] = m.predict_proba(sc.transform(Xv[te]))[:, 1]

    return auroc(p, yv)


def _random_panel(n_sets: int, method: str, seed: int, cfg: dict):
    """Score n_sets random gene sets per cohort, as a panel baseline."""
    from icinull.scoring import mean_z_score

    sig_sizes = pd.read_csv(ROOT / "signatures" / "signature_provenance.csv")
    sizes = sig_sizes.n_genes.to_numpy()

    cohorts = passing_cohorts()
    loaded = [(name, load_cohort(name, config=cfg)) for name in cohorts]

    # The panel's gene sets must be THE SAME in every cohort. Drawing a fresh
    # set per cohort would make column `random_007` a different gene set in
    # each one, so a model trained on the other cohorts would be applied to
    # features that have nothing to do with the ones it was fitted on -- that
    # measures nothing, and it drives the held-out AUROC toward 0.5 for a
    # reason unrelated to the comparison being made. Sets are therefore drawn
    # once, from the intersection of the cohorts' expressed universes, so
    # every column means the same thing across folds.
    shared_pool = set(loaded[0][1].universe)
    for _, c in loaded[1:]:
        shared_pool &= set(c.universe)
    pool = np.asarray(sorted(shared_pool), dtype=object)
    if len(pool) < 100:
        print(f"note: only {len(pool)} genes are expressed in every cohort; "
              f"the random panel is drawn from that shared universe.")

    rng = np.random.default_rng(seed)
    draws = []
    for j in range(n_sets):
        size = int(min(sizes[j % len(sizes)], len(pool) - 1))
        draws.append(rng.choice(pool, size=max(size, 2), replace=False))

    blocks = []
    for name, c in loaded:
        y = c.labels("primary", cfg)
        idx = [s for s in y.index if s in c.expr.columns]
        if len(idx) < 10:
            continue
        cols = {f"random_{j:03d}": mean_z_score(c.expr, g, 1).loc[idx]
                for j, g in enumerate(draws)}
        B = pd.DataFrame(cols)
        B.index = pd.MultiIndex.from_product([[name], idx],
                                             names=["cohort", "sample_id"])
        blocks.append(B)

    return pd.concat(blocks) if blocks else None


if __name__ == "__main__":
    main()
