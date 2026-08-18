#!/usr/bin/env python3
"""Step 09 -- regenerate the manuscript figures from the result tables.

Each function draws one manuscript figure and reads only committed result
CSVs, so a figure can never disagree with the tables. Figures whose
upstream table is absent are skipped with a message rather than drawn from
partial data.

  Figure 1  cohort composition, responder rate, survival availability
  Figure 2  per-signature AUROC heatmap across cohorts
  Figure 3  published signatures vs the size-matched random null
  Figure 4  single-axis collapse (PC1 variance, PC1-GEP, Venet regression)
  Figure 5  the infiltration/purity confound
  Figure 6  forest plots            [drawn by R/07_meta_analysis.R --forest]
  Figure 7  the transferable ceiling and in-cohort optimism
  Figure 8  robustness (scoring method, cancer type, ssGSEA collapse)

Usage:
  python scripts/09_figures.py            # all available
  python scripts/09_figures.py --only 3,5
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
sys.path.insert(0, str(ROOT / "src"))

mpl.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlelocation": "left",
})

BLUE, ORANGE, GREY = "#1f77b4", "#d95f02", "#7f7f7f"


def read(name: str):
    p = RES / name
    if not p.exists():
        return None
    return pd.read_csv(p)


def save(fig, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / name
    fig.savefig(out)
    plt.close(fig)
    print(f"  -> figures/{name}")


# --------------------------------------------------------------------------

def figure1() -> None:
    g = read("cohort_gate.csv")
    if g is None:
        print("fig1: cohort_gate.csv missing, skipped")
        return
    p = g[g.passes_gate].sort_values("n_response_labeled")

    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2),
                             gridspec_kw={"width_ratios": [2, 1.2, 1.2]})
    y = np.arange(len(p))

    ax = axes[0]
    ax.barh(y, p.n_responder, color=BLUE, label="Responder (CR/PR)")
    ax.barh(y, p.n_nonresponder, left=p.n_responder, color="#d9d9d9",
            edgecolor="none", label="Non-responder (SD/PD)")
    for i, (r, n) in enumerate(zip(p.n_responder, p.n_response_labeled)):
        ax.text(n + max(p.n_response_labeled) * 0.015, i, str(int(n)),
                va="center", fontsize=6, color=GREY)
    ax.set_yticks(y, p.cohort)
    ax.set_xlabel("Patients with RECIST response label")
    ax.set_title("Response composition")
    ax.legend(frameon=False, loc="lower right")
    ax.margins(x=0.08)

    # Primary tumour type per cohort, not per-sample biopsy site (see the
    # note in R/01_harmonize.R): Mariathasan is urothelial carcinoma whose
    # clinical column records where each biopsy was taken.
    ctcol = ("primary_cancer_type"
             if "primary_cancer_type" in p and p.primary_cancer_type.notna().all()
             else "cancer_type")
    types = sorted({t for v in p[ctcol].astype(str)
                    for t in v.split("|") if t and t != "nan"})
    cmap = dict(zip(types, plt.cm.tab10.colors))
    ax = axes[1]
    for i, (_, r) in enumerate(p.iterrows()):
        t = str(r[ctcol]).split("|")[0]
        ax.scatter(r.responder_rate, i, s=42, color=cmap.get(t, GREY),
                   edgecolor="white", linewidth=0.6, zorder=3)
    pooled = p.n_responder.sum() / p.n_response_labeled.sum()
    ax.axvline(pooled, ls="--", color=GREY, lw=0.9)
    ax.text(pooled, len(p) - 0.4, f" pooled {pooled:.2f}", fontsize=6,
            color=GREY, va="top")
    ax.set_yticks(y, [""] * len(p))
    ax.set_xlabel("Responder rate")
    ax.set_title("Response rate")
    ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=cmap[t],
                                  label=t, markersize=5) for t in types],
              frameon=False, title="Cancer type", loc="lower right",
              title_fontsize=6)
    ax.margins(x=0.15)

    ax = axes[2]
    for i, (_, r) in enumerate(p.iterrows()):
        if r.n_with_os > 0:
            ax.scatter(r.n_with_os, i, s=42, color=BLUE, edgecolor="white",
                       linewidth=0.6, zorder=3)
            ax.text(r.n_with_os * 1.12, i, str(int(r.n_with_os)),
                    fontsize=6, va="center", color=GREY)
        else:
            ax.text(0.04, i, "no OS", fontsize=6, va="center",
                    color=GREY, style="italic",
                    transform=ax.get_yaxis_transform())
    ax.set_xscale("log")
    ax.set_yticks(y, [""] * len(p))
    ax.set_xlabel("Patients with overall survival (log)")
    ax.set_title("Survival data")

    fig.suptitle(
        f"{len(p)} uniformly-processed public ICI cohorts: "
        f"{int(p.n_response_labeled.sum())} response-labeled patients, "
        f"{int(p.n_with_os.sum())} with OS",
        x=0.01, ha="left", fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "figure1_cohorts.png")


def figure2() -> None:
    perf = read("perf_per_sig_cohort.csv")
    if perf is None:
        print("fig2: perf_per_sig_cohort.csv missing, skipped")
        return
    d = perf[(perf.method == "mean_z") & perf.auroc.notna()]
    piv = d.pivot_table(index="signature", columns="cohort", values="auroc")
    order = piv.mean(axis=1).sort_values(ascending=False)
    n_show = min(40, len(order))
    piv = piv.loc[order.index[:n_show]]
    cols = (d.groupby("cohort").n_pos.max() + d.groupby("cohort").n_neg.max())
    piv = piv[cols.sort_values(ascending=False).index]

    fig, axes = plt.subplots(
        1, 2, figsize=(10, max(4.5, 0.16 * n_show)),
        gridspec_kw={"width_ratios": [3, 1], "wspace": 0.03},
    )
    ax = axes[0]
    im = ax.imshow(piv.to_numpy(), aspect="auto", cmap="RdBu_r",
                   vmin=0.35, vmax=0.65)
    ax.set_xticks(range(piv.shape[1]), piv.columns, rotation=45, ha="right")
    ax.set_yticks(range(piv.shape[0]),
                  [s.replace("_", " ")[:38] for s in piv.index])
    ax.set_title(f"Per-signature response AUROC across {piv.shape[1]} "
                 f"cohorts (top {n_show}, mean-z)")
    cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02)
    cb.set_label("AUROC (0.50 = chance)")

    ax = axes[1]
    m = order.loc[piv.index]
    ax.barh(range(len(m)), m, color=BLUE, height=0.75)
    ax.axvline(0.5, color="black", lw=0.8, ls="--")
    ax.set_yticks(range(len(m)), [""] * len(m))
    ax.set_ylim(len(m) - 0.5, -0.5)
    ax.set_xlim(0.40, max(0.68, float(m.max()) + 0.02))
    ax.set_xlabel("Mean AUROC")
    ax.set_title("Mean across cohorts")
    save(fig, "figure2_signature_auroc.png")


def figure3() -> None:
    null = read("null_results.csv")
    pool = read("null_pooled.csv")
    if null is None or pool is None:
        print("fig3: null_results.csv / null_pooled.csv missing, skipped")
        return

    fig = plt.figure(figsize=(11, 5.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.5, 1], hspace=0.45,
                          wspace=0.28)

    ax = fig.add_subplot(gs[:, 0])
    top = pool.sort_values("mean_observed_auroc", ascending=False).head(32)
    top = top.iloc[::-1]
    y = np.arange(len(top))
    colors = [BLUE if b else "#bdbdbd" for b in top.beats_null_fdr05]
    ax.barh(y, top.mean_observed_auroc, color=colors, height=0.72)
    ax.scatter(top.mean_null_auroc, y, s=13, color=ORANGE, zorder=3,
               label="random null (mean)")
    ax.axvline(0.5, color="black", lw=0.8)
    ax.set_yticks(y, [s.replace("_", " ")[:40] for s in top.signature])
    ax.set_xlim(0.42, float(top.mean_observed_auroc.max()) + 0.02)
    ax.set_xlabel("AUROC (mean across cohorts)")
    ax.set_title("A  Published signatures vs size-matched random null")
    ax.legend(frameon=False, loc="lower right")
    ax.text(0.01, -0.055, "blue = beats null at FDR<0.05; grey = does not",
            transform=ax.transAxes, fontsize=6, color=GREY)

    ax = fig.add_subplot(gs[0, 1])
    bins = np.linspace(0.2, 0.9, 40)
    ax.hist(null.null_mean.dropna(), bins=bins, color=ORANGE, alpha=0.65,
            density=True, label="random null")
    ax.hist(null.observed_auroc.dropna(), bins=bins, color=BLUE, alpha=0.55,
            density=True, label="published")
    ax.axvline(null.null_mean.median(), color=ORANGE, ls="--", lw=0.9)
    ax.axvline(null.observed_auroc.median(), color=BLUE, ls="--", lw=0.9)
    frac = null.beats_null_p05.mean()
    ax.set_xlabel("AUROC")
    ax.set_ylabel("density")
    ax.set_title("B  All signature x cohort tests")
    ax.legend(frameon=False)
    ax.text(0.02, 0.95, f"{frac:.0%} of tests beat null (p<0.05)\n"
            f"median published {null.observed_auroc.median():.3f} vs "
            f"null 95th pct {null.null_q95.median():.3f}",
            transform=ax.transAxes, va="top", fontsize=6)

    ax = fig.add_subplot(gs[1, 1])
    bc = null.groupby("cohort").agg(
        frac=("beats_null_p05", "mean"),
        n=("n_pos", "max"),
    )
    bc["n"] = bc.n + null.groupby("cohort").n_neg.max()
    ax.scatter(bc.n, bc.frac, s=34, color="#333333", zorder=3)
    for c, r in bc.iterrows():
        ax.annotate(c, (r.n, r.frac), fontsize=6,
                    textcoords="offset points", xytext=(4, 1))
    ax.axhline(0.05, ls="--", color=ORANGE, lw=0.9)
    ax.text(bc.n.max(), 0.055, "chance (5%)", fontsize=6, color=ORANGE,
            ha="right")
    ax.set_xlabel("Cohort size (response-labeled n)")
    ax.set_ylabel("Fraction beating null")
    ax.set_title("C  Null-beating tracks cohort, not biology")
    ax.margins(0.12)

    fig.suptitle("Most published ICI-response signatures do not outperform "
                 "size-matched random gene sets", x=0.01, ha="left",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, "figure3_random_null.png")


def figure4(method: str = "mean_z") -> None:
    pca = read(f"axis_pca_{method}.csv")
    ven = read(f"venet_alignment_{method}.csv")
    prov = ROOT / "signatures" / "signature_provenance.csv"
    if pca is None or ven is None:
        print("fig4: axis tables missing, skipped")
        return

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))

    ax = axes[0]
    p = pca.sort_values("pc1_var_explained", ascending=False)
    ax.bar(range(len(p)), p.pc1_var_explained * 100, color=BLUE)
    ax.axhline(p.pc1_var_explained.mean() * 100, ls="--", color="black",
               lw=0.9)
    ax.text(len(p) - 0.5, p.pc1_var_explained.mean() * 100 + 1.5,
            f"mean {p.pc1_var_explained.mean():.0%}", fontsize=6, ha="right")
    ax.set_xticks(range(len(p)), p.cohort, rotation=45, ha="right")
    ax.set_ylabel("PC1 variance explained (%)")
    ax.set_title("A  One axis dominates the\nsignature score space")

    ax = axes[1]
    ax.scatter(pca.pc1_corr_gep, pca.pc1_corr_estimate_immune, s=34,
               color="#333333", zorder=3)
    for _, r in pca.iterrows():
        ax.annotate(r.cohort, (r.pc1_corr_gep, r.pc1_corr_estimate_immune),
                    fontsize=6, textcoords="offset points", xytext=(4, 1))
    ax.set_xlabel("corr(PC1, T-cell-inflamed GEP)")
    ax.set_ylabel("corr(PC1, ESTIMATE immune)")
    ax.set_title("B  PC1 is the inflammation axis")
    ax.margins(0.15)

    ax = axes[2]
    v = ven.dropna()
    canon = set()
    if prov.exists():
        pv = pd.read_csv(prov)
        canon = set(pv.loc[pv.collection == "canonical", "signature"])
    is_c = v.signature.isin(canon)
    ax.scatter(v.loc[~is_c, "mean_corr_with_gep"], v.loc[~is_c, "mean_auroc"],
               s=16, color=GREY, alpha=0.7, label="IOBR gene set")
    ax.scatter(v.loc[is_c, "mean_corr_with_gep"], v.loc[is_c, "mean_auroc"],
               s=26, color=BLUE, label="canonical named", zorder=3)
    if len(v) > 2:
        b = np.polyfit(v.mean_corr_with_gep, v.mean_auroc, 1)
        xs = np.linspace(v.mean_corr_with_gep.min(),
                         v.mean_corr_with_gep.max(), 50)
        ax.plot(xs, np.polyval(b, xs), color=ORANGE, lw=1.4)
        r = np.corrcoef(v.mean_corr_with_gep, v.mean_auroc)[0, 1]
        ax.set_title(f"C  Axis alignment explains {r ** 2:.0%} of\n"
                     f"between-signature performance")
    ax.axhline(0.5, ls=":", color="black", lw=0.8)
    ax.set_xlabel("Signature alignment to axis (corr with GEP)")
    ax.set_ylabel("Mean response AUROC")
    ax.legend(frameon=False, loc="upper left")
    ax.margins(0.06)

    fig.suptitle("Published ICI signatures collapse onto a single "
                 "tumor-inflammation axis", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "figure4_single_axis.png")


def figure5() -> None:
    part = read("partialled_auroc_mean_z.csv")
    pca = read("axis_pca_mean_z.csv")
    null = read("null_results.csv")
    est = read("estimate_scores.csv")
    if part is None or pca is None:
        print("fig5: confound tables missing, skipped")
        return

    pm = part.groupby("signature").agg(
        raw=("auroc_raw", "mean"), part=("auroc_partialled", "mean")
    ).dropna()

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))

    ax = axes[0]
    ax.scatter(pm.raw, pm.part, s=18, color=GREY, alpha=0.75)
    lim = [min(pm.raw.min(), pm.part.min()) - 0.02,
           max(pm.raw.max(), pm.part.max()) + 0.02]
    ax.plot(lim, lim, ls="--", color="black", lw=0.8)
    ax.axhline(0.5, color=ORANGE, lw=0.9)
    ax.axvline(0.5, color=ORANGE, lw=0.9)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Signature AUROC (raw)")
    ax.set_ylabel("AUROC after removing infiltration")
    ax.set_title("A  Signal collapses toward chance\nwhen infiltration removed")

    ax = axes[1]
    bars = {"best 10\nsignatures": pm.raw.nlargest(10).mean(),
            "all\nsignatures": pm.raw.mean(),
            "after removing\ninfiltration": pm.part.mean()}
    if null is not None:
        bars["random\nnull"] = null.null_mean.mean()
    cols = [BLUE, "#9ecae1", ORANGE, "#bdbdbd"][:len(bars)]
    ax.bar(range(len(bars)), list(bars.values()), color=cols)
    for i, v in enumerate(bars.values()):
        ax.text(i, v + 0.004, f"{v:.3f}", ha="center", fontsize=6)
    ax.axhline(0.5, color="black", lw=0.9)
    ax.set_xticks(range(len(bars)), list(bars), fontsize=6)
    ax.set_ylim(0.45, max(bars.values()) + 0.03)
    ax.set_ylabel("Mean AUROC")
    ax.set_title("B  What the signatures are worth")

    ax = axes[2]
    ax.scatter(pca.pc1_corr_estimate_immune, pca.pc1_corr_tumor_purity,
               s=34, color="#333333", zorder=3)
    for _, r in pca.iterrows():
        ax.annotate(r.cohort, (r.pc1_corr_estimate_immune,
                               r.pc1_corr_tumor_purity),
                    fontsize=6, textcoords="offset points", xytext=(4, 1))
    ax.set_xlabel("corr(signature axis, immune infiltration)")
    ax.set_ylabel("corr(signature axis, tumor purity)")
    ax.set_title("C  The shared axis is infiltration,\ninversely purity")
    ax.margins(0.15)

    fig.suptitle("A single infiltration/purity score accounts for the "
                 "published signatures", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "figure5_infiltration_confound.png")


def figure7() -> None:
    ceil = read("ceiling_loco.csv")
    opt = read("optimism.csv")
    if ceil is None:
        print("fig7: ceiling_loco.csv missing, skipped")
        return

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8),
                             gridspec_kw={"width_ratios": [1.4, 1]})

    LABEL = {
        "gep_axis": "T-cell-inflamed GEP (1 sig)",
        "estimate_immune": "ESTIMATE immune (1 axis)",
        "full_panel": "Full signature panel",
        "full_panel_plus_infiltration": "Full panel + infiltration",
        "random_panel": "Random gene-set panel",
    }
    c = ceil.dropna(subset=["loco_auroc"]).sort_values("loco_auroc")
    ax = axes[0]
    colors = [ORANGE if p == "estimate_immune"
              else "#bdbdbd" if p == "random_panel" else BLUE
              for p in c.predictor]
    y = np.arange(len(c))
    ax.barh(y, c.loco_auroc, color=colors, height=0.68)
    right = float(c.loco_auroc.max())
    for i, r in enumerate(c.itertuples()):
        # Value labels go at a common right-hand column so a CI whisker
        # can never overprint the number it belongs to.
        if np.isfinite(r.ci_low):
            ax.plot([r.ci_low, r.ci_high], [i, i], color="black", lw=1.2)
            right = max(right, float(r.ci_high))
        ax.text(1.0, i, f"{r.loco_auroc:.3f}", fontsize=6, va="center",
                ha="right", transform=ax.get_yaxis_transform())
    ax.axvline(0.5, color="black", lw=0.9)
    ax.set_yticks(y, [LABEL.get(p, p) for p in c.predictor])
    ax.set_xlim(0.42, right + 0.055)
    ax.set_xlabel("Leave-one-cohort-out held-out AUROC")
    ax.set_title("A  The transferable ceiling")

    ax = axes[1]
    if opt is not None:
        LAB2 = {
            "best_signature_selected_in_cohort": "best signature\n(in-cohort)",
            "within_cohort_cv_panel": "trained panel\n(in-cohort CV)",
            "transferable_ceiling_gep_axis": "1 axis\n(held out)",
            "transferable_full_panel": "full panel\n(held out)",
        }
        o = opt.dropna(subset=["auroc"])
        cols = ["#d62728" if "in_cohort" in q or "within_cohort" in q
                else BLUE for q in o.quantity]
        ax.bar(range(len(o)), o.auroc, color=cols)
        for i, v in enumerate(o.auroc):
            ax.text(i, v + 0.005, f"{v:.3f}", ha="center", fontsize=6)
        ax.axhline(0.5, color="black", lw=0.9)
        ax.set_xticks(range(len(o)), [LAB2.get(q, q) for q in o.quantity],
                      fontsize=6)
        ax.set_ylim(0.45, float(o.auroc.max()) + 0.04)
        ax.set_ylabel("AUROC")
        ax.set_title("B  In-cohort optimism vs what transfers")
        ax.text(0.02, 0.96, "red = what discovery papers report",
                transform=ax.transAxes, fontsize=6, color="#d62728", va="top")

    fig.suptitle("Bulk-RNA ICI-response prediction saturates at one "
                 "inflammation axis", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "figure7_transferable_ceiling.png")


def figure8() -> None:
    perf = read("perf_per_sig_cohort.csv")
    null = read("null_results.csv")
    gate = read("cohort_gate.csv")
    if perf is None:
        print("fig8: perf_per_sig_cohort.csv missing, skipped")
        return

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5))

    ax = axes[0]
    mz = perf[(perf.method == "mean_z") & perf.auroc.notna()]
    ss = perf[(perf.method == "ssgsea") & perf.auroc.notna()]
    j = (mz.groupby("signature").auroc.mean().to_frame("mean_z")
         .join(ss.groupby("signature").auroc.mean().rename("ssgsea")).dropna())
    ax.scatter(j.mean_z, j.ssgsea, s=16, color=GREY, alpha=0.75)
    lim = [min(j.min()) - 0.02, max(j.max()) + 0.02]
    ax.plot(lim, lim, ls="--", color="black", lw=0.8)
    r = j.mean_z.corr(j.ssgsea)
    rho = j.mean_z.corr(j.ssgsea, method="spearman")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("mean-z per-signature AUROC")
    ax.set_ylabel("ssGSEA per-signature AUROC")
    ax.set_title(f"A  Scoring-method invariance\n(r={r:.2f}, rho={rho:.2f})")

    ax = axes[1]
    if null is not None and gate is not None:
        ct = gate.set_index("cohort").cancer_type.str.split("|").str[0]
        n = null.copy()
        n["cancer_type"] = n.cohort.map(ct)
        by = n.groupby("cancer_type").agg(
            frac=("beats_null_p05", "mean"),
            n=("n_pos", "max"),
        )
        nn = (n.groupby("cancer_type").n_pos.max()
              + n.groupby("cancer_type").n_neg.max())
        by = by.sort_values("frac")
        ax.barh(range(len(by)), by.frac, color=BLUE, height=0.7)
        ax.axvline(0.05, ls="--", color=ORANGE, lw=0.9)
        ax.set_yticks(range(len(by)),
                      [f"{t} (n={int(nn[t])})" for t in by.index], fontsize=6)
        ax.set_xlabel("Fraction of signatures beating null")
        ax.set_title("B  Null-beating holds across\ncancer types")

    ax = axes[2]
    vals, labs = [], []
    for m in ("mean_z", "ssgsea"):
        pca = read(f"axis_pca_{m}.csv")
        ven = read(f"venet_alignment_{m}.csv")
        if pca is None or ven is None:
            continue
        v = ven.dropna()
        r2 = (np.corrcoef(v.mean_corr_with_gep, v.mean_auroc)[0, 1] ** 2
              if len(v) > 2 else np.nan)
        vals.append((pca.pc1_var_explained.mean(), r2))
        labs.append(m)
    if vals:
        x = np.arange(2)
        w = 0.36
        for i, (lab, (pc1, r2)) in enumerate(zip(labs, vals)):
            ax.bar(x + (i - 0.5) * w, [pc1, r2], width=w,
                   color=[BLUE, ORANGE][i], label=lab)
            for xi, v in zip(x + (i - 0.5) * w, [pc1, r2]):
                if np.isfinite(v):
                    ax.text(xi, v + 0.015, f"{v:.2f}", ha="center", fontsize=6)
        ax.set_xticks(x, ["PC1 var\nexplained", "Venet R^2\n(axis->perf)"],
                      fontsize=6)
        ax.set_ylim(0, 1.0)
        ax.set_ylabel("Proportion")
        ax.set_title("C  Collapse robust to scoring method")
        ax.legend(frameon=False)

    fig.suptitle("Robustness: conclusions invariant to scoring method and "
                 "cancer type", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "figure8_robustness.png")


FIGURES = {1: figure1, 2: figure2, 3: figure3, 4: figure4,
           5: figure5, 7: figure7, 8: figure8}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="",
                    help="comma-separated figure numbers, e.g. 3,5")
    args = ap.parse_args()

    want = ([int(x) for x in args.only.split(",") if x.strip()]
            if args.only else sorted(FIGURES))
    for n in want:
        if n not in FIGURES:
            print(f"no such figure: {n}")
            continue
        print(f"figure {n}:")
        FIGURES[n]()
    print("\nFigure 6 (forest plots) is drawn by:")
    print("  Rscript R/07_meta_analysis.R --forest "
          "T_cell_inflamed_GEP_Ayers_et_al,IFNG_6gene_Ayers")


if __name__ == "__main__":
    main()
