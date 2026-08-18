#!/usr/bin/env Rscript
# Step 07 -- random-effects meta-analysis of per-signature AUROC.
#
# Pools each signature across cohorts on the logit-AUROC scale using
# metafor's REML estimator with the Hartung-Knapp-Sidik-Jonkman variance
# adjustment (conservative for the small number of cohorts here), and
# reports I^2, tau^2, the prediction interval, and leave-one-cohort-out
# diagnostics.
#
# Input : results/perf_per_sig_cohort.csv  (step 03)
# Output: results/meta_results.csv         one row per signature
#         results/meta_loco.csv            leave-one-cohort-out estimates
#         figures/forest_<signature>.png   forest plots (--forest)
#
# Usage:
#   Rscript R/07_meta_analysis.R
#   Rscript R/07_meta_analysis.R --forest T_cell_inflamed_GEP_Ayers_et_al,IFNG_6gene_Ayers

suppressPackageStartupMessages({
  library(metafor)
})

ROOT <- "."
RES  <- file.path(ROOT, "results")
FIG  <- file.path(ROOT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
forest_arg <- NULL
if (any(grepl("^--forest", args))) {
  i <- grep("^--forest", args)[1]
  forest_arg <- if (grepl("=", args[i])) sub("^--forest=", "", args[i])
                else args[i + 1]
}

perf_file <- file.path(RES, "perf_per_sig_cohort.csv")
if (!file.exists(perf_file))
  stop("results/perf_per_sig_cohort.csv missing -- run scripts/03_score_signatures.py")

perf <- read.csv(perf_file, stringsAsFactors = FALSE)
perf <- perf[perf$method == "mean_z" &
               is.finite(perf$logit_auroc) &
               is.finite(perf$logit_auroc_se) &
               perf$logit_auroc_se > 0, ]

if (!nrow(perf)) stop("No evaluable signature x cohort rows to meta-analyze.")

sigs <- sort(unique(perf$signature))
message("Meta-analyzing ", length(sigs), " signatures across ",
        length(unique(perf$cohort)), " cohorts ...")

inv_logit <- function(x) 1 / (1 + exp(-x))

rows <- list()
loco_rows <- list()

for (s in sigs) {
  d <- perf[perf$signature == s, ]
  # A single cohort cannot support a random-effects fit; report the
  # cohort estimate as-is rather than a spurious pooled one.
  if (nrow(d) < 2) {
    rows[[length(rows) + 1]] <- data.frame(
      signature = s, collection = d$collection[1], k = nrow(d),
      pooled_auroc = inv_logit(d$logit_auroc[1]),
      ci_low = NA_real_, ci_high = NA_real_,
      pi_low = NA_real_, pi_high = NA_real_,
      tau2 = NA_real_, I2 = NA_real_, QEp = NA_real_,
      p_vs_half = NA_real_, ci_excludes_half = NA,
      loco_min = NA_real_, loco_max = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }

  fit <- try(rma(yi = d$logit_auroc, sei = d$logit_auroc_se,
                 method = "REML", test = "knha"), silent = TRUE)
  if (inherits(fit, "try-error")) next

  pred <- predict(fit)

  # leave-one-cohort-out
  lo <- rep(NA_real_, nrow(d))
  for (j in seq_len(nrow(d))) {
    dj <- d[-j, ]
    if (nrow(dj) < 2) next
    fj <- try(rma(yi = dj$logit_auroc, sei = dj$logit_auroc_se,
                  method = "REML", test = "knha"), silent = TRUE)
    if (!inherits(fj, "try-error")) lo[j] <- inv_logit(as.numeric(fj$b))
  }

  loco_rows[[length(loco_rows) + 1]] <- data.frame(
    signature = s, cohort_left_out = d$cohort,
    pooled_auroc_without = lo, stringsAsFactors = FALSE
  )

  rows[[length(rows) + 1]] <- data.frame(
    signature   = s,
    collection  = d$collection[1],
    k           = nrow(d),
    pooled_auroc = inv_logit(as.numeric(fit$b)),
    ci_low      = inv_logit(fit$ci.lb),
    ci_high     = inv_logit(fit$ci.ub),
    pi_low      = inv_logit(pred$pi.lb),
    pi_high     = inv_logit(pred$pi.ub),
    tau2        = fit$tau2,
    I2          = fit$I2,
    QEp         = fit$QEp,
    p_vs_half   = fit$pval,
    ci_excludes_half = inv_logit(fit$ci.lb) > 0.5,
    loco_min    = if (all(is.na(lo))) NA_real_ else min(lo, na.rm = TRUE),
    loco_max    = if (all(is.na(lo))) NA_real_ else max(lo, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

meta <- do.call(rbind, rows)
meta <- meta[order(-meta$pooled_auroc), ]
write.csv(meta, file.path(RES, "meta_results.csv"), row.names = FALSE)

if (length(loco_rows))
  write.csv(do.call(rbind, loco_rows), file.path(RES, "meta_loco.csv"),
            row.names = FALSE)

# ---- forest plots ---------------------------------------------------------
if (!is.null(forest_arg)) {
  want <- trimws(strsplit(forest_arg, ",")[[1]])
  for (s in want) {
    d <- perf[perf$signature == s, ]
    if (nrow(d) < 2) { message("skip forest (k<2): ", s); next }
    d <- d[order(d$logit_auroc), ]
    fit <- rma(yi = d$logit_auroc, sei = d$logit_auroc_se,
               method = "REML", test = "knha")

    png(file.path(FIG, paste0("forest_", s, ".png")),
        width = 1500, height = 1100, res = 170)
    # `forest.rma` draws the summary polygon itself; label it via mlab
    # rather than adding a second one with addpoly.
    forest(fit, transf = inv_logit, refline = 0.5, slab = d$cohort,
           xlab = "Response AUROC",
           header = c(gsub("_", " ", s), "AUROC [95% CI]"),
           at = c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
           mlab = sprintf("Pooled (REML + Hartung-Knapp): %.3f",
                          inv_logit(as.numeric(fit$b))))
    dev.off()
    message("  forest -> figures/forest_", s, ".png")
  }
}

# ---- summary --------------------------------------------------------------
ok <- meta[!is.na(meta$ci_low), ]
message("\nsignatures meta-analyzed        : ", nrow(ok))
if (nrow(ok)) {
  b <- ok[1, ]
  message(sprintf("best pooled signature          : %s  AUROC %.3f (95%% CI %.3f-%.3f)",
                  b$signature, b$pooled_auroc, b$ci_low, b$ci_high))
  message("pooled CI excludes 0.5         : ",
          sum(ok$ci_excludes_half, na.rm = TRUE), " / ", nrow(ok))
  message(sprintf("median I^2                      : %.1f%%",
                  median(ok$I2, na.rm = TRUE)))
  message(sprintf("median tau^2                    : %.4f",
                  median(ok$tau2, na.rm = TRUE)))
}
message("-> ", file.path(RES, "meta_results.csv"))
