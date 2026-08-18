#!/usr/bin/env Rscript
# Step 01 -- harmonize each cohort MultiAssayExperiment to flat tables.
#
# For every downloaded cohort this writes two files into data/harmonized/:
#
#   <cohort>_expr.tsv.gz   symbol x sample log2-TPM matrix
#   <cohort>_clin.tsv      one row per sample, harmonized clinical fields
#
# and one summary table, results/cohort_inventory.csv.
#
# Decisions applied here, all pre-specified (config/analysis_config.yaml):
#
#   * assay              expr_gene_tpm (log2-TPM, GENCODE gene models)
#   * gene identity      rowData$gene_name (HUGO symbol); rows with no
#                        symbol are dropped
#   * duplicate symbols  collapsed to the row with the highest mean
#                        expression across samples
#   * response           the compendium's own RECIST-derived R/NR label
#                        (R = CR/PR, NR = SD/PD), carried through as-is;
#                        raw `recist` is also retained so the strict
#                        CR/PR-vs-PD sensitivity analysis can be run
#                        downstream without re-reading the MAE objects
#
# Everything downstream reads only these flat tables, so the heavy
# Bioconductor dependency is confined to this one step.
#
# Usage:  Rscript R/01_harmonize.R

suppressPackageStartupMessages({
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})

source(file.path(dirname(sub("^--file=", "",
  c(grep("^--file=", commandArgs(FALSE), value = TRUE), "./R/x")[1])),
  "_root.R"))
ROOT <- repo_root()
RAW  <- file.path(ROOT, "data", "raw")
OUT  <- file.path(ROOT, "data", "harmonized")
RES  <- file.path(ROOT, "results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

ASSAY_PREF <- c("expr_gene_tpm", "expr", "expr_gene_counts")

CLIN_FIELDS <- c(
  "patientid", "sex", "age", "cancer_type", "histo", "tissueid",
  "treatmentid", "treatment", "stage", "recist", "response",
  "response.other.info", "survival_time_os", "event_occurred_os",
  "survival_time_pfs", "event_occurred_pfs", "survival_unit",
  "TMB_raw", "TMB_perMb", "data_set"
)

files <- sort(list.files(RAW, pattern = "\\.rds$", full.names = TRUE))
if (!length(files)) stop("No .rds cohorts in ", RAW,
                         " -- run R/00_download_cohorts.R first.")

inv <- list()

for (f in files) {
  cohort <- sub("^ICB_", "", sub("\\.rds$", "", basename(f)))
  message("== ", cohort)

  mae <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(mae)) { message("   unreadable, skipped"); next }

  exps  <- names(experiments(mae))
  which <- ASSAY_PREF[ASSAY_PREF %in% exps]
  if (!length(which)) {
    message("   no expression assay (has: ", paste(exps, collapse = ", "),
            ") -- skipped")
    next
  }
  se <- experiments(mae)[[which[1]]]

  # ---- expression: to numeric matrix, symbol-indexed -------------------
  a <- assay(se)
  if (!is.matrix(a)) a <- as.matrix(a)
  mode(a) <- "numeric"

  rd  <- rowData(se)
  sym <- if ("gene_name" %in% colnames(rd)) as.character(rd$gene_name)
         else rownames(a)
  keep <- !is.na(sym) & nzchar(sym)
  a <- a[keep, , drop = FALSE]; sym <- sym[keep]

  # collapse duplicate symbols by maximum mean expression
  mm  <- rowMeans(a, na.rm = TRUE)
  ord <- order(sym, -mm)
  a   <- a[ord, , drop = FALSE]; sym <- sym[ord]
  dup <- duplicated(sym)
  a   <- a[!dup, , drop = FALSE]
  rownames(a) <- sym[!dup]
  a <- a[order(rownames(a)), , drop = FALSE]

  # ---- clinical, aligned to expression columns via the sample map -----
  cd <- as.data.frame(colData(mae), stringsAsFactors = FALSE)
  sm <- as.data.frame(sampleMap(mae), stringsAsFactors = FALSE)
  sm <- sm[sm$assay == which[1], , drop = FALSE]

  # map expression colnames -> primary (patient) identifiers
  prim <- sm$primary[match(colnames(a), sm$colname)]
  ok   <- !is.na(prim)
  a    <- a[, ok, drop = FALSE]
  prim <- prim[ok]

  clin <- cd[prim, , drop = FALSE]
  have <- intersect(CLIN_FIELDS, colnames(clin))
  clin <- clin[, have, drop = FALSE]
  clin$sample_id <- colnames(a)
  clin$cohort    <- cohort
  clin <- clin[, c("cohort", "sample_id", have)]

  # ---- write ------------------------------------------------------------
  ef <- file.path(OUT, paste0(cohort, "_expr.tsv.gz"))
  gz <- gzfile(ef, "w")
  write.table(cbind(gene = rownames(a), as.data.frame(a)), gz,
              sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)

  cf <- file.path(OUT, paste0(cohort, "_clin.tsv"))
  write.table(clin, cf, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "")

  resp <- if ("response" %in% colnames(clin)) clin$response else rep(NA, nrow(clin))
  n_lab <- sum(resp %in% c("R", "NR"))
  n_r   <- sum(resp == "R", na.rm = TRUE)
  n_nr  <- sum(resp == "NR", na.rm = TRUE)
  n_os  <- if ("survival_time_os" %in% colnames(clin))
             sum(!is.na(clin$survival_time_os)) else 0L

  inv[[length(inv) + 1]] <- data.frame(
    cohort = cohort, assay = which[1],
    n_samples = ncol(a), n_genes = nrow(a),
    n_response_labeled = n_lab, n_responder = n_r, n_nonresponder = n_nr,
    minority_class = min(n_r, n_nr),
    n_with_os = n_os,
    cancer_type = paste(sort(unique(na.omit(clin$cancer_type))), collapse = "|"),
    treatment   = paste(sort(unique(na.omit(clin$treatment))),   collapse = "|"),
    stringsAsFactors = FALSE
  )

  message(sprintf("   %d genes x %d samples | labeled %d (R=%d NR=%d) | OS %d",
                  nrow(a), ncol(a), n_lab, n_r, n_nr, n_os))
}

inventory <- do.call(rbind, inv)
inventory <- inventory[order(-inventory$n_response_labeled), ]
write.csv(inventory, file.path(RES, "cohort_inventory.csv"), row.names = FALSE)

message("\nHarmonized ", nrow(inventory), " cohorts -> ", OUT)
message("Inventory: ", file.path(RES, "cohort_inventory.csv"))
