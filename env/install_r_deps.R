#!/usr/bin/env Rscript
# R dependencies for steps 00, 01 and 07.
#
# Bioconductor supplies the MultiAssayExperiment / SummarizedExperiment
# classes the ORCESTRA cohort objects are serialized as; metafor supplies
# the random-effects meta-analysis.
#
# Note: GenomeInfoDb requires the GenomeInfoDbData annotation package,
# which some conda builds ship without installing. If library(GenomeInfoDb)
# fails with "GenomeInfoDbData not found", install it explicitly from the
# current Bioconductor release (as done below) rather than from an
# archived release URL.

cran <- c("jsonlite", "data.table", "metafor")
bioc <- c("GenomeInfoDbData", "MultiAssayExperiment", "SummarizedExperiment")

need <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(need))
  install.packages(need, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

need_b <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(need_b))
  BiocManager::install(need_b, ask = FALSE, update = FALSE)

for (p in c(cran, bioc[-1]))
  cat(sprintf("%-24s %s\n", p,
              as.character(utils::packageVersion(p))))
cat("\nR:", R.version.string, "\n")
