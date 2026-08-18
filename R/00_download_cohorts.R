#!/usr/bin/env Rscript
# Step 00 -- download the PredictioR/ORCESTRA ICB compendium.
#
# The 21 public ICI-treated cohorts are distributed by ORCESTRA as
# MultiAssayExperiment .rds objects hosted on Zenodo, each with its own
# DOI. This script resolves the current download manifest from the
# ORCESTRA API, records it (name, DOI, URL, checksum) for provenance, and
# fetches any cohort not already present.
#
# Nothing here is analysis: this step only puts bytes on disk. Re-running
# it is a no-op for cohorts already downloaded, so it is safe to resume
# after an interrupted download.
#
# Usage:  Rscript R/00_download_cohorts.R [--refresh-manifest]

suppressPackageStartupMessages({
  library(jsonlite)
  library(tools)
})

source(file.path(dirname(sub("^--file=", "",
  c(grep("^--file=", commandArgs(FALSE), value = TRUE), "./R/x")[1])),
  "_root.R"))

ROOT <- repo_root()

RAW_DIR  <- file.path(ROOT, "data", "raw")
CFG_DIR  <- file.path(ROOT, "config")
MANIFEST <- file.path(CFG_DIR, "orcestra_manifest.json")
API      <- "https://www.orcestra.ca/api/clinical_icb/canonical"

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
refresh <- "--refresh-manifest" %in% args

# ---- manifest -------------------------------------------------------------
# The committed manifest is the version the published analysis used. Use
# --refresh-manifest to pull the current one; if ORCESTRA has re-versioned
# a cohort the DOIs will differ, which is a deviation worth recording.
if (refresh || !file.exists(MANIFEST)) {
  message("Fetching cohort manifest from ORCESTRA ...")
  api <- fromJSON(API, simplifyDataFrame = FALSE)
  man <- lapply(api, function(e) list(
    name         = e$name,
    doi          = e$doi,
    download_url = e$downloadLink,
    datatypes    = paste(vapply(e$availableDatatypes,
                                function(d) d$name, character(1)),
                         collapse = ";")
  ))
  man <- man[order(vapply(man, function(x) x$name, character(1)))]
  write(toJSON(man, auto_unbox = TRUE, pretty = TRUE), MANIFEST)
  message("Wrote ", MANIFEST, " (", length(man), " cohorts)")
}

man <- fromJSON(MANIFEST, simplifyDataFrame = FALSE)
message("Manifest lists ", length(man), " cohorts.")

# ---- download -------------------------------------------------------------
rows   <- list()
failed <- character(0)
for (e in man) {
  dest <- file.path(RAW_DIR, paste0(e$name, ".rds"))

  # A cohort counts as present only if it loads. Size alone is not enough:
  # an interrupted download leaves a large partial file that a later run
  # would otherwise accept, and the failure would surface much later as a
  # corrupt object in step 01.
  have <- file.exists(dest) && file.size(dest) > 1e6 &&
    !inherits(try(readRDS(dest), silent = TRUE), "try-error")

  if (have) {
    message(sprintf("  [have] %-20s %6.0f MB", e$name,
                    file.size(dest) / 1e6))
  } else {
    if (file.exists(dest))
      message(sprintf("  [redo] %-20s (present but unreadable)", e$name))
    message(sprintf("  [get ] %-20s %s", e$name, e$doi))

    # Download to a temporary name and move into place only after the
    # object loads, so an interrupted run never leaves a file that looks
    # complete.
    tmp <- paste0(dest, ".part")
    ok <- tryCatch({
      download.file(e$download_url, tmp, mode = "wb", quiet = TRUE)
      if (!file.exists(tmp) || file.size(tmp) < 1e6)
        stop("downloaded file is missing or implausibly small")
      if (inherits(try(readRDS(tmp), silent = TRUE), "try-error"))
        stop("downloaded file is not a readable .rds object")
      TRUE
    }, error = function(err) {
      message("        FAILED: ", conditionMessage(err))
      FALSE
    })

    if (!ok) {
      unlink(tmp)
      failed <- c(failed, e$name)
      next
    }
    file.rename(tmp, dest)
    message(sprintf("        ok    %6.0f MB", file.size(dest) / 1e6))
  }

  rows[[length(rows) + 1]] <- data.frame(
    cohort    = e$name,
    doi       = e$doi,
    url       = e$download_url,
    bytes     = file.size(dest),
    md5       = as.character(md5sum(dest)),
    retrieved = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
}

prov <- do.call(rbind, rows)
out  <- file.path(ROOT, "results", "cohort_download_provenance.csv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(prov, out, row.names = FALSE)

message("\nDownloaded/verified ", nrow(prov), " of ", length(man), " cohorts.")
message("Provenance (DOI + md5 per cohort): ", out)

if (length(failed)) {
  message("\n", length(failed), " cohort(s) did not download: ",
          paste(failed, collapse = ", "))
  message("Re-run this step to retry them -- completed cohorts are skipped, ",
          "so a retry only fetches what is missing. Downstream steps will ",
          "silently analyze the smaller set, so do not proceed until this ",
          "reports 0 failures.")
  quit(status = 1)
}
