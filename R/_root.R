# Repository-root resolution, shared by all R steps.
#
# Sourced by R/00, R/01 and R/07 so that every step finds the repo the same
# way whether it is run as `Rscript R/<step>.R` from the repo root (the
# documented way), from inside R/, or via source() in an interactive
# session. Failing loudly beats silently writing results to the wrong tree.

repo_root <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  self <- if (length(a)) sub("^--file=", "", a[1]) else NULL
  if (is.null(self)) {
    fr <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (!is.null(fr)) self <- fr
  }
  cands <- character(0)
  if (!is.null(self))
    cands <- c(cands, dirname(dirname(normalizePath(self, mustWork = FALSE))))
  cands <- c(cands, ".", "..")
  for (p in cands)
    if (dir.exists(file.path(p, "config")) &&
        file.exists(file.path(p, "pyproject.toml")))
      return(normalizePath(p, mustWork = FALSE))
  stop("cannot locate the repository root (no config/ + pyproject.toml ",
       "found). Run from the repository root, e.g. ",
       "Rscript R/01_harmonize.R", call. = FALSE)
}
