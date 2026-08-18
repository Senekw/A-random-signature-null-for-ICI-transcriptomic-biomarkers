# Reference run

These tables are the output of a full-compendium run:

* all 21 ORCESTRA cohorts downloaded from their Zenodo DOIs
  (`cohort_download_provenance.csv` records per-file md5s);
* 15 passing the pre-specified gate (`cohort_gate.csv`, with the 6
  exclusions and their failing criteria);
* 102 signatures, 1,000 null draws per signature x cohort, seed 20260705;
* `reconciliation_report.csv` / `docs/reproduction_report.md` compare every
  number here against the submitted manuscript.

They are committed so a reader can inspect the results and check the
reconciliation without a multi-hour run. They are not a substitute for
running the pipeline: `make clean-results && make all` regenerates
everything, and `make check` re-verifies it against
`config/manuscript_claims.yaml`.

Per-cohort score matrices (`results/scores/`) are not committed -- they are
large and every table here derives from them. Step 03 regenerates them.

No patient-level data is present in these files: they are per-signature and
per-cohort summary statistics. `estimate_scores.csv` is per-sample but
contains only derived enrichment scores from public expression data.
