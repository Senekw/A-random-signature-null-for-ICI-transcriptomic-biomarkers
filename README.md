# A random-signature null for ICI transcriptomic biomarkers

Re-runnable null-calibration harness for

> **Most published immunotherapy-response gene signatures do not outperform
> random gene sets, and collapse onto a single tumor-inflammation axis**
> S. Lee. *Submitted, Journal for ImmunoTherapy of Cancer.*

Dozens of transcriptomic signatures are proposed to predict response to
immune-checkpoint inhibitors. They are rarely benchmarked against a
calibrated null, so it is unclear whether their apparent performance
reflects specific biology or the generic behaviour of any large
immune-weighted gene set. This repository is the benchmark, and it is
built so that a new signature can be put through the same test in a few
lines.

---

## The test, in five lines

```python
from icinull import load_cohort, random_null_test

cohort = load_cohort("Mariathasan")          # harmonized ICI cohort
my_genes = ["CD8A", "GZMB", "PRF1", "IFNG"]  # your candidate signature

random_null_test(cohort.expr, cohort.labels(), my_genes,
                 universe=cohort.universe, n_draws=1000, seed=20260705)
# NullResult(observed_auroc=..., null_mean=..., null_q95=...,
#            p_empirical=..., z_vs_null=...)
```

The question it answers is not "is my signature associated with
response". It is: **does my signature beat gene sets of the same size,
drawn at random from the same expressed universe, scored the same way, in
the same cohort?**

If `p_empirical > 0.05`, the signature has not cleared the bar that
matters.

---

## What the pipeline does

| Step | Script | Output |
|------|--------|--------|
| 00 | `R/00_download_cohorts.R` | ICB cohorts from ORCESTRA (+ DOI/md5 provenance) |
| 01 | `R/01_harmonize.R` | `data/harmonized/*_expr.tsv.gz`, `*_clin.tsv`, `results/cohort_inventory.csv` |
| 01b | `scripts/01_build_signatures.py` | `signatures/signatures.json`, `signature_provenance.csv` |
| 02 | `scripts/02_cohort_gate.py` | `results/cohort_gate.csv` |
| 03 | `scripts/03_score_signatures.py` | `results/perf_per_sig_cohort.csv`, per-cohort score matrices |
| 04 | `scripts/04_null_calibration.py` | `results/null_results.csv`, `null_pooled.csv` |
| 05 | `scripts/05_axis_and_confound.py` | `results/axis_pca_*.csv`, `venet_alignment_*.csv`, `estimate_scores.csv`, `partialled_auroc_*.csv` |
| 06 | `scripts/06_transferable_ceiling.py` | `results/ceiling_loco.csv`, `optimism.csv` |
| 07 | `R/07_meta_analysis.R` | `results/meta_results.csv`, `meta_loco.csv`, forest plots |
| 08 | `scripts/08_headline_numbers.py` | `results/headline_numbers.json` |
| 09 | `scripts/09_figures.py` | `figures/figure*.png` |
| 10 | `scripts/10_check_manuscript.py` | reconciliation report vs the paper |

Every number quoted in the manuscript is regenerated into
`results/headline_numbers.json`, each entry tagged with the CSV and column
it came from. Figures read only those committed CSVs, so a figure cannot
disagree with a table.

Step 10 closes the loop: `config/manuscript_claims.yaml` lists all 37
numbers the paper quotes, each with its path in `headline_numbers.json` and
a per-claim tolerance, and `make check` reports any that have moved beyond
it. "Does this repository reproduce the paper" therefore has a mechanical
answer — `make check` exits non-zero if it does not.

**Cohort count matters for the ceiling.** Steps 02–05 and 07 are
per-signature or per-cohort and behave sensibly on any subset. Step 06
(leave-one-cohort-out) does not: with fewer than three cohorts each fold
trains on one cohort, and the trained-panel rows are fold noise that can
fall below chance and order arbitrarily. The script warns when this
applies. Do not read a directional ceiling result off a subset run.

---

## Install and run

```bash
# Python 3.11+
conda env create -f env/environment.yml && conda activate icinull
pip install -e .

# R 4.4 with Bioconductor
Rscript env/install_r_deps.R

# 1. fetch cohorts (~0.9 GB). Re-run until it reports 0 failures; it skips
#    what is already present and exits non-zero while any are missing.
make download

# 2. fast end-to-end check (100 draws) before committing to the full run
make smoke

# 3. the real thing, ending in a reconciliation report against the paper
make all
```

`make test` runs the unit tests (29 tests, no download required — they
build small synthetic matrices with known answers, including a check that
the null is *calibrated*: random gene sets do not beat it above the
nominal rate).

Full pipeline runtime is dominated by step 04 (1000 random draws ×
102 signatures × cohorts). On 8 cores expect a few hours; use
`--draws 100` for a fast smoke run and `--jobs N` to parallelize.

---

## Data provenance

Cohorts are the public **PredictioR / ORCESTRA** ICB compendium — 21
uniformly-processed ICI-treated transcriptomic cohorts, each distributed
as a `MultiAssayExperiment` with its own Zenodo DOI.
`config/orcestra_manifest.json` pins the name, DOI and download URL of
every cohort as retrieved for this analysis;
`R/00_download_cohorts.R` records a per-file md5 at download time into
`results/cohort_download_provenance.csv`.

Signature membership is **assembled from upstream sources at build time,
never transcribed here**: canonical predictors from the curated CSVs in
[`bhklab/SignatureSets`](https://github.com/bhklab/SignatureSets), the 90
immune/immuno-oncology sets from
[IOBR](https://github.com/IOBR/IOBR)'s `signature_collection` restricted
to `config/signature_roster_iobr.txt`.
`scripts/verify_roster.py` re-derives the roster from the current upstream
release and fails loudly if IOBR has renamed or removed a set.

The ESTIMATE 141-gene immune and stromal signatures are committed as
reference data (`config/estimate_gene_sets.csv`) with their source and
symbol-update provenance in the file header.

No patient-level data is redistributed in this repository. `data/` is
populated by step 00 from the DOIs above and is git-ignored.

---

## Pre-specification

[`docs/methods_prespec.md`](docs/methods_prespec.md) records every analysis
decision that was fixed before performance metrics were computed — the
cohort gate, the response definition, the scoring methods, direction
assignment, size matching, the seed — plus a deviations section and an
explicit statement of what the design cannot rule out.

The machine-readable form is
[`config/analysis_config.yaml`](config/analysis_config.yaml). Changing a
value there invalidates the pre-specification and must be recorded as a
deviation.

For evaluating your own signature against the null, see
[`docs/using_the_harness.md`](docs/using_the_harness.md).

---

## Provenance of this code

This harness was written to the manuscript's Methods specification and
verified against the reported cohort composition and per-cohort AUROCs
(e.g. the T-cell-inflamed GEP in the Kim gastric cohort: 31 labeled
samples, 13 responders, AUROC 0.79 — matching the manuscript's forest
plot). It is a clean-room implementation, not a copy of a historical
script: results regenerate from this pipeline and are internally
consistent, but they are not byte-identical to a run of different code.

`results/` and `figures/` in a fresh clone are empty. Run the pipeline to
populate them; that is the point.

---

## Citation

See [`CITATION.cff`](CITATION.cff). Please cite the upstream data and
signature resources as well — PredictioR/ORCESTRA for the cohorts, IOBR
and SignatureSets for the gene sets, and the primary publication of each
canonical signature (PMIDs in `signatures/signature_provenance.csv`).

## Licence

Code: MIT (`LICENSE`). Assembled signature membership and cohort data
remain under their upstream licences.
