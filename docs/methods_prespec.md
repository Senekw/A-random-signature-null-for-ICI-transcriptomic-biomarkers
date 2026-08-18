# Pre-specified analysis plan

This document records the analysis decisions that were fixed **before any
performance metric was computed**, and every deviation from them. It is
the audit trail for the claim that the null result is not the product of
analytic choices made after seeing the answers.

The machine-readable form of everything below is
[`config/analysis_config.yaml`](../config/analysis_config.yaml). Where a
value appears in both, the config file is authoritative and this document
explains it.

---

## 1. Why a pre-specification is needed here

The paper's central claim is negative: most published ICI-response
signatures do not beat size-matched random gene sets. A negative result
from a benchmark is only as credible as the benchmark's resistance to
tuning. Almost every knob in this pipeline — which cohorts count, how
response is dichotomized, how a signature is scored, what "the same size"
means, which universe random sets are drawn from — could be turned to
make the published signatures look better or worse.

So each knob was set once, from a stated principle, and recorded here.

---

## 2. Cohort inclusion (the gate)

Cohorts come from the PredictioR/ORCESTRA ICB compendium
(`config/orcestra_manifest.json` pins name, DOI and download URL for all
21 cohorts as retrieved). Three criteria, applied before any signature was
scored:

| ID | Criterion | Threshold | Rationale |
|----|-----------|-----------|-----------|
| G1 | Response-labeled samples | ≥ 20 | Below this, an AUROC's standard error exceeds the effect sizes under test; a cohort contributes noise, not evidence. |
| G2 | Minority response class | ≥ 5 | An AUROC with < 5 in either class is determined by a handful of pairs. |
| G3 | Measured genes | ≥ 10,000 | Published signatures and random gene sets must be drawn from a comparable universe. A targeted panel measures the signature's genes but cannot supply a fair random draw. |

Both `n_response_labeled` and the gene count are computed from the
harmonized matrix, not from the publication's reported N, so a sample
dropped in processing is dropped from the gate too.

Every cohort's criterion values — pass or fail — are written to
`results/cohort_gate.csv`. Exclusions are reported with their failing
criterion rather than omitted.

**A note on G3 in practice.** At least one compendium cohort
(`ICB_VanDenEnde`) is distributed as a sparse assay: ~31,000 rows of which
1,365 are protein-coding, so most signature gene sets have near-zero
coverage. It passes G1/G2 and is retained for the per-signature and null
analyses (each test stands alone, and the coverage guard in §5 handles
low-coverage tests), but it is reported and excluded from the multi-
signature **panel** analysis, where a shared feature set across cohorts is
required and a sparse cohort would silently redefine the predictor. That
exclusion is printed by `scripts/06_transferable_ceiling.py` on every run.

---

## 3. Response definition

Primary endpoint: the compendium's own binary `response` field (R / NR),
which encodes the RECIST dichotomy CR+PR versus SD+PD as curated
per-cohort by the compendium authors.

Using the compendium's field rather than re-deriving from raw RECIST
strings is deliberate: it keeps the labels identical to what any other
user of the same compendium would obtain, and avoids this analysis
inventing a private response definition.

Two sensitivity definitions are supported by the code
(`response_definitions` in the config):

* `strict` — CR+PR versus PD only, dropping SD;
* `recist` — re-derived from the raw `recist` column.

Overall survival is a secondary endpoint, reported as a Cox
log-hazard-ratio per standard deviation of score
(`os_loghr_per_sd` in `results/perf_per_sig_cohort.csv`), for the cohorts
that carry OS.

---

## 4. Expression harmonization

* **Assay**: `expr_gene_tpm` (log2-TPM) where present, else the cohort's
  single expression assay.
* **Gene identity**: HUGO symbol from `rowData$gene_name`. Rows with no
  symbol are dropped.
* **Duplicate symbols**: collapsed to the row with the highest mean
  expression across samples. (Fixed a priori; the alternative — averaging
  — was not tried after the fact.)
* **Expressed universe**: genes detected above the matrix floor in ≥ 20%
  of samples in that cohort. The universe is cohort-specific, because the
  random null must be drawn from the genes that cohort actually measured.

---

## 5. Signature library

102 signatures: 12 canonical named ICI predictors + 90 immune/immuno-
oncology gene sets, disjoint by construction.

Membership is **assembled from upstream sources at build time**, never
transcribed into this repository:

* canonical sets from the curated CSVs in `bhklab/SignatureSets`
  (the signature companion of the PredictioR compendium), plus three
  single-gene checkpoint markers (CD274, PDCD1, CTLA4);
* the 90 from IOBR's `signature_collection`, restricted to the roster in
  `config/signature_roster_iobr.txt` — the union of 13 named `sig_group`
  categories, minus the canonical names.

`scripts/verify_roster.py` re-derives that union from the current upstream
release and fails loudly on drift. The roster is part of the
pre-specification: it is not to be edited to match a newer upstream
without recording the change in §9.

**Direction** is assigned a priori per signature (+1 = higher predicts
response, −1 = higher predicts resistance) from the source publication,
and applied as a sign on the score. Direction is **never** chosen to
maximize AUROC — that would guarantee AUROC ≥ 0.5 and destroy the test.

**Coverage guard**: a signature × cohort test requires ≥ 50% of the
signature's genes present, and ≥ 3 genes (single-gene markers excepted).
Failing combinations are written to `results/perf_per_sig_cohort.csv` with
`excluded_reason`, not silently dropped.

---

## 6. Scoring

* **Primary**: mean of per-gene z-scores across samples within a cohort
  (`mean_z`). Deterministic, invariant to monotone per-gene transforms,
  and the method most signature papers use.
* **Sensitivity**: single-sample GSEA (Barbie et al.), α = 0.25.

Both are computed for every signature × cohort. The manuscript's
conclusions are reported for the primary method with the sensitivity
method as a robustness check (Figure 8A); neither was selected after
comparing outcomes.

---

## 7. The random-signature null

For each signature × cohort:

1. Take the signature's **realized** size — the number of its genes
   actually present in that cohort — not its published size. Matching on
   published size would compare a 15-gene score against 18-gene random
   sets.
2. Draw `n_draws = 1000` gene sets of that size, uniformly without
   replacement, from that cohort's expressed universe.
3. Score each by `mean_z` and compute its response AUROC.
4. One-sided empirical p-value: `(#{null AUROC ≥ observed} + 1) / (n + 1)`.
   The +1 keeps p bounded away from zero and makes the smallest reportable
   p-value explicit (1/1001).

The alternative is one-sided by design. A signature that predicts
response *worse* than random is not evidence for the signature.

Per-signature pooling across cohorts is Stouffer's method on the
per-cohort z-versus-null statistics, with Benjamini-Hochberg FDR control
across signatures.

**Seed**: 20260705, fixed before the first null run. Seed stability is
checked at 1 and 42 (`seed_stability_check`); the reported conclusions do
not depend on the seed, and the check is part of the pipeline rather than
a post-hoc reassurance.

---

## 8. Downstream analyses

* **Single-axis collapse**: PCA on the standardized samples × signature-
  score matrix per cohort. PC1 sign is arbitrary in PCA, so PC1 is
  oriented to correlate positively with the T-cell-inflamed GEP before
  any correlation is reported — an orientation convention, not a
  sign-flip to improve a number.
* **The Venet test**: regression of each signature's mean AUROC on its
  mean correlation with the reference axis. The R² is the fraction of
  between-signature performance variance explained by axis alignment
  alone.
* **Infiltration confound**: ESTIMATE's fixed 141-gene immune and stromal
  signatures (`config/estimate_gene_sets.csv`, provenance in that file's
  header), scored by ssGSEA on the cohort's expressed universe. Each
  signature's AUROC is recomputed after linearly partialling the immune
  score out of its score.
* **Meta-analysis**: random-effects REML on the logit-AUROC scale with
  Hartung-Knapp-Sidik-Jonkman adjustment (conservative for a small number
  of cohorts), AUROC standard errors by Hanley-McNeil, plus I², τ²,
  prediction intervals and leave-one-cohort-out estimates.
* **Transferable ceiling**: leave-one-cohort-out validation. Train on all
  cohorts but one, predict the held-out cohort, pool held-out predictions,
  compute one AUROC. L2-logistic, C = 0.3, standardization fitted on
  training cohorts only. Compared predictors: the single GEP axis, a
  single ESTIMATE immune score, the full 102-signature panel, the panel
  plus infiltration, and a size-matched random-gene-set panel. CI by
  cohort-level (cluster) bootstrap, 2000 resamples.
* **In-cohort optimism**: the best signature *selected within* each
  cohort, and a within-cohort cross-validated panel, reported alongside
  what transfers — the gap is the quantity discovery papers report.

---

## 9. Deviations from the pre-specification

None to date.

Any change to §2–§8 after results existed belongs here, with the date, the
change, the reason, and the effect on the headline numbers. An empty
section is a claim; it is checkable, because every number in the
manuscript is regenerated into `results/headline_numbers.json` from the
committed tables by `scripts/08_headline_numbers.py`.

---

## 10. What this design cannot rule out

Stated plainly, because a pre-specification that only lists strengths is
not useful:

* **Cohort-level confounding.** These are observational cohorts of
  different tumor types, agents, lines of therapy and sequencing
  platforms. Pooling by meta-analysis handles between-cohort heterogeneity
  statistically but cannot make the cohorts exchangeable.
* **The null tests discrimination, not biology.** A signature can be
  mechanistically correct and still fail this test — the test asks whether
  the signature discriminates responders better than random genes do, in
  these cohorts, at these sample sizes.
* **Power.** With 20–300 labeled samples per cohort, an AUROC difference
  of 0.05 is not reliably detectable in a single cohort. That is why the
  central claim is about the *distribution* across 102 signatures × the
  passing cohorts, and why per-signature conclusions are pooled before
  FDR control.
* **Bulk RNA only.** The ceiling reported here is a ceiling for bulk tumor
  transcriptomes with these labels. It says nothing about multimodal,
  single-cell, or spatial predictors.
* **Reimplementation, not re-run of the original code.** See the README
  note on provenance: this harness was written to the Methods
  specification. Numbers regenerate from the committed pipeline; they are
  not byte-identical to a historical run of different code.
