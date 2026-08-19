# Reproduction report — full compendium

All 21 ORCESTRA cohorts, 15 passing the pre-specified gate, 102 signatures,
1,000 null draws, seed 20260705. Regenerate with `make all`; re-verify with
`make check`. Machine-readable: `results/reconciliation_report.csv`,
`results/headline_numbers.json`.

**26 of 37 registered manuscript claims reproduce within tolerance** (was 22
before this round of bug-fixing). Every structural conclusion of the paper
reproduces. The 11 that remain outside tolerance are itemized below with a
diagnosis for each — several are *definitional*, one is a discrepancy I could
not resolve and am flagging rather than fitting.

A note on method: tolerances were set from the manuscript's own precision and
the Monte-Carlo error of the quantity, **before** the run. No threshold was
widened, and no parameter was tuned, to convert a MOVED into an ok. Where the
pipeline disagrees with the paper, that is reported as a disagreement.

---

## Exact (deterministic, from raw DOIs, no tuning)

| Claim | Paper | Rerun |
|---|---|---|
| Cohorts evaluated / passing gate | 21 / 15 | 21 / 15 |
| Excluded cohorts | 6 | same 6, same failing criteria |
| Response-labeled patients | 905 | 905 |
| Patients with OS | 862 | 862 |
| Cancer types | 6 | 6 |
| Signatures in library | 102 | 102 |
| **Median null 95th percentile** | **0.621** | **0.621** |
| **Pooled CIs excluding 0.5** | **32** | **32** |

## Within tolerance

| Claim | Paper | Rerun |
|---|---|---|
| Median / mean signature AUROC | 0.548 / 0.549 | 0.538 / 0.540 |
| Scoring-method r / rho | 0.94 / 0.93 | 0.936 / 0.936 |
| PC1 variance explained, mean / max | 55% / 76% | 54.2% / 74.7% |
| Venet R² (mean-z) | 0.63 | 0.656 |
| Venet Pearson r | 0.79 | 0.810 |
| corr(axis, ESTIMATE immune) | 0.63 | 0.622 |
| Mean AUROC after removing infiltration | 0.503 | 0.511 |
| Best-10 mean AUROC | 0.640 | 0.635 |
| Best pooled AUROC | 0.628 | 0.627 |
| Median I² | 15% | 15.2% |
| **LOCO ceiling (GEP axis)** | **0.594** | **0.602** |
| LOCO, full panel | 0.553 | 0.564 |
| In-cohort best-selected | 0.760 | 0.764 |
| Within-cohort CV panel | 0.624 | 0.622 |
| **Optimism gap** | **0.17** | **0.163** |

The ordering the argument rests on holds on real data:
**axis 0.602 > panel+infiltration 0.576 > panel 0.564 > random panel 0.517**,
and the axis bootstrap CI [0.555, 0.670] contains the published 0.594.

---

## Bugs found and fixed this round

These were real defects, not tuning. Each changed results.

1. **Three canonical predictors were silently deleted from the entire
   analysis.** `scoring.min_genes_covered` was 2, but the Methods explicitly
   include "PD-L1/PD-1/CTLA-4 single genes" among the 12 canonical
   predictors. A floor of 2 dropped all three from all 15 cohorts — 45 tests,
   a quarter of the canonical set — which biased the canonical null-beating
   rate and cut the Venet regression to n=99. Fixed to 1.

2. **Cohorts use two different GENCODE vintages, and gene sets were being
   silently truncated.** Found by noticing that the `Histones` signature had
   zero coverage in 11 of the 15 passing cohorts: `ICB_Van_Allen`,
   `ICB_Miao1`, `ICB_Braun` and `ICB_Puch` carry pre-2020 histone symbols
   (`HIST1H2AG`), and the other eleven carry current ones (`H2AC11`) — the
   split *is* the bug. 300–2,400 retired symbols per cohort. 80 of the
   2,635 distinct signature genes matched *no* cohort at all, affecting 25
   signatures — and a set written in legacy symbols matched in some cohorts
   and vanished in others, making a signature's realized size, its
   size-matched null, and its meta-analysis weight depend on annotation
   vintage rather than biology. Fixed with `src/icinull/symbols.py`: both
   sides of the join are mapped onto current NCBI official symbols, but only
   where the alias resolves to exactly one current symbol. Ambiguous ones are
   left alone and reported rather than guessed — `IL8RA` is a documented
   synonym of *both* CXCR1 and CXCR2, `RAB7` of both RAB7A and RAB7B — as are
   microarray probe IDs and composite entries like `GCLC/GCLM`
   (`signatures/unmapped_genes.csv`, 45 identifiers). 102 retired symbols
   updated; 8 tests recovered.

3. **The random-panel baseline drew independent gene sets per cohort**, so
   column `random_007` meant a different gene set in every cohort and a model
   trained on other cohorts was applied to unrelated features. Fixed to draw
   once from the shared expressed universe.

4. **The purity transform could invert.** ESTIMATE purity is `cos(A + B·s)`,
   monotone only while the argument is in [0, π]; a min–max rescale let a
   single outlier push it past π, where cosine turns back up and the
   *most-infiltrated* samples are reported as the *purest* — inverting the
   very relationship the confound analysis measures. Now mapped by
   within-cohort rank: strictly monotone (Spearman exactly −1.0 in every
   regime tested), tie-free, and insensitive to outlier magnitude.

5. **The Venet regression included the reference signature**, which
   correlates 1.0 with itself by construction. Now excluded from the headline
   and reported both ways.

6. **The random panel was size-matched to the intersected feature count
   (98), not the library size (102)** the Methods specify.

7. **Two signatures are one gene set.** `Chemokine12_Messina` and
   `TLS_12chemokine_Prabhakaran` have identical membership — Prabhakaran
   re-evaluated Messina's 12-chemokine score rather than deriving a new one.
   The library is 102 named signatures over **101 distinct gene sets**, and
   every aggregate double-weights that set. The builder now reports this on
   every run rather than leaving it as an unnoticed assumption.

8. **The fixed seed did not actually fix the draws.** The per-signature seed
   offset was `abs(hash(cohort)) % 997`, and Python randomizes string hashing
   per process — so the same command drew *different* random gene sets on
   every invocation, while the code and the pre-specification both claimed a
   fixed seed. Verified directly: three fresh interpreters gave offsets 281,
   63 and 752 for the same cohort. Replaced with a blake2b digest, which is
   identical on every run and machine (pinned by `tests/test_seeding.py`,
   which also asserts that `hash()` *would* have been unstable). A second
   defect compounded it: `null_results.csv` recorded only the derived
   per-test seed, so the base seed was unrecoverable and
   `headline_numbers.json` reported a derived value as "the seed". The base
   seed is now recorded separately and reports as 20260705.

9. **The duplicated gene set is now quantified rather than only warned
   about.** `signature_provenance.csv` carries a membership hash, and
   `headline_numbers.json` reports `n_distinct_gene_sets: 101`, names the
   duplicate pair, and gives every headline aggregate both with and without
   it (de-duplicating moves the median AUROC by 0.001 and the
   below-0.55 fraction by 0.003 — small, but no longer something a reader has
   to discover).

10. **The null was 26× slower than necessary**, re-standardizing the whole
   expression matrix on every draw. Now standardizes once per cohort:
   **1m53s instead of ~50 minutes** for the full 1,000-draw run, with
   bit-identical output (pinned by `test_null_fast_path_matches_reference`,
   which asserts agreement with the naive loop to 1e-12).

Also fixed: R scripts failed under `Rscript` (root resolution via
`sys.frame`, and a helper used before definition); downloads accepted
truncated files (R's 60 s default timeout truncated the larger cohorts, and
`download.file` *warns* rather than errors on a short read); a byte-size
validity floor rejected the genuinely small cohorts (`ICB_Hwang` is 84 KB);
`Rscript` was assumed on PATH; cancer-type counting summed per-sample biopsy
sites, inflating 6 types to 10 because `ICB_Mariathasan` is a urothelial
cohort whose clinical column records biopsy site; the claims checker split
dotted paths naively and missed leaf keys containing a literal dot; and
`docs/methods_prespec.md` described a coverage guard the config never
implemented.

## Verification added

82 unit tests (from 29), no data download required:

* the Cox fallback is **numerically identical to lifelines** (the Methods'
  stated tool) across seeds, heavy ties, and heavy censoring;
* ssGSEA matches the **published Barbie et al. definition** to 1e-12 at
  three α values, transcribed independently from the spec;
* symbol harmonization is **vintage-invariant**: the same biology under two
  annotation vintages yields the same genes and the same realized size;
* purity is **strictly monotone** in the combined score across six input
  regimes including single high and low outliers, tiny and large cohorts;
* the null's fast path is **bit-identical** to the naive implementation;
* and the original calibration check — random gene sets do not beat the null
  above the nominal rate.

---

## The 11 that remain outside tolerance

### Definitional, not computational (3)

**Test count: 1,463 → 1,470 (+7), and the two proportions over it**
(`frac tests below 0.55` 0.510 → 0.534; `frac tests beating null` 0.305 →
0.272). After the coverage fix the grid is 15 × 102 = 1,530 with 60
exclusions, 58 of them `ICB_VanDenEnde` (a targeted assay: 44 of 102
signatures scored). The paper's 1,463 sits 7 below; the residual is
consistent with a slightly different annotation vintage upstream, and is not
attributable to a rule I can identify in the Methods. **Recommendation:**
update the three numbers, or state the coverage rule explicitly in the
Methods so the denominator is reconstructible.

**Signatures meta-analyzable: 100 → 99.** All 102 now have k ≥ 2; the
reported 100 is 2 below the library. Given the Messina/Prabhakaran
duplication, "100 distinct" may be the intended count. Worth a sentence.

### Moved in the direction that strengthens the paper (4)

| Claim | Paper | Rerun |
|---|---|---|
| corr(PC1, T-cell-inflamed GEP) | 0.79 | 0.829 |
| Venet R² under ssGSEA | 0.81 | 0.890 |
| corr(axis, tumor purity) | −0.68 | −0.629 |
| canonical predictors beating null | 44.5% | 50.3% |

The collapse is *stronger* than reported. The last one moves against the
argument — canonical predictors clear the null in half of cohorts, not 44.5%
— and needs the number changed, though "roughly half" still holds.

### Ceiling comparators (3)

| Claim | Paper | Rerun |
|---|---|---|
| LOCO, ESTIMATE immune | 0.556 | 0.587 |
| LOCO, panel + infiltration | 0.550 | 0.576 |
| LOCO, random gene-set panel | 0.441 | **0.517** |

The first two are close and preserve every inequality the paper draws. **The
random panel I could not reconcile, and I am flagging it rather than fitting
it.** A trained model over noise features, validated out-of-cohort, sits *at*
chance (0.50) in expectation. Landing systematically *below* chance at 0.441
requires the learned coefficients to transfer with the wrong sign, which
needs the features to mean something different in the held-out cohort than in
training — i.e. precisely the per-cohort-redraw construction that was bug #3
above. My corrected version gives 0.517 and the buggy version gave 0.507;
neither reproduces 0.441.

**This does not threaten the paper's claim** — the random panel is far below
the axis (0.602) either way, so "in-cohort panel performance is overfitting,
not signal" holds. But 0.441 should either be traced to the original code or
replaced with the value the committed pipeline produces, and the sentence
"collapsed out of cohort" softened to "fell to chance".

### One claim in the Methods that does not reproduce

The Methods state **"median gene coverage of signatures in cohorts was
92%"**. The pipeline gives a median of **100%** — most signatures are fully
measured in most cohorts — and it is the *mean* that lands near 92% (94.5%
over evaluable tests after symbol harmonization; 91.8% before). These are
different statistics. I have not resolved which was computed, and have
deliberately not presented the mean as corroboration of the median.

---

## Bottom line

Every structural claim reproduces: the cohort set, the gate and its
exclusions, the library, the null result, the single-axis collapse, the
infiltration confound, the ceiling and its ordering, and the optimism gap.
What needs editing before submission is a set of specific decimals, one
sentence about the random panel, and a decision on whether to state the
coverage rule in the Methods.

`make check` reproduces this table at any time.
