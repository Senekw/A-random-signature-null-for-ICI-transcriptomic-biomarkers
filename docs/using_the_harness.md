# Testing your own signature against the null

The manuscript's argument is that a new ICI-response signature should have
to clear one bar before it is proposed as a biomarker: **beat gene sets of
the same size, drawn at random from the same expressed universe, scored
the same way, in the same cohorts.** This document is how to do that with
this repository.

---

## 1. One cohort, one signature

```python
from icinull import load_cohort, random_null_test

cohort = load_cohort("Mariathasan")
result = random_null_test(
    cohort.expr,                 # genes x samples, log2-TPM, symbol index
    cohort.labels(),             # 1 = responder, 0 = non-responder
    ["CD8A", "GZMB", "PRF1", "IFNG", "CXCL9"],
    universe=cohort.universe,    # the cohort's expressed genes
    n_draws=1000,
    seed=20260705,
)
print(result)
```

`NullResult` fields:

| field | meaning |
|-------|---------|
| `observed_auroc` | your signature's response AUROC in this cohort |
| `null_mean`, `null_sd` | the random-set AUROC distribution |
| `null_q95` | 95th percentile of the null — the bar to clear |
| `p_empirical` | one-sided `(#{null ≥ observed} + 1) / (n + 1)` |
| `z_vs_null` | standardized distance from the null mean |
| `realized_size` | how many of your genes this cohort measured |

**Read `realized_size` first.** If your 30-gene signature has a realized
size of 11 in this cohort, you are testing an 11-gene signature, and the
null is matched to 11 — which is correct, but it means the cohort is not
really evaluating what you proposed.

---

## 2. All cohorts, with pooling

```python
import numpy as np
from icinull import list_cohorts, load_cohort, random_null_test, stouffer

genes = ["CD8A", "GZMB", "PRF1", "IFNG", "CXCL9"]
rows = []
for name in list_cohorts():
    c = load_cohort(name)
    r = random_null_test(c.expr, c.labels(), genes,
                         universe=c.universe, n_draws=1000, seed=20260705)
    rows.append((name, r))
    print(f"{name:<16} AUROC {r.observed_auroc:.3f}  "
          f"null q95 {r.null_q95:.3f}  p {r.p_empirical:.4f}")

z, p = stouffer([r.z_vs_null for _, r in rows])
print(f"\npooled across {len(rows)} cohorts: z = {z:.2f}, p = {p:.3g}")
```

A signature that beats the null in one cohort and not the others has not
been validated; it has been discovered. The pooled statistic is the claim
worth reporting.

---

## 3. Is your signature actually new?

Beating the null is necessary but not sufficient. The paper's second
finding is that published signatures collapse onto one axis — tumor immune
infiltration. A signature that beats the null *by being a good
infiltration score* adds nothing to the T-cell-inflamed GEP.

Two checks:

```python
import numpy as np
from icinull import load_cohort, load_signatures, mean_z_score, auroc
from icinull.estimate import estimate_scores

c = load_cohort("Mariathasan")
sigs = load_signatures()

mine = mean_z_score(c.expr, genes, 1)
gep = mean_z_score(c.expr, sigs["T_cell_inflamed_GEP_Ayers_et_al"]["genes"], 1)
immune = estimate_scores(c.expr.loc[c.universe])["immune_score"]

print("corr with GEP            :", np.corrcoef(mine, gep)[0, 1].round(3))
print("corr with ESTIMATE immune:", np.corrcoef(mine, immune)[0, 1].round(3))

# incremental discrimination after removing infiltration
y = c.labels()
idx = [s for s in y.index if s in mine.index]
b = np.polyfit(immune.loc[idx], mine.loc[idx], 1)
resid = mine.loc[idx] - np.polyval(b, immune.loc[idx])
print("AUROC raw       :", round(auroc(mine.loc[idx], y.loc[idx].to_numpy()), 3))
print("AUROC partialled:", round(auroc(resid, y.loc[idx].to_numpy()), 3))
```

If `AUROC partialled` falls to ~0.5, your signature is an infiltration
proxy. That may still be useful clinically — but it is not new biology,
and it should be compared against a single fixed infiltration score, not
against chance.

---

## 4. Does it transfer?

The third finding is a ceiling: under leave-one-cohort-out validation, the
whole enterprise saturates at roughly the performance of one inflammation
axis, and in-cohort selection inflates apparent performance substantially.

To place your signature on that scale, add it to the panel and re-run
step 06:

```python
# after scoring your signature into results/scores/<cohort>_mean_z.tsv.gz
# as an extra column, or by adding it to signatures/signatures.json
```

Then compare your held-out number against the `gep_axis` row of
`results/ceiling_loco.csv`. Beating the axis on held-out data across
cohorts is the result that would matter.

---

## 5. Reporting checklist

If you use this harness, report:

1. realized size per cohort (not just published size);
2. per-cohort AUROC **and** the null's 95th percentile in that cohort;
3. the one-sided empirical p-value and the number of draws;
4. the pooled statistic across cohorts, with FDR context if you tested
   more than one candidate;
5. correlation with the T-cell-inflamed GEP and with an infiltration
   score, plus AUROC after partialling infiltration out;
6. the seed.

Item 6 is not pedantry — with 1000 draws the empirical p-value has
Monte-Carlo error of order 0.01, so a p of 0.04 is not stably below 0.05.
Increase `n_draws` if your claim depends on the threshold.

---

## 6. Pitfalls

* **Don't pick direction to maximize AUROC.** Assign it a priori from your
  hypothesis. Choosing the sign after seeing the answer guarantees
  AUROC ≥ 0.5 and voids the test.
* **Don't draw the null from all genes.** Use the cohort's expressed
  universe (`cohort.universe`). Random sets containing undetected genes
  are scored on noise and make the null too easy to beat.
* **Don't compare across scoring methods.** Score your signature and the
  null identically. `mean_z` versus ssGSEA changes AUROC by a few points.
* **Don't test on a cohort that fails the gate.** ≥ 20 labeled samples,
  ≥ 5 in the minority class, ≥ 10,000 measured genes
  (`docs/methods_prespec.md` §2). Below that the AUROC's standard error
  swamps the effect.
