# What the CPTAC covariance-preserving experiment means

This is the plain-language guide to the completed **CPTAC-TUMOR-COV-NULL-001** experiment. The [technical report](../results/cptac_tumor_covariance_null_20261009/REPORT.md), [frozen protocol](../audit/cptac_covariance_null/PROTOCOL.md), and [raw 499 results](../results/cptac_tumor_covariance_null_20261009/499_draws.json) contain the full record. This experiment was conducted **after** the original preregistration; it does not replace the original study.

## Where is the experiment?

- [Implementation](../audit/cptac_covariance_null/experiment.py) and [numerical tests](../audit/cptac_covariance_null/test_experiment.py) — code that computes the statistic and verifies its invariants.
- [499-draw GitHub Actions run](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/37959356460) and [seeded replay](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/37959908996) — executions; the ordered null values were reproduced exactly.
- [Machine-readable output](../results/cptac_tumor_covariance_null_20261009/499_draws.json) — every null statistic, metadata, input hashes, and summary.
- [Methodological correction](../paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md) — important qualifications to some historical claims.

## What question does it ask?

A persistent-homology **H1 feature** is a loop-like pattern in a point cloud. Its *persistence* records how long the feature lasts as the geometric connection scale changes. The experiment measures the **longest finite H1 interval** in a fifty-dimensional principal component analysis (PCA) representation of cancer proteomics samples. It does not image literal rings in tumors.

The question is whether that longest interval is unusually large compared with data that retain the **entire cross-protein covariance matrix**. If covariance and its corresponding PCA spectrum can produce similarly persistent H1 features under the specified reference model, a large real-data statistic does not, by itself, isolate topology beyond that second-order structure.

## What did we compute?

We used **110 primary-tumor samples** from the Clinical Proteomic Tumor Analysis Consortium clear-cell renal cell carcinoma (CPTAC-CCRCC) proteomics dataset, study PDC000127. The analysis starts from 2,000 features selected and standardized on the *original mixed 194-sample cohort*, then restricts the matrix to the 110 tumor rows. This detail matters: the experiment did **not** select a new feature panel using only tumor samples.

We computed the observed maximum H1 persistence on the tumor-only PCA(50) coordinates. We then generated **499 sample-space rotations** and measured the same statistic on each. By construction, these reference matrices keep the tumor-feature means, **all feature-by-feature covariance values**, and PCA eigenvalue spectrum fixed. They do **not** keep individual patients' positions, feature marginal distributions, batch labels, or other clinical dependencies fixed.

This is a *conditional covariance-preserving reference*, not 499 newly sampled tumors or an independent biological replication.

## What does the actual distribution show?

![Histogram of all 499 covariance-preserving null H1 persistence values, with the observed statistic and the 45-draw right tail marked.](../results/cptac_tumor_covariance_null_20261009/null_distribution.svg)

*This figure was generated from the [archived 499 values](../results/cptac_tumor_covariance_null_20261009/499_draws.json), not from the conceptual circles shown in the chat. The colored tail begins exactly at the observed value; every amber draw counts toward the empirical tail probability.*

| Measure | Recorded result |
| --- | ---: |
| Observed maximum H1 persistence | **5.3543** |
| Null average | **4.4693** |
| Null sample standard deviation | **0.6428** |
| Null draws at least as high as observed | **45 of 499** |
| One-sided Monte Carlo p-value, with plus-one correction | **(1 + 45) / (499 + 1) = 0.092** |
| Prespecified sensitivity criterion | **p ≤ 0.05** |

The observed value is above the null **mean**, but 45 of the 499 covariance-preserving controls were at least as large. A bar chart showing only “5.354 versus 4.469” hides the width and right tail of the reference distribution; the histogram is the appropriate visual comparison.

## What does p = 0.092 mean?

The result **did not meet** the prespecified 0.05 threshold. It is a *non-rejection* under the specified reference model, **not proof that covariance fully explains the data**. The p-value is not the probability that biology is absent, that the null hypothesis is true, or that the observed loop is fake.

The interpretation is conditional on the experiment's exchangeable-row, matrix-normal reference. Patient/sample identity, batch effects, selection on the mixed cohort, and non-Gaussian measurement structure were not independently eliminated. The analysis establishes neither a biological mechanism nor the absence of one.

## Why is the older p = 0.002 not a contradiction?

The earlier [CPTAC matched-transform sensitivity](../results/cptac_resid_null_sensitivity_20261009/REPORT.md) reported **p = 0.002**. But it used **194 tumor-plus-normal samples after class-mean residualization** and independently permuted protein values **within each class**, destroying cross-protein dependence.

The newer experiment uses **110 tumors only**, no class-mean residualization, and a covariance-preserving rotation. **Two variables changed: the population and the null model.** Comparing the resulting p-values does not identify which change caused the difference. Neither result supersedes the other's archived calculation.

## What is the most informative next test?

Hold the **same 110 tumor samples, feature selection, PCA(50), H1 statistic, and code version fixed**, and compare an independent-feature permutation reference with the covariance-preserving reference. This would isolate the effect of the null construction **within this tumor-only analysis**. It would still not establish that any remaining structure is biological without checking patient-level independence and technical batches.

## How can another researcher regenerate the chart?

From a clean checkout with Python 3.11, run:

```bash
python audit/cptac_covariance_null/plot_null_distribution.py --check
```

The check generates the expected chart deterministically from the committed JSON and fails if it differs from the checked-in SVG. To regenerate the checked-in visualization deliberately, omit `--check`. This does **not** repeat the original 499-draw computation. For full numerical reruns and assumptions, use the [technical report](../results/cptac_tumor_covariance_null_20261009/REPORT.md) and [usage guide](USAGE.md).
