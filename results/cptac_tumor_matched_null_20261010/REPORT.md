# CPTAC tumor-only matched-null comparison — verified result (2026-10-10)

**Experiment:** `CPTAC-TUMOR-MATCHED-NULL-002`  
**Status:** **EXECUTED AND INDEPENDENTLY REPLAYED**; post-preregistration, outcome-informed sensitivity study. This report was written after the new result was observed. The [analysis contract](../../audit/cptac_matched_null/PROTOCOL.md) was frozen before the new null draws. The original July 2026 preregistration and the already archived October covariance-null observations remain unmodified.

## Answer to the stated question

**On an identical tumor-only cohort and analysis statistic, the observed topological persistence is significant against independent per-protein permutations, but it does not cross the 0.05 threshold under the covariance-preserving Haar-rotation null.**

| Exact matching condition | Both references |
| --- | --- |
| Source | CPTAC PDC000127 clear-cell renal carcinoma proteomics |
| Patients/observations | Same **110 tumor observations** from 194 mixed tumor/normal observations |
| Feature panel | Same **2,000 proteins**, selected and standardized once on the original mixed population |
| Computed statistic | Largest finite one-dimensional persistent homology (H1) interval |
| Scores and distance | PCA(50), **Euclidean**, exact singular-value decomposition equivalent |
| Software | Identical `ripser` Vietoris–Rips maximum-H1 computation |
| Observed topological statistic | **5.354297637939453** |

| Quantity | Independently shuffled protein features (new) | Covariance-preserving Haar rotations (archived) |
| --- | ---: | ---: |
| Null samples | 499 | 499 |
| Mean max-H1 | **2.8410664216** | **4.4693335569** |
| Sample SD of max-H1 | 0.5083809031 | 0.6427938924 |
| Null maximum | **4.5843772888** | **6.9796943665** |
| Nulls with H1 >= 5.3542976379 | **0** | **45** |
| One-sided plus-one Monte Carlo p | **0.002** | **0.092** |
| 95% exact binomial Monte Carlo interval for conditional *exceedance probability* | **[0, 0.0073653]** | **[0.0665387, 0.1188068]** |
| Reference threshold `p <= 0.05` | **Reject reference** | **Do not reject reference** |

The independent-feature reference reached the **499-draw Monte Carlo resolution floor** (`(1+0)/500 = 0.002`), not a proof that the idealized tail probability equals exactly 0.002. The 95% intervals describe **simulation uncertainty** for these specific reference distributions, not uncertainty about an underlying biological effect.

## What changed — and what was held fixed?

1. **Feature-permutation null:** independently shuffle the 110 observed values **within each of the same 2,000 protein columns**, using `numpy.default_rng(seed=20261011)`. This preserves every protein's complete tumor-sample empirical marginal distribution, including its exact mean and variance, but destroys its alignment/dependence with other proteins. **Refit PCA(50)** for every shuffled matrix; evaluate the same H1 statistic.
2. **Haar covariance null:** reuse **without recomputation or reseeding** the 499 full-rank sample-space rotations already executed under [`CPTAC-TUMOR-COV-NULL-001`](../cptac_tumor_covariance_null_20261009/REPORT.md), seed `20261010`. These rotations preserve the entire sample-centered cross-protein covariance matrix, feature means and PCA singular-value spectrum, but do **not** preserve each protein's empirical marginal distribution.
3. **Observed analysis:** identical actual tumor samples, feature set, preprocessing lineage, PCA(50), Euclidean point-cloud metric, topological backend, and numerical statistic in both rows. Original mixed-cohort feature selection and standardization are deliberately unchanged, because re-selecting proteins on tumor-only data would introduce a separate change.

The earlier **194-sample mixed tumor/normal class-residualized experiment with `p=0.002` is not part of this matched comparison**. Its sample set and statistic differ; it cannot be substituted for this new tumor-only feature-permutation p-value.

## What conclusion survives critical scrutiny?

**SUPPORTED (conditional comparison):** The observed H1 signal is highly unusual for the null that destroys cross-protein dependence. The same observed statistic is not unusually large at the `0.05` threshold for the distinct null that retains the complete sample-centered covariance. This result is **consistent with ordinary multivariate second-order structure being sufficient to produce strong-looking persistent H1** under the conditional Haar reference.

**NOT ESTABLISHED:** that covariance is the **sole cause**, that true tumor biology contains or lacks topological loops, or that all alternative null models would make the same prediction. The two null families necessarily change **more than covariance alone**: the feature-permutation null retains empirical protein marginals but destroys dependence; Haar rotations retain cross-protein covariance but alter marginals and row distributions. A difference between their p-values is not a mathematical proof of an exclusive causal mechanism.

**Not a preregistered discovery:** The frozen 499-draw seed and `p<=0.05` criterion were decided after seeing the earlier Haar `p=0.092` result. This is therefore a **post-hoc but pre-execution-specified sensitivity comparison**, not a blind confirmatory study. We did not try additional seeds to obtain significance.

## Evidence and independent verification

- **Frozen design:** [`audit/cptac_matched_null/PROTOCOL.md`](../../audit/cptac_matched_null/PROTOCOL.md), committed before generating new tumor-only featurewise nulls.
- **Executable source and synthetic controls:** [`experiment.py`](../../audit/cptac_matched_null/experiment.py) and [`test_experiment.py`](../../audit/cptac_matched_null/test_experiment.py).
- **New primary numerical record:** [all 499 full-precision draws and old-reference summary](499_draws.json); **file SHA-256** `7f64f411c840561fb9468a7edf887836b4a9b93b5d838ad386f43211f5e87b09`. The exact **ordered 499-value feature-null array** serialized as compact decimal JSON has SHA-256 `b3b93df2c69eb62b95964c422dd6b88886cd99df8b5398d2b27fa4a67b945cb4`.
- **Original matched covariance-preserving comparator:** [archived 499 Haar values](../cptac_tumor_covariance_null_20261009/499_draws.json), independently hash-bound at `e45385676be011c6f2df7f6dbbe30202d5a28418bc6e490e8af176e77b89fcfa`; this source was not edited.
- **Primary execution:** [GitHub Actions run #38099020457](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38099020457); primary job PASSED. Its ZIP artifact `11686333294` SHA-256 `cab926478259800659afc8e74f73a01d0c5003575aa4d89150321fda485eaa02`; primary scientific calculation took **9.0579 seconds** excluding runner/package setup.
- **Independent deterministic replay:** [GitHub Actions run #38099306068](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38099306068); all **499** replay values were **numerically and order-wise identical** to the first execution, down to the recorded full-precision floating-point values. Replay ZIP SHA-256 `042569765104465a5353f448462b91e9f92479dd0040570b02cf73ea5512ab9e`. JSON metadata such as CI checkout commit and timing differs between executions.
- **Unit and numerical controls:** nine tests passed: exact single-feature multiset conservation, seeds, covariance disruption in a synthetic correlated matrix, Gram-eigen/PCA full-SVD pairwise-distance equivalence, invalid-input rejection, Monte Carlo intervals, synthetic circle/line topology, orthogonal-score-space H1 invariance, and hash-bound prior evidence. Observed-cloud maximum pairwise-distance discrepancy versus independently fitted full SVD was **1.71e-13**, and first permuted cloud **4.62e-13**, far below the `1e-8` tolerance.
- **Source data input SHA-256:** CPTAC log2-ratio CSV `d7d81d7297c0ebdc7ce7f2542e14d939f33d0853873c13e6f3bbb43cb6b238a3`; sample labels `5cc62f00952d608141e2e9e706e7f1cf271c5d443cdfcf7ddb35978bc1b77f6e`. Input checksums and `110 x 2000` post-preprocessing shape were checked in the executables.
- **Compute environment:** standard GitHub-hosted Ubuntu 24.04, Python 3.11 and pinned `requirements.txt`, single-threaded BLAS/OpenMP. No GPU, nonstandard large runner or historical 30-shard run.

For exact reproduction, from the repository root install its pinned requirements, then run:

```bash
sha256sum -c checksums.sha256
python -m unittest discover -s audit/cptac_matched_null -p 'test_*.py' -v
python audit/cptac_matched_null/experiment.py --mode confirm --output /tmp/cptac_tumor_matched_null_499.json
```

The final command is **not needed** to consume the already archived result; it explicitly reruns 499 nulls with the fixed seed `20261011`.

## Remaining uncertainties and publishing boundary

Patient case/aliquot independence, within-cohort exchangeability, acquisition batch effects, non-Gaussian feature distributions and the original mixed-population selection effects remain unverified. No independent source cohort was run under these same two nulls, and no biological mechanism was proven. The stronger manuscript conclusions about class-residualized classifier significance were already withdrawn in the [methodological correction](../../paper/METHODOLOGICAL_CORRECTIONS_2026-10-09.md) and remain withdrawn.

**Archiving:** This follow-up was executed **after** the published [GitHub v0.2.0 release](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0) and its [Zenodo DOI 10.5281/zenodo.23272317](https://doi.org/10.5281/zenodo.23272317). **Do not claim these new matched-null draws are in that immutable v0.2.0 Zenodo ZIP.** The completed study has [separate v0.3.0 release notes](../../release/v0.3.0/RELEASE_NOTES.md). Use the v0.3.0 fixed source tag for new work, and cite its *new version DOI* only after independently verifying the Zenodo deposit and file list. Independent-cohort validation is a [separate unexecuted phase](../../docs/INDEPENDENT_COHORT_VALIDATION_PHASE.md).

**Decision:** This matched sensitivity experiment has reached its stated stop condition. More null reruns solely to lower a p-value are not justified. Independent patient/batch provenance and a separately replicated cohort would be the next evidence needed for any broader biological claim.
