# CPTAC-TUMOR-MATCHED-NULL-002 — frozen post-hoc comparison protocol

**Frozen:** 2026-10-10, before any new tumor-only feature-permutation null draw is computed.  
**Original source:** `main` commit `2dc7e1bab46c1f5d083253d90becc95efabc881f`.  
**Classification:** post-preregistration, outcome-informed follow-up. The existing covariance-preserving result (`p=0.092`) was already known when this comparison was designed. **This is not a fresh confirmatory trial.**

## Question

If the cohort, samples, preprocessing, selected protein features, principal component analysis (PCA) dimension, Euclidean distance and max-H1 persistent-homology statistic remain identical, does the *same observed tumor-only cloud* exceed (A) a reference that destroys cross-protein dependence but retains each protein's empirical values, while failing to exceed (B) the previously executed reference that preserves the complete sample-centered covariance matrix?

This adjudicates two specified *null-reference models on a matched observed analysis*, **not** whether covariance is the only cause of biological structure. The references have inherently different preserved invariants and distributional assumptions.

## Fully matched observed data and statistic

- Real data: the committed CPTAC clear-cell renal carcinoma (CCRCC) PDC000127 log2-ratio proteomics matrix and tumor/normal metadata. Check source SHA-256 exactly: `d7d81d7297c0ebdc7ce7f2542e14d939f33d0853873c13e6f3bbb43cb6b238a3` for the CSV matrix and `5cc62f00952d608141e2e9e706e7f1cf271c5d443cdfcf7ddb35978bc1b77f6e` for labels.
- Reuse `audit/cptac_covariance_null/experiment.py::load_tumor()` unchanged: original mixed 194-case zero-imputation, **per-sample** median-centering, mixed-population top-2,000 variable proteins, mixed-population `StandardScaler` fitted before selecting 110 tumor samples. This is deliberate historical lineage; no independent tumor-only feature selection is claimed.
- Input point cloud shape: **110 tumor samples × 2,000 selected standardized proteins**; no normal observations in either reference, no label-conditioned class residualization, no CV AUC analysis.
- Observed scores: exact full-rank singular-value decomposition PCA with 50 components applied **after tumor-only restriction**, then Euclidean Vietoris–Rips filtration computed with `ripser==0.6.14`; statistic `T` is the longest finite dimension-one persistence interval. The observed value must agree to within `2e-5` with the archived `5.354297637939453`. Both controls compare to this **same single observed value**.
- Fixed source covariance result: `results/cptac_tumor_covariance_null_20261009/499_draws.json` (source-file SHA-256 `e45385676be011c6f2df7f6dbbe30202d5a28418bc6e490e8af176e77b89fcfa`), containing 499 Haar rotations with seed `20261010`, 45 exceedances, one-sided plus-one `p=0.092`. Do **not** rerun or change these values solely to obtain a preferred finding.

## New reference A: tumor-only independent-feature permutations

For each draw, start with the **identical 110×2,000 preprocessed tumor-only matrix `X`**. For each protein column `j`, independently and uniformly permute its values over the same 110 tumor sample positions using a fresh `rng.permutation(n)`; keep the same 2,000 columns and their order.

- This reference retains **every column's exact 110-value empirical multiset** (hence its mean/variance and the original fixed feature selection), but generally destroys pairwise cross-protein covariance and higher-order cross-feature dependence.
- There is **no second standardization or sample-median normalization** after tumor restriction in either branch. Permutation operates only at the analysis stage after the original preprocessing.
- **Refit PCA(50) from scratch on every shuffled-feature draw.** Efficient implementation may eigendecompose the `110×110` centered sample Gram matrix; its top-component PCA score-cloud pairwise Euclidean distances must be independently compared against `sklearn.decomposition.PCA(svd_solver="full")` on both synthetic matrices and the first real tumor-only permutation, tolerance `1e-8`.
- Apply the **same ripser max-finite-H1 code to the new PCA(50) scores**, with no metric changes, smoothing, rescaling, prior label knowledge or class-based adjustments. The only change relative to the existing Haar reference is how a null point cloud is generated.

## New prospective-within-follow-up execution rules

- Synthetic/unit tests (prior to any real-data null outcome): preservation of every per-column multiset; fixed size/features; reproducible independent feature permutations; nonconservation of off-diagonal sample-centered covariance in a nondegenerate synthetic example; numerical PCA score-distance equivalence to independently computed exact SVD; invariance of H1 under score-space orthogonal rotation; positive H1 on synthetic ring and zero on collinear synthetic points.
- **Pilot:** `B=24` feature-permutation draws, seed `20261012`. Run only after tests pass. The pilot determines feasibility and catches coding bugs; it is **not** pooled with the primary null and should not change the analysis or its thresholds.
- **Primary:** `B=499` independent draws, seed `20261011`. Use the exact same source and code as pilot. One-sided Monte Carlo estimate `p=(1 + #{T_b≥T_observed})/500`. Frozen rejection criterion `p≤0.05` for *this* feature-permutation reference; descriptive z-score `(T_obs-mean(null))/SD(null)` reported with full precision and all 499 draws.
- Covariance-preserving comparison uses its archived `B=499`, not a pilot, with identical plus-one p-value arithmetic. State empirical exceedance counts and binomial Monte Carlo uncertainty for both; do not claim that merely comparing two p-values tests a causal mechanism.
- **Decision table**, labels decided before the new 499-draw result:
  - Feature-null rejects, Haar null does not: a multivariate-dependence-sensitive H1 excess survives an independent-feature null, **but no significant excess beyond the covariance-preserving Haar reference was established**.
  - Both reject: each null is rejected conditionally; additional batch, independent-cohort, and biological attribution controls are still needed.
  - Neither rejects: no excess detected against either at this threshold; non-rejection is not equivalence.
  - Haar rejects but feature-null does not: different reference families yield opposite support; investigate assumptions rather than choose a preferred p-value.
- The earlier 194-case class-residualized `p=0.002` is **not part of this matched comparison**. It changes the observed sample set and statistic.

## Verification and resource gates

- Fix Python 3.11, `requirements.txt` pins, standard GitHub-hosted Ubuntu CPU, `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1` and finite numerical checks.
- Source data, prior Haar JSON, and historical release SHA-256 verified before new draws.
- Confirm PCA mathematical and independent numerical parity; recompute the observed H1 and ensure it is equal in both null rows.
- Bound each job, avoid new paid instances/GPUs and old multi-shard workflows. If processing exceeds runner cap, stop and report—not quietly reduce the 499-draw primary contract.
- Before interpretation, store machine-readable JSON with **all** 499 unrounded null draws, source hashes, seeds, backend, environment versions, observed value, tail p, count, runtime, GitHub run ID/SHA. Independently replay the full seed once if budget permits and compare every null value in order.
- Publish a plain-language report with uncertainty, explicit assumptions and a data-derived comparison figure; preserve the historical July preregistration, `v0.2.0` Zenodo archive, and all previous negative evidence.

## Statistical scope

Both reference families assume that tumor rows are exchangeable to the extent required by their specified transformations. That is **not independently established** for CPTAC patient identity, sample pairing, acquisition batch, or proteomic measurement processes. Both analyses also condition on the previously mixed-population feature panel and standardization. The matched comparison tests observed sensitivity to *these two choices of null*, **not** an exclusive mechanistic claim that covariance explains the tumor biology.

**Stop condition:** full-count data and checks pass, outcomes reported with bounded language; do not iterate seeds until a desired significance is obtained.
