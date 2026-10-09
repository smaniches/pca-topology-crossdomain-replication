# CPTAC covariance-preserving null experiment — frozen protocol (2026-10-09)

**ID:** `CPTAC-TUMOR-COV-NULL-001`  
**Category:** new *post-preregistration* sensitivity, not a confirmatory re-analysis of the July 2026 program  
**Freeze point:** before viewing the outcome of any CPTAC covariance-preserving null draw  
**Base:** repository `217ab89824743ed916b2bd055a9656e4faf3b03d`

## Question and competing explanations

**Question:** Does the maximum one-dimensional Vietoris–Rips persistence of CPTAC-CCRCC *tumor-only* PCA(50) scores exceed a conditional reference that keeps **all sample-centered cross-protein covariance** fixed?

- **H1:** the ordering/configuration of individual tumor sample points carries an H1 feature unusually long compared with conditional second-order geometry.
- **H0:** a sample-centered matrix-normal model with independent, exchangeable rows and arbitrary fixed feature covariance suffices; the observed H1 value is typical of Haar-rotated sample-space matrices that preserve the sample mean and feature cross-products.

These are conditional hypotheses. Even rejecting this null does **not** establish biological loops: finite-sample heavy tails, dependencies, batch structure, and non-Gaussian distributions remain explanations.

## Data and exact scope

Use the existing cached CPTAC PDC000127 log2-ratio proteomics matrix and labels. Verify the input hashes exactly as in `audit/cptac_resid_null/experiment.py`. Reuse `confound_audit_cptac_ccrcc.preprocess` without mutation: impute zero, subtract per-sample median, select top 2,000 proteins by variance in the **mixed** original population, standardize the selected features on the **mixed** original population. **Then subset tumor samples** (110 from 194). This preserves the original exploratory feature-selection lineage but is not a tumor-only de novo feature selection. The choices may carry cross-class information, and their effect is not tested here. In particular do **not** apply label-dependent residualization or evaluate a classifier.

All reference draws use the **same** fixed feature set and standardization; they are not new independent cohorts.

## Conditional orthogonal-null construction

Let `X` be the tumor-only `n × p` standardized matrix, with `n=110`, `p=2000`, `mu` its column means, and `C=X-1 mu`. Let `U` be a fixed orthonormal `n × (n-1)` Helmert basis perpendicular to the all-ones vector; let `Z=U.T @ C`. Compute the exact reduced singular-value decomposition `Z=L diag(s) V.T` once. The observed *exact* PCA(50) scores are `U @ L[:,:50] @ diag(s[:50])`, which are distance-equivalent to the full-SVD PCA implementation.

For each draw `b`, generate a Haar-uniform `(n-1)×(n-1)` orthogonal matrix `Q_b` using QR of iid standard-normal entries with diagonal-sign correction. The surrogate tumor data are `X_b = 1 mu + U Q_b Z`. The **exact** surrogate PCA(50) scores can be computed without rebuilding the full matrix as `U Q_b L[:,:50] @ diag(s[:50])`. This is mathematically equivalent to fitting full-SVD PCA on `X_b` (up to arbitrary component sign conventions). The code MUST test this equivalence on synthetic matrices.

**Fixed invariants for *every* surrogate:**
1. Same sample count, feature count, tumor mean vector, and singular values.
2. **Exactly the same full sample-centered feature cross-product:** `C_b.T @ C_b = C.T @ C`, hence the same feature covariance and correlation wherever defined.
3. Same PCA eigenvalue spectrum, retained 50 singular values, and total retained variance.
4. Same H1 backend (`ripser` with maximum finite H1 interval) applied to observed and all null score point clouds.

The null *does not* preserve the empirical per-feature distributions, individual sample-feature values, sample-pair labels, batch labels, or rowwise biological covariates. It changes the sample-space orientation conditional on all second-order feature statistics.

## Exchangeability and validity domain

Under an iid matrix-normal model of tumor rows with a shared arbitrary feature covariance and mean, an orthogonal rotation on the sample-centered row subspace is distribution-preserving conditional on the feature scatter matrix. **This property is not generally true for real patient data.** The resulting tail area is only calibrated for this explicitly assumed row-exchangeable Gaussian/elliptical reference and the selected, globally standardized feature matrix. The source sample/case/aliquot independence and acquisition batch structure are not independently verified. The original preprocessing's mixed-population variance-based feature selection is part of the fixed observed analysis and may violate unconditional selection exchangeability; the conditional null is treated as *conditional on the selected matrix*, not a complete generative null for that selection stage.

Do **not** use the classical data-label permutation p-value to test this result; labels are used only to select the tumor subset.

## Pilot gates and prospective decision rules

**Before CPTAC null execution:**
- deterministic tests: Helmert orthogonality; centered feature covariance and mean preservation, fixed-eigenvalue spectrum; Haar-QR orthogonality and reproducibility; exact full-SVD PCA distance equivalence, including after applying a rotation; finite-value protections; synthetic Gaussian sanity case;
- input SHA-256 and class counts verified;
- numerical parity of `ripser` statistics when the same score cloud is reoriented by an **orthogonal transformation in feature-score space**;
- standard GitHub public-repository CPU runner only, bounded thread counts, no paid GPU or larger runner.

**Pilot:** 24 independent Haar draws using seed `20261009`. Pilot is for invariant checks, runtime and null-scale inspection; it is **not** the primary inferential report. If it fails any invariant, abort; do not adapt the significance rule to the observed outcome.

**Primary execution:** 499 independent draws using a distinct predeclared seed `20261010` from the frozen same implementation. Report observed max-H1, *all* unrounded null values, null mean and sample standard deviation, exceedances, and the plus-one upper-tail Monte Carlo p-value `(1 + #{T_b >= T_obs}) / 500`, with uncertainty acknowledged. The prespecified sensitivity signal criterion is one-sided `p <= 0.05`; use no retrospective threshold adjustment. It is a **post-hoc study with prospectively frozen within-study thresholds**, not preregistered confirmation. A result above 0.05 is a non-rejection, not proof of the null. A low p-value is conditional evidence only, not biological attribution.

## Resource/verification gate

Use the existing pinned Python dependencies, compute a single exact reduced SVD, and rotate 50-dimensional scores per draw. The pilot must pass and runtime be measured before confirm is permitted. Cancel if a synthetic invariant is violated, if H1 is not finite, or if cohort/hash verification fails. Store structured machine-readable outputs and GitHub run identity; commit a scientific report and raw null list without rewriting historical outputs.

**Potential falsifiers of strong original interpretation:** if observed H1 is ordinary under this covariance-preserving null, then the existing original featurewise-permutation excess can plausibly be explained by **second-order structure alone** under this null's assumptions. If unusually large, it excludes that particular conditional matrix-normal reference but does not distinguish biological topology from other deviations.

**Stop condition:** both real-versus-null computation and exact invariant checks verified, with the limited conclusion and all draws archived. Do not iterate until a preferred p-value appears.
