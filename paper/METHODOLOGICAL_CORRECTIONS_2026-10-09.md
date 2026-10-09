# Methodological correction and claim boundaries — 2026-10-09

**Status:** Post-analysis correction notice. This is **not** a change to the locked 2026-07-08 preregistration, a new confirmatory result, or a replacement of recorded null draws.

**Read this before citing** `paper/paper.pdf`, `paper/sections/*.tex`, `results/final_verdict.md`, or the original manuscript's confound-independence conclusions. Those files were prepared before the issues below were established. The original PDF and numerical outputs remain in the repository as historical evidence. This notice does not update any external archival deposit; readers of a previously downloaded PDF may not see it.

## 1. The reported residualized-classifier permutation p-value is not valid

The original manuscript reports `p=0.005` for residualized classifier area under the receiver operating characteristic curve (AUC) in GSE81089, GSE146889, and TCGA-LUAD. **These significance claims are withdrawn.** The observed numerical AUCs remain records of the historical computation, but are **not** valid estimates of held-out discrimination or demonstrated confound removal.

The relevant execution sequence in `code/confound_attribution_audit/confound_audit_common.py` is:

1. `residualize_class_mean(X_std, tumor_mask)` computes a different class mean using the *true label of every sample*, including samples later held out by cross-validation.
2. `PCA.fit_transform` uses the full label-residualized matrix.
3. `cv_auc(X_pca_resid, y)` cross-validates a classifier on that already label-informed representation.
4. `auc_collapse_permutation_test` permutes `y` while keeping the representation computed from the **original** `y` fixed.

The holdout fold has already influenced its own features through its true class membership. The permutation procedure also fails to recompute the label-dependent transformation for each label assignment. Its observed statistic and null statistics therefore do not have the claimed exchangeability. A markedly below-chance AUC (for example the historical GSE81089 value `0.094`) does **not** establish that tumor/normal information was eliminated; a reversed discriminant ranking is a possible alternative. Do not reinterpret `1-AUC` as a newly validated classifier score.

The independently executed GSE146889 audit had already identified the permutation problem and withheld that p-value; see `audit/gse146889_clean_lineage/RECONCILIATION.md`. This correction extends the same code-level warning to the shared historical implementation across cohorts.

**What would be required for a predictive claim:** an evaluation design in which every operation that uses labels is fit on training observations only and applied to held-out samples **without their true labels**, with a valid full-pipeline null/refitting scheme. This correction notice does not claim such a re-evaluation has been performed.

## 2. Historical post-residualization nulls were not pipeline matched

The historical `residualization_control` in the same shared module applies class-mean subtraction to the **observed** matrix, then computes H1 on its refitted principal component analysis (PCA) projection. Its `gaussian_null_maxh1` and `permutation_null_maxh1_pca_refit` reference draws are **not** each subjected to the same label-conditioned class-mean subtraction. Historical post-residualization z-scores and p-values therefore cannot by themselves establish persistence after a matched confound-removal procedure.

This qualification applies to reported residualized-space significance across the shared cohort audit. It does **not** retroactively recalculate, erase, or invalidate separately generated preregistered raw/PCA real-versus-null numbers, within-class observed H1 magnitudes, or independently reported sample-bootstrap values. Those have their own assumptions and must be assessed on their own definitions.

The independent CPTAC-CCRCC **post-preregistration** sensitivity analysis at `results/cptac_resid_null_sensitivity_20261009/REPORT.md` corrected that processing asymmetry for one specific null: independently shuffle each feature **within** each class, then class-center, refit PCA(50), and compute maximum H1 using the same implementation on observed and surrogate matrices. Its 499-draw empirical one-sided p-value was `0.002` (0 exceedances). **This is not a replacement p-value for any old AUC or Gaussian-null claim.** The surrogate destroys cross-protein covariance and may break sample-level dependencies, so it does not rule out ordinary second-order covariance, technical batch, or other nonbiological structure.

## 3. TCGA-LUAD methylation does not demonstrate loss caused by PCA itself

In `code/replication_TCGA_LUAD/02_methylation_fetch_preprocess_and_results_table.py`, the methylation cohort uses 36 observations and the requested PCA dimension is capped at 35. For 36 centered samples, retaining all 35 sample-centered components preserves all **Euclidean** pairwise distances, up to floating-point error, so Euclidean Vietoris–Rips persistence cannot be lost by this full-rank PCA transformation.

The manuscript's reported zero maximum H1 is instead obtained after computing **a spectral graph-Laplacian distance** on the PCA coordinates, as implemented in `code/replication_TCGA_LUAD/_common.py`. That spectral metric was selected via the protocol's intrinsic-dimension gate. The **metric-gated failure under that specified analysis** remains the reported numerical outcome, but language claiming that PCA itself destroyed the methylation loop is unsupported. Moreover, comparing raw Euclidean persistence with PCA-space spectral persistence mixes a representation change with a distance-metric change; it is not a causal estimate of PCA's isolated effect.

See `results/replication_TCGA_LUAD/tcga_luad_report.md`, including its Euclidean PCA(35) transparency comparison, and `prereg/PREREGISTRATION.md` Section 5 for the metric gate.

## 4. What remains established, what remains open

- **Historical execution records retained:** the original raw/PCA persistence calculations, registered feature/metric choices, numerical null draws where archived, follow-up GSE146889 clean-lineage execution, and full-count CPTAC re-execution remain available. Passing verification means those computations match their declared inputs and controls; it is not proof of biological mechanisms.
- **Cohort dependence remains descriptive:** within-class and within-stratum results vary by cohort and can inform follow-up; they do not alone certify independence from every relevant confound.
- **No biological-topology mechanism has been established:** the principal original Gaussian and independently shuffled-feature null families remove cross-feature covariance. Their rejection is compatible with ordinary covariance structure. The next discriminating sensitivity is a covariance-preserving null with a frozen contract and validated invariants, ideally on a single-class subset to avoid label-conditioned residualization and patient-pairing complications.
- **Original manuscript is historical:** `paper/paper.pdf` and `paper/sections/*.tex` contain wording predating this notice. Their original numbers are not silently rewritten here. The original preregistration and recorded results are unchanged. Any future revised manuscript must carry a version identifier and distinguish those revisions from the original archival PDF.

## 5. Subsequent covariance-preserving sensitivity (run after this correction was frozen)

The follow-up [CPTAC tumor-only covariance-preserving experiment](../results/cptac_tumor_covariance_null_20261009/REPORT.md) used 110 tumor samples, fixed original mixed-cohort feature selection, exact PCA(50) scores, and 499 Haar-rotated sample-space nulls preserving all tumor-only cross-protein covariance. The observed max-H1 was 5.3543 and 45/499 null values were at least as large (plus-one p=0.092), so the declared one-sided 0.05 threshold was **not** met. All 499 draw values and a deterministic independent replay are archived. This is an **outcome appended after the correction notice**, not a change to the pre-outcome correction rationale.

It must not be contrasted as a causal experiment against the earlier p=0.002 from a **mixed tumor/normal class-residualized** independent-feature null; both the sample subset and null differ. The new result is a bounded conditional non-rejection, not a claim that biological structure is absent.

## Code and audit trail

- Class-mean residualization, class-conditioned AUC and historical null code: [`code/confound_attribution_audit/confound_audit_common.py`](../code/confound_attribution_audit/confound_audit_common.py)
- Existing independently identified GSE146889 AUC null problem: [`audit/gse146889_clean_lineage/RECONCILIATION.md`](../audit/gse146889_clean_lineage/RECONCILIATION.md)
- Methylation implementation: [`code/replication_TCGA_LUAD/02_methylation_fetch_preprocess_and_results_table.py`](../code/replication_TCGA_LUAD/02_methylation_fetch_preprocess_and_results_table.py) and [`_common.py`](../code/replication_TCGA_LUAD/_common.py)
- CPTAC matched-transform experiment, source, and raw 499 draws: [`results/cptac_resid_null_sensitivity_20261009/REPORT.md`](../results/cptac_resid_null_sensitivity_20261009/REPORT.md)
- Original frozen protocol: [`prereg/PREREGISTRATION.md`](../prereg/PREREGISTRATION.md)

**Decision consequence:** cite the original manuscript together with this correction notice. Do not reuse the invalid AUC p-values as support for biological confound-independence, and do not attribute the methylation spectral collapse to PCA alone.
