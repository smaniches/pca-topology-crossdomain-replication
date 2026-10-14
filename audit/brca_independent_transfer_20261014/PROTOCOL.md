# BRCA-TRANSFER-001 — pre-outcome independent-cohort transfer protocol

**Status:** PRE-OUTCOME SPECIFICATION. No new breast proteome null-test result was available to this study when this file was committed.  
**Date:** 2026-10-14. **Scientific investigator:** Santiago Maniches, TOPOLOGICA LLC.  
**Implementation assistance:** AI tooling under the investigator's direction; this is not independent external peer review.

## Separate question and immutable historical boundary

Does the **directional reference-model sensitivity** seen in the already published renal cancer comparison (CPTAC-TUMOR-MATCHED-NULL-002, 110 CCRCC tumors, 2000 proteins, PCA(50): feature permutation p=0.002, covariance-preserving Haar p=0.092) also occur in an **independently collected, breast-cancer proteomics cohort**?

This experiment is a new research branch, not a rerun, amendment, completion obligation or replacement of GitHub/Zenodo v0.3.0. It must not alter original preregistrations, null draws, manuscripts, Git tags, deposits, or claims. Results in either direction will be retained.

**Transfer outcome:** independently report both one-sided conditional Monte Carlo tail probabilities and whether the same descriptive threshold pattern (independent-feature null p<=0.05 and full-covariance Haar null p>0.05) appears. This is NOT a test proving covariance causality, equivalence, biologically meaningful loops, or population-wide transport of an effect. Non-rejection of the Haar reference is not evidence for absence of topology.

## Frozen independent data acquisition

- New cohort: 2016 TCGA breast cancer tumor iTRAQ proteome; primary data provenance **PDC000173**, **PXD015901**, Mertins et al., Nature 2016, DOI 10.1038/nature18003. No patient or tumor matrix is reused from the original renal PDC000127 analysis.
- Download-only redistribution mirror: `DominicJSykes/BreastCancerProteomics`, commit `5863a2b27256555793b2e726d4a1cd33fbe01d82` on GitHub. Clinical file: `clinical_data_breast_cancer.csv`, expected Git blob SHA-1 `33fc811ca6f664f101def2ac8861bb6df1962bc4`. Proteome file: `77_cancer_proteomes_CPTAC_itraq.csv`, expected Git blob SHA-1 `7c66480a7691c3e112d0d4bbd7c33cad11524f15`. The mirror is **secondary**, not a substitute for verifying official PDC sample mapping, batch information and public data license when making biological generalization claims.
- Download exactly these immutable files (git blob verification **before** any analysis) and record their SHA-256 in the run manifest. Use files ephemerally, not committed or redistributed in this project.
- Select columns whose headers encode TCGA primary tumor aliquots as `<patient-code>.01TCGA` and whose patient codes correspond to a unique `Complete TCGA ID` in the external clinical file. Exclude the three normal/reference breast tissue columns and all metadata fields. **Fail closed** on duplicate patient identifiers, ambiguous labels, fewer than 60 matched distinct tumors, nonnumeric abundance data, or unresolvable sample IDs; do not silently substitute a different cohort.
- Include all qualifying tumor samples. Do not use subtype, stage, outcomes, or survival to pick samples or features. Report counts and missingness. The original publication reported 105 tumors and a 77-sample high-quality analytical subset; the 77-sample mirror is used, but the exact recovered row count is a *data-quality gate*, not assumed.

## Deterministic preprocessing and topology statistic

1. Read log2 iTRAQ protein-ratio data. Match and sort tumors by canonical TCGA patient ID. Treat each unique RefSeq accession as one feature and keep the first occurrence in source row order if duplicates occur (report how many were suppressed).
2. Remove features with less than 80% finite measurements across included tumors. Fill remaining missing values with the **observed within-feature median**. Reject features with non-finite or zero across-tumor variance.
3. Median-center **each tumor sample** across its measured/imputed proteins. Rank features by descending across-tumor population variance on this centered matrix, breaking exact ties by accession then source row; select **exactly 2000**; fail if insufficient. Standardize each retained feature using its across-tumor mean and sample standard deviation. No tumor class labels enter preprocessing.
4. Work on the same selected and standardized n-by-2000 matrix for the observed statistic and both reference generators. Use exactly **50 PCA components** (require n>50), centering columns. Compute PCA scores by deterministic, sorted sample-Gram eigen-decomposition; confirm observed and first-reference score-space pairwise distances against independent full SVD at tolerance 1e-7.
5. Distance: Euclidean on PCA(50) scores. Primary statistic T: **longest finite H1 persistent-homology interval** from ripser Vietoris–Rips, with no scalar normalization or downstream hyperparameter tuning.

## Prespecified conditional null references

- **A, independent-feature marginal null**: 499 draws. Independently and uniformly permute the n tumor values in every one of the 2000 protein columns of the frozen analysis matrix, refit PCA(50), compute T. RNG `numpy.random.default_rng(2026101401)`. Exact column multisets must be invariant; cross-protein covariance need not be.
- **B, covariance-preserving Haar null**: 499 draws. Rotate the mean-centered n-by-2000 data in the (n-1)-dimensional sample subspace orthogonal to the all-ones vector by a random Haar orthogonal matrix (Gaussian QR and sign normalization), then add feature means back. Refit PCA(50), compute T. RNG `numpy.random.default_rng(2026101402)`. Assert conservation of column means, cross-protein covariance and PCA singular-value spectrum (relative Frobenius tolerance 1e-9 on the first draw).
- Strictly the same observed matrix, 2000 features, PCA components, Euclidean metric and T in both arms. **Do not resample or reselect features inside either null**. Conditional reference models have different invariants; no causal attribution from their p-value difference.
- Compute tail estimates `(1 + count(T_null >= T_observed)) / 500` separately, exact binomial Clopper–Pearson 95% intervals for exceedance fractions, all 499 ordered null draws, seeds, library versions, source hashes, runtime, and numerical invariant checks. Report p<=0.05 *descriptively* for each named reference; this two-reference pattern is not a new familywise-controlled discovery test.

## Nonnegotiable validation and stop conditions

Before viewing real-data T/null outcomes: commit this protocol and the exact source hashes and algorithm. Executable synthetic tests must cover no-loop collinear control, ring H1 positive control, independent column-permutation marginal invariants, orthogonal covariance preservation, patient ID uniqueness, basic malformed-input rejection and PCA distance parity. CI first downloads and authenticates exact files, then checks sample metadata; a nonmatching source is a source-contract failure, never a scientific null finding.

**Hard scientific limitation:** the data mirror does not alone supply authoritative iTRAQ analytical-batch labels, replicate resolution, instrument channel and acquisition covariates. Thus even a perfectly reproduced directional pattern is *transfer evidence in this data representation*, not batch-adjusted replication, confirmed biological topology or independent third-party scientific replication. Resolve the original batch/aliquot mapping from PDC before promoting any biological/general conclusion. Document sample IDs and source identity without distributing protected raw records.

**Outcomes remain public either way:** report the four combinations of null rejection/nonrejection, nonfinite/error cases, missing input, and any runtime/capability failure. If data acquisition, sample identity, feature sufficiency or batch-aware interpretation fails, mark the experiment `BLOCKED` or `PARTIAL`, do not tune criteria or invent substitute labels/data. Publication, Zenodo deposit and promotion to published scientific conclusion require subsequent owner review.
