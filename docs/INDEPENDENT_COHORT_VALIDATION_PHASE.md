# Independent-cohort validation — separate future research phase

**Status:** NOT STARTED. No new cohort acquired, selected, analyzed, or used to generate a scientific result under this phase. This file records the next scientific question; it is **not** a continuation, outstanding requirement, or post-hoc amendment of the completed `CPTAC-TUMOR-MATCHED-NULL-002` study.

## Completed finding and boundary

The completed [110-tumor matched CPTAC comparison](../results/cptac_tumor_matched_null_20261010/REPORT.md) held samples, selected proteins, PCA(50), Euclidean distance, and max-H1 statistic fixed, and obtained `p=0.002` against independent protein permutations versus `p=0.092` against a covariance-preserving Haar reference (499 draws each). Its source and numerical result have been executed twice and committed; **do not reopen it to obtain a different p-value**.

That contrast is evidence of **null-model sensitivity**, not proof that tumor topology is biologically meaningful, that covariance alone explains it, or that the result transfers across cohorts. Both reference transformations have different invariants.

## New falsifiable research questions

1. Does a separately sourced cancer proteomics cohort, with independently ascertained tumor samples, show a similar difference between these two specified reference nulls?
2. Do sample pairing, repeated patients or aliquots, batch labels, acquisition platform, or clinical covariates account for the apparent topology or the divergence between null models?
3. Is the contrast stable under a genuinely independently selected feature panel and pre-specified sample inclusion, rather than the earlier mixed tumor/normal feature-selection lineage?

## Promotion gates before computing new outcomes

- **Independence:** establish provenance of each new cohort, patient identifiers, repeated aliquots, tumor-label definitions, measured proteins, technical batches, acquisition platform and relevant covariates. Prove it is not a restatement or subset of the original CPTAC data.
- **Preoutcome design:** separately freeze an accession, exact data and metadata checksums, exclusions, missing-data policy, normalization/feature selection procedure, component count or fixed explained-variance policy, distance metric, H1 statistic, two null transformations, number of draws, seeds, confidence intervals, and success/failure criteria. Original CPTAC numbers may motivate but **must not** be treated as independently prospectively discovered hypotheses.
- **Null assumptions:** assess exchangeability within the new patient/batch structure. Where naive permutation would break a known constraint, predeclare patient- or batch-conditioned alternatives and explicitly mark loss of strict comparability.
- **Negative controls:** verify synthetic positive and negative shapes; confirm exact preservation of declared null invariants and check PCA / persistent-homology numerical parity. Negative results stay published.
- **Statistical scope:** control for multiple tested cohorts, feature panels and transformations; report Monte Carlo resolution and uncertainty; do not infer a biological mechanism from two unequal null-model p-values.
- **Execution:** begin only after cohort availability, permissions, data integrity, budget and reproducibility gates are satisfied. Preserve all raw draws and exact-source provenance before interpreting conclusions.

## Completion criteria for this *future* phase

A genuinely independent cohort, frozen analysis contract, covariate/batch audit, complete numerical evidence, exact rerun and bounded scientific interpretation must all exist. Until then the phase status remains **UNKNOWN / NOT EXECUTED**. A successful v0.3.0 GitHub/Zenodo release **does not** satisfy any of these requirements.

## Scope and versioning

Do not modify the [July preregistration](../prereg/PREREGISTRATION.md), original published manuscripts, [v0.2.0 source snapshot](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0), completed [v0.3.0 matched study](../results/cptac_tumor_matched_null_20261010/REPORT.md), or their raw evidence. Independent validation receives a new experiment identifier and independent release only after evidence is ready.
