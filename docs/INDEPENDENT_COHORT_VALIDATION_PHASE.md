# Independent-cohort validation — separate future research phase

**Status update, 2026-10-14:** a distinct, externally sourced but internally executed conditional transfer check, `BRCA-TRANSFER-001`, now has two complete deterministic runs on a pinned public TCGA breast-cancer proteomics cohort (77 patients; both ordered 499-draw null arrays identical across runs). See its [results and limits](../audit/brca_independent_transfer_20261014/RESULTS.md). The cohort is separately sourced, but the check was designed, executed and interpreted by the same investigator-directed, AI-assisted workflow; it is **not** a replication by independent external researchers. Its original batch/channel/aliquot reconciliation is outstanding, its source is a secondary public mirror, and its design was informed by the known renal outcome. It therefore does not meet the promotion gates and completion criteria below. The **broader batch-reconciled, external independent-validation phase remains NOT COMPLETED.** BRCA-TRANSFER-001 is not part of the v0.3.0 release or any Zenodo version.

**Status as recorded in the v0.3.0 release snapshot:** NOT STARTED. No new cohort acquired, selected, analyzed, or used to generate a scientific result under this phase. That statement was accurate for v0.3.0 and is retained unchanged in its archive. This file records the next scientific question; it is **not** a continuation, outstanding requirement, or post-hoc amendment of the completed `CPTAC-TUMOR-MATCHED-NULL-002` study.

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

A genuinely independent cohort, frozen analysis contract, covariate/batch audit, complete numerical evidence, exact rerun and bounded scientific interpretation must all exist. Until then the phase status remains **UNKNOWN / NOT EXECUTED**; the narrower BRCA-TRANSFER-001 check dated above does not change that. A successful v0.3.0 GitHub/Zenodo release **does not** satisfy any of these requirements.

## Scope and versioning

Do not modify the [July preregistration](../prereg/PREREGISTRATION.md), original published manuscripts, [v0.2.0 source snapshot](https://github.com/smaniches/pca-topology-crossdomain-replication/releases/tag/v0.2.0), completed [v0.3.0 matched study](../results/cptac_tumor_matched_null_20261010/REPORT.md), or their raw evidence. Independent validation receives a new experiment identifier and independent release only after evidence is ready.
