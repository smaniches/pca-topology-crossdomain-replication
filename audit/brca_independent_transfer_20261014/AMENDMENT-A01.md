# BRCA-TRANSFER-001 — Pre-outcome schema amendment A01

**Date:** 2026-10-14. **Status:** FROZEN BEFORE ANY REAL-DATA PCA, H1 OR NULL-REFERENCE COMPUTATION.
**Original protocol commit:** b7359b748b79908896da933b9919963e7858891e.
**Reason for amendment:** the first CI data-eligibility run failed closed at sample mapping, reporting only 3 matched columns. No scientific statistic, comparison, or reference draw ran (the confirm job was skipped). This is a source-schema correction, **not** a post-result tuning choice.

## Header-only evidence (no abundance values analyzed)

The authenticated mirror at commit 5863a2b27256555793b2e726d4a1cd33fbe01d82 has:
- 80 tumor **measurement columns** matching \`<TCGA-patient-code>.<two-digit-iTRAQ-experiment-number>TCGA\`;
- 77 unique tumor patient codes, all matched to the clinical table's \`Complete TCGA ID\`;
- exactly three patients with two measurements each: AO-A12D (runs 01/05), AO-A12B (01/34), C8-A131 (01/32);
- three normal-tissue columns ending \`.CPTAC\`, which must remain excluded.

The original protocol assumed \`.01TCGA\` was a unique primary-tumor suffix. In this mirror the two-digit suffix is an experiment/channel identifier, so that test wrongly selected only three tumor measurements. The original PDC000173 metadata independently documents that three TCGA tumor samples were measured in duplicate for QC. The original first-stage QC exception was \`Expected between 60 and 110 distinct matched tumor samples, found 3\`. Its 499-draw confirmation job was not executed.

## Replacement sample-mapping contract

1. Match all columns using \`^([A-Z0-9]{2}-[A-Z0-9]{4})\\.([0-9]{2})TCGA$\`, ignoring case, and match the patient part against canonicalized clinical TCGA IDs. Do **not** interpret the two-digit suffix as tumor type or sample quality.
2. Group measurements by canonical patient code. Expect precisely **80 tumor measurement columns mapping to 77 unique patients and three doublets** from this hash-pinned source. Fail closed if a patient has >2 measurements or the source differs from this frozen cohort structure.
3. For each protein and patient, take the arithmetic mean of the available **log2 iTRAQ abundance ratios** across that patient's measurements, ignoring non-finite replicates; if both are absent, retain a missing value to be handled by the unchanged 80%-coverage and median-imputation gates. This produces one independent observation per patient. Averaging log ratios is a fixed geometric-mean-like summarization of original abundances; it is not a blinded assessment of replicate comparability.
4. Sort the resulting **77 patient rows** lexicographically by canonical code; drop all three normal columns. Do not select one favorable technical replicate, use stage/subtype labels or include 80 non-independent observations.
5. All feature eligibility, 2000-feature selection, per-sample median centering, standardization, PCA(50), ripser/H1 statistic, two 499-draw null distributions, seeds, rejection thresholds and four-way interpretation **remain exactly as frozen in the original protocol**.

## Remaining limitations and historical preservation

This mirror does not provide authoritative iTRAQ batch/channel and replicate-quality flags. Combining two different acquisitions (including a documented pre/post-maintenance replicate) may leave residual technical heterogeneity; any result stays conditional on this fixed aggregation and cannot be elevated to a batch-adjusted scientific replication. Investigate official PDC batch and published QC selection only **after** completing the frozen comparison, in a separately preregistered follow-up, not by quietly changing current inclusion.

Preserve the failed eligibility run, original protocol and amendment in history. Never rewrite the source or fabricate an outcome to make the replication agree with the 2026-10-10 renal study. Do not merge or mint a new DOI until the investigator reviews the completed scientific evidence and limitations.
