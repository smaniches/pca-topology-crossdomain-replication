# BRCA-TRANSFER-001: frozen-reference transfer to an independently sourced breast cohort

**Status:** EXECUTED ONCE; exact null-draw array and source-integrity checks reviewed; awaiting independent second-run confirmation and original-study batch/aliquot audit. **Draft research, not a new published scientific result or third-party replication.**

**Scientific investigator:** Santiago Maniches, TOPOLOGICA LLC. AI-assisted code and internal review, not unaffiliated replication.  
**Execution:** [GitHub Actions run 38614749781](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38614749781), 2026-10-14.  
**Exact evidence:** [QC and 998 full null draws, workflow artifact 11918044893](https://github.com/smaniches/pca-topology-crossdomain-replication/actions/runs/38614749781/artifacts/11918044893), stored for 90 days by GitHub Actions; ZIP SHA-256 `bdd97e720ad5ba3c9b9d5a020bc4a7575eb2227bd0be346b565db11554535d60`. JSON source-file SHA-256 `c6514a6278378d2867759ecb2f7385101dc07421778fd00709015d79872a3a6d`; QC JSON SHA-256 `0b2442a133e378b19110cb6bcab42a7acb991bb0f5e7bac07505bc21fc2233dd`.

## Question and status of design

The prior renal tumor experiment `CPTAC-TUMOR-MATCHED-NULL-002` was published in [v0.3.0/Zenodo 23290337](https://doi.org/10.5281/zenodo.23290337); its observed H1 was unusual under independent-feature permutation but not under a covariance-preserving Haar reference, conditional on each reference model.

Does the **same direction of null-model sensitivity** appear in an independently collected TCGA breast-cancer cohort?

The original frozen BRCA protocol was committed as `b7359b748b79908896da933b9919963e7858891e` before any BRCA PCA/H1 outcome. The first QC run correctly rejected a mistaken `.01TCGA` patient-code assumption **without executing any real-data topology/null computation**. The pre-outcome, schema-only [amendment A01](AMENDMENT-A01.md), committed as `19eefd7bf3af1bff1cd206fae33531e7b5e8eda7`, explicitly reinterprets the two-digit suffix as the iTRAQ experiment and aggregates three duplicated patients. No sample/feature selection, null, seed, threshold or inference rule changed after results became visible.

## Observed cohort and numerical outcome

Data: *Mertins et al.* (2016), *Nature*, [doi:10.1038/nature18003](https://doi.org/10.1038/nature18003), TCGA breast cancer iTRAQ tumor proteomics, original study [PDC000173](https://proteomic.datacommons.cancer.gov/pdc/study/PDC000173). Pinned publicly distributed mirror at `DominicJSykes/BreastCancerProteomics` commit `5863a2b27256555793b2e726d4a1cd33fbe01d82`. The analyzed source comprises **80 tumor measurement columns representing 77 distinct patients** (three duplicate pairs averaged in log2 abundance space), plus **three excluded normal columns**. The 77-by-2000 feature matrix uses preregistered missingness, median imputation, per-sample centering, variance selection and feature scaling. Its exact floating-point SHA-256 is `9b3fb769a50fe78de3ce0cea00739f729307ec18f62a4a1db6661255ba32cd1d`.

- Observed longest finite Euclidean Vietoris–Rips H1 interval after PCA(50): **4.963970184326172**.
- Both references use the same **77 patients, 2000 features, observed preprocessing, PCA(50) and H1 statistic**, 499 new draws apiece, independently seeded as specified.

| Reference (conditional) | BRCA-TRANSFER-001 breast cohort | Previously published renal cohort |
|---|---:|---:|
| Independently permute each protein's 77 observed values | **7/499** null intervals ≥ observed; `p=(1+7)/500=0.016` | **0/499**, `p=0.002` |
| Covariance-preserving sample-space Haar rotations | **74/499** null intervals ≥ observed; `p=(1+74)/500=0.150` | **45/499**, `p=0.092` |

The feature-permutation null distribution has mean **3.4236486278220504**, sample SD **0.602682317997741**, and exact binomial Monte Carlo 95% interval for its *unadjusted exceedance fraction* **[0.005658113650, 0.028688565738]**. The Haar reference has mean **4.257772342475478**, SD **0.7183417507997607**, and analogous exceedance-fraction interval **[0.118273217165, 0.182558103274]**. These are **simulation uncertainty intervals**, not confidence intervals for a biological effect or for the plus-one p-estimator.

**Prespecified descriptive pattern:** The featurewise null rejects at threshold `p<=0.05`, whereas the full-covariance-preserving Haar reference does not. **The same directional pattern appears in both distinct cohorts**, but this does not mean the p-values, underlying biological signal or mechanism replicated quantitatively.

## Executed verification and provenance

- The read-only external-data fetch verified pinned Git blob identity **before** analysis. Proteomics CSV SHA-256 `d8464512d296b07636d8760876b79f34e098b604fd164f4532ef4cdbb6048df2`; clinical CSV SHA-256 `243e0aa23197ffc68ce6198c89c939312438a86bc55db32a79d6ad90b4162350`.
- New protocol synthetic method/failure-gate tests all passed before the real-data workflow executed.
- Observed PCA score-cloud pairwise-distance parity against separately fitted full SVD: max absolute difference **1.56319e-13**. First independent-feature null: **3.12639e-13**; first Haar null: **2.13163e-13**.
- First-feature-null empirical column multisets were identical; Haar first-draw cross-protein covariance relative Frobenius deviation **5.72016e-16**, below frozen `1e-9` tolerance.
- On a separate local analysis of the archived **998 full-precision draw values** (read-only, without rerunning the experimental generator), the 7 and 74 exceedances, both plus-one p-values, sample means/SDs, matrix dimensions, and missingness/QC counts were recomputed and matched. Full ordered-draw JSON array SHA-256 values (compact decimal arrays): feature arm `f3c3d92a5dd810d3f79179da7d5eb49e0a23890cc730761282e536e145a09ea2`, Haar arm `854c7a714796ccf64dbd3678b9c57b349963f416069fbb240e6465e5c5bf19f9`. This verifies **arithmetic consistency**, not independent experimental replication.
- The source has 12,553 protein rows; 10,078 pass the frozen 80%-finite eligibility filter; 13,610 residual missing protein-by-patient values across eligible proteins were median-imputed. The exact selected-protein accession-list SHA-256 is `faf025e51d7c89d8019ce748b51e0c9f4d4905e6d860652ec1103c828ca2705b`. Three clinical-matched technical replicate patient groups were averaged, not counted as independent subjects.

## Limits: what this does not show

1. **Batch/channel/aliquot confounding is not addressed.** The source files lack a canonical, independently reconciled full instrument batch, maintenance, labeling-channel, and replicate-quality mapping; the original [PDC study](https://proteomic.datacommons.cancer.gov/pdc/study/PDC000173) documents common iTRAQ references, three duplicated tumors and a pre/post-maintenance acquisition. Both conditional reference models could react to technical structure or cohort heterogeneity. The source is a pinned **secondary public mirror**, not raw canonical PDC acquisition evidence.
2. **This is not an unaffiliated scientific replication.** The same investigator-directed and AI-assisted workflow constructed and interpreted the new test. The subjects, tumor type and study are distinct; the analysis is computationally executed internally.
3. **The two references retain different invariants.** Column permutations retain empirical protein marginal distributions but destroy cross-protein alignment; Haar rotations retain sample-centered protein covariance and PCA spectrum while generally changing marginals and specimen identity distributions. Their tail-probability difference **cannot isolate covariance as a causal mechanism** or certify that biological topological loops exist.
4. **Outcome-informed design:** The original renal p-values were known before BRCA-TRANSFER-001 was designed. Pre-outcome freezing for breast computations was honored, but this is a targeted transfer check, **not a blind preregistered novel discovery**.
5. **No production- or clinical-use implication** and no publication of altered previous Zenodo/GitHub v0.3.0 scientific records. External batch-aware validation and independent third-party replication remain open.

## Promotion and repeatability

The current run completed all executable gates. **Next gate:** a second exact deterministic run on the frozen source/code, independent post-run digest comparison and original PDC batch/aliquot QC reconciliation. After these, external scientific review is still needed for any assertion beyond a finite two-cohort null-reference sensitivity observation.

**Keep PR as a draft; no scientific publication or new DOI from this run without the investigator's review.**