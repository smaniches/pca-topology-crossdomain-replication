# PCA vs. persistent homology: a cross-domain replication test

**Status: Manuscript complete, adversarially reviewed, and repo-completeness-gate audited.**
Pilot, cross-dataset replication (3 confirmatory cohorts), ablation sweep, confound-attribution
audit (generalized to all 4 cohorts), and a LaTeX manuscript (`paper/paper.pdf`) are all complete.
The manuscript has been through two rounds of a 5-referee adversarial review pass; see `paper/`
for the current draft and `results/` for every underlying numeric result. The repository has
passed a 7-gate completeness audit (file inventory, standard files, no hardcoded paths, end-to-end
reproduction, number provenance, portability, checksum registry) for the pilot-cohort code path;
see `MANIFEST.md` / `checksums.sha256` and the "Known gaps" section below for the one remaining,
explicitly disclosed scope limitation.

## What this project tests

Does PCA(50) inflate persistent-homology (H1) signal in a way that reflects genuine biological
structure, or is the inflation a generic geometric artifact of PCA that appears even on pure or
permuted noise? A second, related question: when a genuine H1 signal is found, is it a proxy for
tumor/normal class separability, or a geometrically distinct signal? See
[`prereg/PREREGISTRATION.md`](prereg/PREREGISTRATION.md) for the full locked hypothesis, primary
statistic, hyperparameters, and decision rules (locked 2026-07-08, before any new-dataset analysis).

## Current state

- **Pilot analysis (GSE81089, NSCLC bulk RNA-seq, n=218):** complete, exploratory (motivated the
  locked protocol; not one of the confirmatory tests). See `results/pilot_GSE81089/report.md`.
- **Cross-dataset replication (3 confirmatory cohorts, pre-registered):** complete. GSE146889
  (colorectal/endometrial/ovarian bulk RNA-seq), CPTAC-CCRCC (renal proteomics), TCGA-LUAD (lung
  adenocarcinoma RNA-seq and methylation). 14/16 BH-FDR-corrected confirmatory tests pass; see
  `results/cross_dataset_checkpoint.md` and `results/final_verdict.md`.
- **Ablation sweep (18 configurations, pilot cohort):** complete. $H_1^{\text{stat}}$ holds
  robustly at the pre-registered operating point but fails at very low PC-count and at one
  high-gene-count/high-PC-count corner; see `results/ablation_sweep/`.
- **Confound-attribution audit (all 4 cohorts):** complete. Confound-independence is
  dataset-dependent: genuine and tumor-scoped in GSE81089 and CPTAC-CCRCC, partial/weaker in
  GSE146889, and failing specifically in the pre-registered PCA(50) space for TCGA-LUAD; see
  `results/confound_attribution_audit/` and `results/final_verdict.md` §5-6.
- **Manuscript (`paper/paper.tex` / `paper/paper.pdf`):** complete, compiles cleanly (0 errors,
  0 undefined references), through 2 rounds of 5-referee adversarial review.

## Repository layout

```
prereg/PREREGISTRATION.md         Locked pre-registration (hypotheses, stats, hyperparameters, thresholds)
results/pilot_GSE81089/           Pilot report, figures, results tables (GSE81089 NSCLC cohort)
results/replication_GSE146889/    Confirmatory replication #1 (colorectal/endometrial/ovarian)
results/replication_CPTAC_CCRCC/  Confirmatory replication #2 (renal proteomics)
results/replication_TCGA_LUAD/    Confirmatory replication #3 (RNA-seq + methylation)
results/ablation_sweep/           18-configuration hyperparameter-sensitivity sweep
results/confound_attribution_audit/  Confound-attribution audit, all 4 cohorts
results/cross_dataset_checkpoint.md  Phase-1 replication synthesis (BH-FDR family)
results/final_verdict.md          Full-program synthesis: what is and isn't established
code/                             Statistics/figure-generation scripts (pilot cohort only -- see below)
paper/                            LaTeX manuscript source, figures, and compiled PDF
reproduce.py                      Single entry point: regenerates the pilot figures (and,
                                  with --full, the full pilot statistics) from scratch
requirements.txt                  Pinned package versions used throughout this project
LICENSE                           MIT license (code only; see the file for scope)
CITATION.cff                      Citation metadata (CFF 1.2.0)
MANIFEST.md                       Full file inventory with sizes, by directory
checksums.sha256                  SHA-256 for every tracked file (verify: `sha256sum -c checksums.sha256`)
```

## Reproducing the pilot cohort

```bash
pip install -r requirements.txt
python3 reproduce.py           # fast: regenerate the 2 pilot figures from bundled
                                # checkpoints (results/pilot_GSE81089/checkpoints/*.pkl), ~1 min
python3 reproduce.py --full    # slow: also re-fetch GSE81089 from GEO and re-run the
                                # full statistics (2000+500 permutation/Gaussian-null
                                # persistent-homology computations), ~20-40 min
```

Both code paths have been verified to run end-to-end from a clean copy of this repository
(exit code 0, output files present, and --full's reproduced max-H1-persistence statistic
verified against the previously-reported value to 1e-3 tolerance).

## Known gaps (disclosed, not silently omitted)

- [ ] `code/` currently contains only the pilot cohort's scripts. Analysis code for the 3
      confirmatory replications, the 18-config ablation sweep, and the 4-cohort
      confound-attribution audit is not yet packaged here (the numeric results and full
      methodology are documented in the corresponding `results/` markdown reports, but the
      code itself has not been added to this repository). This is the one remaining gap
      after the repo-completeness-gate audit -- the gate applies fully to the pilot-cohort
      code path (which passes all 7 gates) but the audit does not extend coverage claims to
      code that does not exist in the repository.
- [ ] No automated tests or CI (beyond the manual gate-audit reproduction check above).

## Provenance

Generated as part of the TOPOLOGICA project (Santiago Maniches, ORCID 0009-0005-6480-1987).
