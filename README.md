# PCA vs. persistent homology: a cross-domain replication test

**Status: WORK IN PROGRESS (pilot phase complete, cross-dataset replication in progress).**
This is a living repository — it will be updated as the replication, ablation, confound-audit,
and paper phases complete. It is NOT yet a finished, packaged deliverable and has not been run
through the project's `repo-completeness-gate` audit. See TODO below.

## What this project tests

Does PCA(50) inflate persistent-homology (H1) signal in a way that reflects genuine biological
structure, or is the inflation a generic geometric artifact of PCA that appears even on pure or
permuted noise? See [`prereg/PREREGISTRATION.md`](prereg/PREREGISTRATION.md) for the full locked
hypothesis, primary statistic, hyperparameters, and decision rules (locked 2026-07-08, before any
new-dataset analysis).

## Current state

- **Pilot analysis (GSE81089, NSCLC bulk RNA-seq, n=218):** complete. See
  [`results/pilot_GSE81089/report.md`](results/pilot_GSE81089/report.md) for the full write-up.
  Headline: real data shows genuine H1 signal beyond both a Gaussian-noise null and a pipeline-
  symmetric permutation null (z = 5.2-26.1), but PCA-driven persistence inflation on pure/permuted
  noise (+1.5 to +2.0) exceeds the PCA-driven inflation on real data (+0.78) -- meaning "PCA
  increased persistence" is not by itself evidence PCA revealed more real structure.
  **This pilot was exploratory** (hyperparameters were chosen while looking at the data) and is
  NOT one of the confirmatory tests governed by the pre-registration; it motivated the locked
  protocol used for all subsequent datasets.
- **Cross-dataset replication:** not yet started as of this commit.
- **Ablation sweep, confound audit, paper, final packaging:** not yet started.

## Repository layout

```
prereg/PREREGISTRATION.md    Locked pre-registration (hypotheses, stats, hyperparameters, thresholds)
results/pilot_GSE81089/      Pilot report, figures, and results tables (GSE81089 NSCLC cohort)
code/                        Analysis code extracted from the pilot session
```

## Known gaps (TODO before this repo can pass repo-completeness-gate)

- [ ] `code/` scripts are extracted from an interactive session and currently assume prior
      in-memory state (e.g. `nulls`, `gauss_null_raw` dicts) rather than being clean, standalone,
      run-from-scratch scripts. A single `reproduce.py` / Makefile entry point is needed.
- [ ] No `requirements.txt` / environment lockfile yet.
- [ ] No `LICENSE` / `CITATION.cff` yet.
- [ ] No automated tests or CI.
- [ ] Cross-dataset replication (2nd bulk RNA-seq cohort, TCGA multi-omics, CPTAC proteomics),
      ablation sweep, confound-attribution audit, LaTeX paper, and adversarial referee review are
      all pending -- see the project plan for the full remaining scope.

## Provenance

Generated as part of the TOPOLOGICA project (Santiago Maniches, ORCID 0009-0005-6480-1987).
