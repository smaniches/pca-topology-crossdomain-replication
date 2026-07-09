# Confound-attribution audit (4 cohorts)

Reproduces `results/confound_attribution_audit/confound_audit_{GSE81089,CPTAC_CCRCC,
TCGA_LUAD,GSE146889}.md` and their supporting `tableN_*.csv` files.

## What this tests

The core TOPOLOGICA finding is that omics point clouds show statistically significant
max-H1-persistence (a topological loop) beyond matched Gaussian and permutation nulls.
But every cohort's tumor/normal classes are also separable at near-ceiling AUC by a
simple classifier -- so does the H1 signal reflect genuine within-tumor topology, or is
it an artifact of the tumor/normal class-mean-shift confound? This audit answers that
with 4 controls, applied identically across all 4 cohorts:

1. **Within-class decomposition** -- recompute max-H1 for the tumor-only subset alone
   (no normal samples, so no class-mean-shift direction exists to exploit). If tumor-only
   reproduces the mixed-set magnitude, the signal is not confound-driven.
2. **Within-stratum control** -- bin samples by a continuous confound score (the linear
   class-mean-shift projection) into quartiles; check whether the z-score against a
   matched null survives within each narrow-confound-range stratum.
3. **Residualization** -- explicitly remove the class-mean-shift direction (shift each
   sample by its class mean minus the pooled grand mean), refit PCA, recompute max-H1 and
   classifier AUC; check whether signal and AUC collapse together (would indicate the
   signal WAS confound-driven) or the signal survives while AUC collapses (genuine).
4. **Block-bootstrap CI** -- resample the tumor-only subset with replacement, get a 95%
   CI on (observed - null mean); check whether it excludes 0.

## Files

- `confound_audit_common.py` -- shared statistical machinery (all 4 controls + the
  max-H1-persistence statistic under two backends, `ripser` and `gudhi`). Cohort-agnostic:
  takes a preprocessed data matrix + tumor/normal mask + null-generation parameters.
- `confound_audit_gse81089.py` -- GSE81089 pilot cohort (218 NSCLC samples). Reproduces
  the pre-registered default pipeline (top-2000 HVG, zero-fill imputation, PCA(50),
  ripser). **Verified: mixed and tumor-only max-H1 match the committed report exactly
  (3.887 both).**
- `confound_audit_cptac_ccrcc.py` -- CPTAC-CCRCC cohort (194 samples, proteomics via the
  PDC GraphQL API). Median-centered, top-2000-variance-protein preprocessing; **gudhi**
  backend (matches what was actually run for this cohort -- confirmed by tracing the
  sub-agent transcript's `run_audit.py` launch to its printed real-data output, which
  matches the committed report to 4 significant figures). **Verified: mixed and
  tumor-only max-H1 match exactly (3.833 both).**
- `confound_audit_tcga_luad.py` -- TCGA-LUAD cohort (116 paired samples via the GDC API).
  **This is the cohort where mixed != tumor-only is the reported finding** (mixed=8.327,
  tumor-only=4.596), i.e. this cohort's mixed-set signal IS partly confound-driven,
  unlike the other 3. Ships `tcga_luad_file_manifest.json`, a pinned list of the exact
  116 GDC file_ids used in the original run -- 6 of 58 paired cases have more than one
  STAR-Counts file in the current GDC index with no documented tie-break rule, so a
  fresh dynamic API query is not guaranteed to reproduce the exact 116-sample cohort;
  the pinned manifest is the default and gives an exact match to the committed report.
  **Verified: mixed=8.327, tumor-only=4.596, residualized (PCA-space) z_gaussian~=0.65
  (report: ~0.50, near-null collapse) using the pinned manifest.**
- `confound_audit_gse146889.py` -- GSE146889 cohort (176 samples). **Reconstructed from a
  delegated sub-agent's conversation transcript, not from code lineage** (see limitation
  below). Reuses the cohort's Phase-1 replication run's mixed-set null distributions
  (fetched by artifact id) rather than regenerating them, matching what the original run
  did. **Verified: mixed=7.564 and tumor-only=5.880 match the committed report closely
  (near-exact); residualized z_gaussian in the right range at reduced draw counts.**

## Data provenance

- GSE81089: re-fetched from GEO FTP (`GSE81089_FPKM_cufflinks.tsv.gz`) from scratch.
- CPTAC-CCRCC: re-fetched from the CPTAC PDC GraphQL API (study PDC000127) from scratch.
- TCGA-LUAD: re-fetched from the GDC API using a pinned file-id manifest (see above);
  optionally re-derivable dynamically with `--data-dir ... ` and no manifest, at the cost
  of exact reproducibility for 6/116 samples.
- GSE146889: re-fetched from GEO FTP (`GSE146889_GeneCount.tsv.gz`) from scratch, but the
  mixed-set null distributions are reused from a prior run's artifacts rather than
  regenerated (by design, matching the original audit).

## Running

```bash
python confound_audit_gse81089.py --n-gauss 500 --n-perm 2000 --n-boot 2000       # full (~hours)
python confound_audit_cptac_ccrcc.py --n-gauss 300 --n-perm 300 --n-boot 200      # full (~hours)
python confound_audit_tcga_luad.py --n-gauss 500 --n-perm 2000 --n-boot 1000      # full (~hours)
python confound_audit_gse146889.py --n-gauss 300 --n-perm 300 --n-boot 2000       # full (~hours)

# Quick correctness smoke test for any cohort (real-data statistics are exact regardless
# of draw count; only null-derived z/p-values become approximate):
python confound_audit_gse81089.py --quick 20
```

## Verification performed

Ran all 4 drivers locally with reduced draw counts (`--quick 15`-`20`, disclosed
explicitly in each run's output and here) due to compute-time constraints on a 4-CPU
sandbox. Verified real-data max-H1 statistics (deterministic given data + preprocessing +
PCA seed, independent of null draw count) against each cohort's committed report:

| Cohort | Mixed (expected) | Mixed (computed) | Tumor-only (expected) | Tumor-only (computed) |
|---|---|---|---|---|
| GSE81089 | 3.887 | 3.887 (exact) | 3.887 | 3.887 (exact) |
| CPTAC-CCRCC | 3.833 | 3.833 (exact) | 3.833 | 3.833 (exact) |
| TCGA-LUAD | 8.327 | 8.3265 | 4.596 | 4.5960 |
| GSE146889 | 7.564 | 7.5643 | 5.880 | 5.880 |

A full-precision run (pre-registered draw counts: 500 Gaussian / 2000 permutation /
1000-2000 bootstrap per cohort) was NOT executed in this session due to compute-time
constraints (would take many hours across 4 cohorts) -- this is the disclosed reduction.
The deterministic real-data statistics above do not depend on this choice; the z-scores
and p-values reported by each script at reduced draw counts are directionally correct
(all match the reported pattern -- e.g. TCGA-LUAD residualized z_gaussian collapsed
toward the reported near-null ~0.50) but not to full numerical precision.

## Known limitations

- **GSE146889 reconstruction is transcript-based, not lineage-based.** This cohort's
  confound audit was originally run by a delegated sub-agent whose code is not reachable
  through this project's code-lineage mechanism (the report artifact's lineage graph has
  only 1 node -- itself). This script was reconstructed by manually reading that
  sub-agent's ~95-cell conversation transcript and reassembling the pipeline in
  dependency order. This is inherently more error-prone than the lineage-based extraction
  used for the other 3 cohorts: some implementation details (exact RNG seed sequencing in
  a few null-generation branches, the precise per-quartile permutation-null construction)
  were inferred from partial transcript evidence rather than pulled as exact producing
  code. The script's real-data statistics were verified to match the committed report
  closely (mixed and tumor-only max-H1 match to 3-4 significant figures), which gives
  reasonable confidence the reconstruction is correct in substance, but this cohort's
  script carries materially more reconstruction risk than the other 3 and should be
  reviewed line-by-line before being treated as an authoritative replication artifact.
- **TCGA-LUAD file selection for 6/58 paired cases is not uniquely determined by the GDC
  API alone** (see above) -- this script ships a pinned manifest to guarantee exact
  reproduction; a fresh dynamic fetch (no manifest) will very likely select a different
  file for those 6 cases and give slightly different real-data statistics.
- Draw counts were reduced for local testing in this session (see table above); the
  4 scripts default to the pre-registered draw counts and should be run with those
  defaults (no `--quick` flag) for a publication-grade full-precision replication.
