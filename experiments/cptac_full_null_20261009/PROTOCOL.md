# EXP-CPTAC-FULLNULL-20261009: complete null-distribution reproduction

**Type:** one independent re-execution of a previously reported, pre-registered cohort. Not a new pre-registration, new biological sample, or an independent scientific replication.

**Question:** Can the packaged CPTAC-CCRCC script re-execute its archived conclusions at the full pre-registered Monte Carlo counts, rather than merely reproduce the deterministic statistics using reduced null draws?

**Frozen source:** `smaniches/pca-topology-crossdomain-replication`, base commit `70f7f6db91f59839eadc66bcbc1151e03e248edc`. Original script: `code/replication_CPTAC_CCRCC/01_fetch_preprocess_and_results_table.py`. Do not modify the 2026-07-08 locked preregistration, historical results, manuscript, or data.

**Dataset:** already-committed public PDC PDC000127 log2-ratio proteomics (194 samples, 110 tumor, 84 normal). Verify SHA-256 using unchanged `checksums.sha256` before running.

| Frozen input | SHA-256 |
|---|---|
| Proteomics matrix | `d7d81d7297c0ebdc7ce7f2542e14d939f33d0853873c13e6f3bbb43cb6b238a3` |
| Labels | `5cc62f00952d608141e2e9e706e7f1cf271c5d443cdfcf7ddb35978bc1b77f6e` |
| Original analysis script | `767d32b200e9db3a4a9a9f2f038d51ac31fbbe9e93a58a569e6ea0deb9aa0700` |
| Dependency list | `f4b6f1925182c3b664b4f131e44fc55220d01405625beca1d7384f6eaf7f5957` |

**Parameters:** 2,000 HVGs; PCA(50); Gaussian null N=500; pipeline-symmetric feature-permutation null N=2,000; fixed seed 42; same GUDHI persistent homology method. Use original preprocessing; do not change metrics or thresholds to make results agree.

**Preflight:** Verify repository checksums, sample identity and exact label ordering (historical reader forcibly overwrites the label index), cohort counts, dependency installation, and two Gaussian/two permutation smoke draws with strict deterministic observed-statistic checks. Abort full computation if these invariants fail.

**Interpretation:** Report all four observed-vs-null comparison rows, their null means, standard deviations, z scores and Monte Carlo p values, and explicitly compare each against the archived CPTAC report. Retain disagreements and negative outcomes; do not treat scientifically valid hypothesis rejection as workflow failure. Assess H0 (all four z>3) and H1 (PCA null inflation >= real PCA inflation) without reinterpreting the locked hypotheses.

**Prior gate:** The original CPTAC report records PCA(50) intrinsic dimension MLE 10.22 and TwoNN 9.58, ratio <0.3; Euclidean VR was permitted. This packaged rerun independently recomputes MLE but not TwoNN and cannot establish the metric cascade boundary afresh. Do not present it as a complete independent implementation.

**Operational ceiling:** one standard, free GitHub-hosted Ubuntu runner; max 120 minutes; no GPU, paid larger runners, data downloads, parallel experiments, publication, or release. Archive full stdout, source/input hashes, versions, results, and output checksum manifest for later review. Changes must remain on an isolated research PR.
