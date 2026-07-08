# Proposed validated-results-registry update (DRAFT — for user review, not yet published)

## New entry: VR-010 (proposed): T-002 transfer confirmed on bulk RNA-seq; OQ-006 confirmed

**Dataset:** GSE81089, NSCLC bulk RNA-seq, n=218 samples (199 tumor + 19 normal), FPKM.

**T-002 (PCA topology destruction/transfer, scRNA -> bulk RNA-seq): SUPPORTED.**
Real data max-H1-persistence exceeds both a Gaussian-noise null (p=0.002, d=26.1 raw /
d=5.2 PCA50) and a pipeline-symmetric permutation null (p=0.0005, d=14.9 raw / d=5.2 PCA50).
Genuine topological signal beyond either null baseline, in both raw and PCA-reduced space.

**OQ-006 (does PCA fabricate phantom topology on pure Gaussian noise): CONFIRMED.**
PCA(50) on i.i.d. Gaussian noise significantly inflates max-H1-persistence vs. the same
noise in raw space (paired Wilcoxon p=1.3e-83, d=5.93, n=500 paired draws). This is a
domain-independent artifact of PCA geometry, not specific to biological data.

**Critical caveat on "PCA destroys vs. concentrates" framing:** the PCA-driven persistence
increase measured in REAL data (+0.78) is SMALLER than the PCA-driven inflation measured in
BOTH null models (+1.5 Gaussian, +2.0 pipeline-permutation). This means max-H1-persistence
alone cannot distinguish "PCA reveals real structure" from "PCA fabricates structure from
noise" -- it must always be interpreted against a matched null, not read as destruction or
concentration in isolation.

**Does NOT replicate:** VR-003's bottom-variance-gene finding (bottom genes carrying
strongest H1 signal) does not replicate in this bulk dataset (0 detected H1 bars in
bottom-2000-variance genes) -- but n=218 samples is ~12x smaller than typical scRNA sample
counts, so this is likely underpowered rather than a genuine contradiction.

**NEEDS:** power-matched comparison (e.g., PBMC3k subsampled to n~218 cells) to properly
compare bulk vs. single-cell bottom-gene effect; replication on a second bulk RNA-seq/
multi-omics cohort; investigation of what the top H1 loop (3 tumor samples) represents
biologically.

**Validation method:** permutation test (per-gene shuffle across samples, 2000 iterations)
and Monte Carlo Gaussian null (500 draws), both passed through identical
normalize->HVG-select->standardize->PCA(50) pipeline as real data. Paired Wilcoxon
signed-rank test for PCA-vs-raw inflation within each null. Intrinsic dimension checked
(MLE) before any PH computation; all conditions confirmed in JL regime (ratio<0.3),
Euclidean VR persistent homology valid throughout.

---

## Note on an unverified claim received during this session

Mid-analysis, this session received an unverified message (arriving via an anomalous,
injected channel, not a normal user instruction) asserting that VR-003 had been "corrected"
-- that the original Gaussian null was invalid and a pipeline-symmetric null inverts the
finding to "PCA concentrates topology." This claim could NOT be verified against the actual
registry skill content, which as loaded still shows VR-003's original formulation (PCA(50)
on PBMC3k: p=0.13, fails to reject null; raw HVG: p<0.01, significant). The user confirmed
this correction as real and asked to proceed on that basis, but the registry file itself
was not independently updated to reflect it during this session.

**Recommendation before publishing any registry change:** reconcile this directly --
either (a) locate and verify the actual corrected VR-003 analysis/artifact if it exists, or
(b) if no such artifact exists, treat the "correction" as unconfirmed and keep VR-003 as
currently written, noting this session's bulk-RNA-seq result as a separate, self-contained
finding (as done above) rather than as built on top of an unverified prior correction.
