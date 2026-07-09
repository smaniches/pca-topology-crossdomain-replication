# Cross-dataset replication figure (Figure 2)

`04c_figure_cross_dataset_replication.py` regenerates
`paper/figures/fig2_cross_dataset_replication.png` from data already
committed elsewhere in this repository:

- Panel (a) (H0 z-scores) and panel (b) (H1 PCA-deltas): `results/cross_dataset_BH_family.csv`,
  `results/pilot_GSE81089/final_results_table.csv`,
  `results/replication_GSE146889/final_results_table_GSE146889.csv`,
  `results/replication_CPTAC_CCRCC/cptac_h1_pca_delta_table.csv`,
  `results/replication_TCGA_LUAD/tcga_luad_h1_pca_delta_summary.csv`.
- Panel (c), CV AUC: `results/confound_attribution_audit/table3_residualization_control[_<COHORT>].csv`
  and `results/replication_TCGA_LUAD/tcga_luad_confound_cv_auc.csv`.

Recovered via `host.lineage` on the session that originally produced this
figure (previously it existed only as a static PNG with no committed
generation script); verified to reproduce the committed figure's data/layout
exactly.

**Disclosed gap:** panel (c)'s dominant-H1-loop tumor-enrichment numbers
(samples touched / total, Fisher's exact p) are NOT present in any committed
CSV -- they exist only as prose in each cohort's own `results/` markdown
report:
- GSE81089: `results/pilot_GSE81089/report.md` ("touches only 3 of 218 samples ... p=1.0")
- GSE146889: `results/replication_GSE146889/report_GSE146889.md` ("45 of 176 samples ... p = 2.87e-09")
- CPTAC-CCRCC: `results/replication_CPTAC_CCRCC/cptac_ccrcc_report.md` (PCA(50) row: 4 loop vertices / 194, p=0.135)
- TCGA-LUAD: `results/replication_TCGA_LUAD/tcga_luad_report.md` ("112/116 samples ... p = 0.1185")

The script transcribes these four (n, N, p) triples as literal constants
rather than computing them, because recovering them as structured,
re-derivable data would require re-running each cohort's cocycle trace
against its raw persistence diagrams, and those raw diagrams are not all
retained in the repository (see the manuscript's Data/Code Availability
section for which intermediates are and are not checkpointed). This is
disclosed here rather than left silent.
