# Repository Manifest

Full file inventory for smaniches/pca-topology-crossdomain-replication.
SHA-256 checksums for every file are in `checksums.sha256` (verify with `sha256sum -c checksums.sha256` from the repo root).

## (root)/

| File | Size |
|---|---|
| CITATION.cff | 662 B |
| LICENSE | 1.3 KB |
| MANIFEST.md | 6.1 KB |
| README.md | 7.9 KB |
| reproduce.py | 4.4 KB |
| requirements.txt | 145 B |

## code/

| File | Size |
|---|---|
| 03_statistics_and_results_table.py | 10.0 KB |
| 04a_figure_persistence_diagrams.py | 5.7 KB |
| 04b_figure_barcode_null_distributions.py | 7.9 KB |

## code/ablation_sweep/

| File | Size |
|---|---|
| 00_download_and_preprocess.py | 4.1 KB |
| 01_make_configs.py | 2.8 KB |
| 02_run_sweep.py | 11.2 KB |
| 03_aggregate_results.py | 11.3 KB |
| README.md | 5.3 KB |
| verify_default.py | 2.7 KB |

## code/confound_attribution_audit/

| File | Size |
|---|---|
| README.md | 8.3 KB |
| confound_audit_common.py | 14.8 KB |
| confound_audit_cptac_ccrcc.py | 9.5 KB |
| confound_audit_gse146889.py | 14.3 KB |
| confound_audit_gse81089.py | 8.7 KB |
| confound_audit_tcga_luad.py | 13.2 KB |
| tcga_luad_file_manifest.json | 19.0 KB |

## code/replication_CPTAC_CCRCC/

| File | Size |
|---|---|
| 01_fetch_preprocess_and_results_table.py | 14.6 KB |
| README.md | 2.9 KB |

## code/replication_CPTAC_CCRCC/data/

| File | Size |
|---|---|
| cptac_ccrcc_labels.csv | 5.8 KB |
| cptac_ccrcc_log2ratio_raw.csv | 12.43 MB |

## code/replication_GSE146889/

| File | Size |
|---|---|
| 01_fetch_preprocess_and_results_table.py | 14.0 KB |
| README.md | 2.8 KB |

## code/replication_TCGA_LUAD/

| File | Size |
|---|---|
| 01_rnaseq_fetch_preprocess_and_results_table.py | 14.9 KB |
| 02_methylation_fetch_preprocess_and_results_table.py | 18.1 KB |
| README.md | 5.3 KB |
| _common.py | 5.6 KB |

## paper/

| File | Size |
|---|---|
| .gitignore | 79 B |
| paper.pdf | 1.53 MB |
| paper.tex | 1.7 KB |
| references.bib | 3.7 KB |

## paper/figures/

| File | Size |
|---|---|
| fig1_pilot_motivation.png | 246.8 KB |
| fig2_cross_dataset_replication.png | 311.0 KB |
| fig3_ablation_sensitivity.png | 394.1 KB |
| fig4_confound_audit.png | 267.1 KB |

## paper/sections/

| File | Size |
|---|---|
| abstract.tex | 4.0 KB |
| discussion.tex | 11.9 KB |
| introduction.tex | 6.7 KB |
| methods.tex | 12.1 KB |
| results.tex | 23.3 KB |

## prereg/

| File | Size |
|---|---|
| PREREGISTRATION.md | 7.1 KB |

## results/

| File | Size |
|---|---|
| cross_dataset_BH_family.csv | 2.1 KB |
| cross_dataset_checkpoint.md | 11.4 KB |
| final_verdict.md | 16.5 KB |

## results/ablation_sweep/

| File | Size |
|---|---|
| ablation_report.md | 20.7 KB |
| ablation_sweep_full_table.csv | 7.2 KB |
| imputation_ablation_summary.csv | 666 B |
| imputation_overlap_table.csv | 246 B |
| interaction_check_table.csv | 631 B |
| permutation_convergence_table.csv | 1.1 KB |

## results/ablation_sweep/figures/

| File | Size |
|---|---|
| interaction_and_imputation_overlap.png | 121.9 KB |
| sweep_sensitivity.png | 232.8 KB |

## results/confound_attribution_audit/

| File | Size |
|---|---|
| confound_audit_CPTAC_CCRCC.md | 18.0 KB |
| confound_audit_GSE146889.md | 15.4 KB |
| confound_audit_GSE81089.md | 13.8 KB |
| confound_audit_TCGA_LUAD.md | 19.3 KB |
| cross_dataset_confound_audit_summary.csv | 1.3 KB |
| table1_within_class_decomposition.csv | 423 B |
| table1_within_class_decomposition_CPTAC_CCRCC.csv | 723 B |
| table1_within_class_decomposition_GSE146889.csv | 427 B |
| table1_within_class_decomposition_TCGA_LUAD.csv | 419 B |
| table2_within_stratum_control.csv | 485 B |
| table2_within_stratum_control_CPTAC_CCRCC.csv | 491 B |
| table2_within_stratum_control_GSE146889.csv | 877 B |
| table2_within_stratum_control_TCGA_LUAD.csv | 607 B |
| table3_residualization_control.csv | 238 B |
| table3_residualization_control_CPTAC_CCRCC.csv | 373 B |
| table3_residualization_control_GSE146889.csv | 476 B |
| table3_residualization_control_TCGA_LUAD.csv | 483 B |
| table4_block_bootstrap_ci.csv | 326 B |
| table4_block_bootstrap_ci_CPTAC_CCRCC.csv | 243 B |
| table4_block_bootstrap_ci_GSE146889.csv | 520 B |
| table4_block_bootstrap_ci_TCGA_LUAD.csv | 325 B |
| table5_confound_spectrum.csv | 241 B |
| table5_confound_spectrum_CPTAC_CCRCC.csv | 264 B |
| table5_confound_spectrum_GSE146889.csv | 3.5 KB |
| table5_confound_spectrum_TCGA_LUAD.csv | 2.4 KB |

## results/confound_attribution_audit/figures/

| File | Size |
|---|---|
| bootstrap_ci_and_auc_check.png | 115.0 KB |
| bootstrap_ci_and_auc_check_CPTAC_CCRCC.png | 122.5 KB |
| bootstrap_ci_and_auc_check_GSE146889.png | 80.2 KB |
| bootstrap_ci_and_auc_check_TCGA_LUAD.png | 155.4 KB |
| within_class_null_comparison.png | 103.0 KB |
| within_class_null_comparison_CPTAC_CCRCC.png | 112.9 KB |
| within_class_null_comparison_GSE146889.png | 61.1 KB |
| within_class_null_comparison_TCGA_LUAD.png | 119.7 KB |
| within_stratum_and_residualization.png | 103.7 KB |
| within_stratum_and_residualization_CPTAC_CCRCC.png | 118.9 KB |
| within_stratum_and_residualization_GSE146889.png | 88.9 KB |
| within_stratum_and_residualization_TCGA_LUAD.png | 115.1 KB |

## results/pilot_GSE81089/

| File | Size |
|---|---|
| final_results_table.csv | 2.1 KB |
| permutation_test_summary.csv | 893 B |
| proposed_registry_update.md | 3.7 KB |
| report.md | 14.5 KB |

## results/pilot_GSE81089/checkpoints/

| File | Size |
|---|---|
| diagrams.pkl | 57.8 KB |
| gaussian_null_dist.pkl | 23.9 KB |
| null_distributions.pkl | 188.2 KB |

## results/pilot_GSE81089/figures/

| File | Size |
|---|---|
| barcode_and_null_distributions.png | 202.8 KB |
| persistence_diagrams_real_vs_gaussian.png | 160.1 KB |

## results/replication_CPTAC_CCRCC/

| File | Size |
|---|---|
| bh_correction_family_table.csv | 669 B |
| cptac_ccrcc_report.md | 12.8 KB |
| cptac_final_results_table.csv | 756 B |
| cptac_h1_pca_delta_table.csv | 340 B |

## results/replication_CPTAC_CCRCC/figures/

| File | Size |
|---|---|
| null_distributions.png | 105.4 KB |
| pca_delta_and_persistence_diagram.png | 141.6 KB |

## results/replication_GSE146889/

| File | Size |
|---|---|
| final_results_table_GSE146889.csv | 1.1 KB |
| report_GSE146889.md | 14.1 KB |

## results/replication_GSE146889/figures/

| File | Size |
|---|---|
| barcode_and_pca_delta_GSE146889.png | 94.3 KB |
| persistence_diagrams_GSE146889.png | 120.7 KB |

## results/replication_TCGA_LUAD/

| File | Size |
|---|---|
| tcga_luad_confound_cv_auc.csv | 196 B |
| tcga_luad_h1_pca_delta_summary.csv | 674 B |
| tcga_luad_primary_stats_summary.csv | 1.2 KB |
| tcga_luad_report.md | 12.6 KB |

## results/replication_TCGA_LUAD/figures/

| File | Size |
|---|---|
| tcga_luad_h1_pca_delta.png | 77.8 KB |
| tcga_luad_null_comparison.png | 167.3 KB |
