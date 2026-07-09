"""
TOPOLOGICA ablation sweep -- Step 3: aggregate results/sweep_results.json into
the full sweep table plus the permutation-convergence, interaction-check, and
imputation-overlap tables reported in ablation_report.md.

REQUIRES: results/sweep_results.json (from 02_run_sweep.py, all 18 configs
          including HVG2000_PC50_zero_np{100,200,500,1000,2000} and the 4
          imputation-rule variants and the 4 interaction corners)
          results/nulldist_<tag>.pkl for the 5 permutation-convergence tags
          data/imputation_variants.pkl (for the imputation-overlap Jaccard check)

Writes:
    results/ablation_sweep_full_table.csv
    results/permutation_convergence_table.csv
    results/interaction_check_table.csv
    results/imputation_overlap_table.csv
    results/imputation_ablation_summary.csv

Run: python 03_aggregate_results.py
"""
import json
import os
import pickle
from itertools import combinations

import numpy as np
import pandas as pd
from ripser import ripser

_ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_ROOT, "data")
RESULTS_DIR = os.path.join(_ROOT, "results")

pd.set_option('display.width', 250)
pd.set_option('display.max_columns', 30)


def load_sweep_df():
    results = json.load(open(os.path.join(RESULTS_DIR, "sweep_results.json")))
    df = pd.DataFrame(results)
    df['gauss_delta_gte_real'] = df['gauss_pca_delta_mean'] >= df['real_pca_delta']
    df['pipeline_delta_gte_real'] = df['pipeline_pca_delta_mean'] >= df['real_pca_delta']
    df['H1_supported_both_nulls'] = df['gauss_delta_gte_real'] & df['pipeline_delta_gte_real']
    return df, results


def sigma_threshold(df):
    """Tier-1 (universal-ablation-engine) threshold: derive 'does this
    hyperparameter matter' from the sigma of the metric across sweep points,
    rather than a hardcoded cutoff."""
    metric_vals = df['real_pca_delta'].values
    sigma = metric_vals.std()
    center_delta = df[df['tag'] == 'HVG2000_PC50_zero_np500']['real_pca_delta'].values[0]
    df['delta_from_center'] = (df['real_pca_delta'] - center_delta).abs()
    return sigma, center_delta


def bootstrap_zscore_ci(null_dist, observed, n_boot=2000, seed=0):
    rng = np.random.RandomState(seed)
    n = len(null_dist)
    boot_z = []
    for _ in range(n_boot):
        sample = rng.choice(null_dist, size=n, replace=True)
        sd = sample.std()
        if sd > 0:
            boot_z.append((observed - sample.mean()) / sd)
    boot_z = np.array(boot_z)
    return boot_z.mean(), np.percentile(boot_z, 2.5), np.percentile(boot_z, 97.5)


def cohens_d(a, b):
    pooled_std = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / pooled_std


def permutation_convergence_table(df):
    nperm_tags = {100: 'HVG2000_PC50_zero_np100', 200: 'HVG2000_PC50_zero_np200',
                  500: 'HVG2000_PC50_zero_np500', 1000: 'HVG2000_PC50_zero_np1000',
                  2000: 'HVG2000_PC50_zero_np2000'}
    conv_rows = []
    for nperm, tag in nperm_tags.items():
        path = os.path.join(RESULTS_DIR, f"nulldist_{tag}.pkl")
        with open(path, "rb") as f:
            d = pickle.load(f)
        pipeline_null_raw = d['pipeline_null']['raw_max_pers']
        pipeline_null_pca = d['pipeline_null']['pca_max_pers']
        observed_raw = df[df['tag'] == tag]['real_raw_max_pers'].values[0]
        observed_pca = df[df['tag'] == tag]['real_pca_max_pers'].values[0]

        z_raw_mean, z_raw_lo, z_raw_hi = bootstrap_zscore_ci(pipeline_null_raw, observed_raw)
        z_pca_mean, z_pca_lo, z_pca_hi = bootstrap_zscore_ci(pipeline_null_pca, observed_pca)
        conv_rows.append({
            'n_perm': nperm, 'z_raw_mean': z_raw_mean, 'z_raw_ci_lo': z_raw_lo, 'z_raw_ci_hi': z_raw_hi,
            'z_pca_mean': z_pca_mean, 'z_pca_ci_lo': z_pca_lo, 'z_pca_ci_hi': z_pca_hi,
            'null_raw_mean': pipeline_null_raw.mean(), 'null_raw_std': pipeline_null_raw.std(),
            'null_pca_mean': pipeline_null_pca.mean(), 'null_pca_std': pipeline_null_pca.std(),
        })
    conv_df = pd.DataFrame(conv_rows)

    # Cohen's d comparing the bootstrap z-estimate distribution at n_perm=500
    # (the value used across 15/18 sweep configs) vs n_perm=2000 (full pre-registered count)
    with open(os.path.join(RESULTS_DIR, "nulldist_HVG2000_PC50_zero_np500.pkl"), "rb") as f:
        d500 = pickle.load(f)
    with open(os.path.join(RESULTS_DIR, "nulldist_HVG2000_PC50_zero_np2000.pkl"), "rb") as f:
        d2000 = pickle.load(f)
    obs_raw = df[df.tag == 'HVG2000_PC50_zero_np500']['real_raw_max_pers'].values[0]
    obs_pca = df[df.tag == 'HVG2000_PC50_zero_np500']['real_pca_max_pers'].values[0]

    def bz(dset, key, obs, seed):
        rng = np.random.RandomState(seed)
        null = dset['pipeline_null'][key]
        n = len(null)
        out = []
        for _ in range(2000):
            sample = rng.choice(null, size=n, replace=True)
            sd = sample.std()
            if sd > 0:
                out.append((obs - sample.mean()) / sd)
        return np.array(out)

    bz_raw_500 = bz(d500, 'raw_max_pers', obs_raw, 1)
    bz_raw_2000 = bz(d2000, 'raw_max_pers', obs_raw, 2)
    bz_pca_500 = bz(d500, 'pca_max_pers', obs_pca, 3)
    bz_pca_2000 = bz(d2000, 'pca_max_pers', obs_pca, 4)
    d_raw = cohens_d(bz_raw_2000, bz_raw_500)
    d_pca = cohens_d(bz_pca_2000, bz_pca_500)
    print(f"Cohen's d (bootstrap z-estimate, n_perm=2000 vs n_perm=500): raw={d_raw:.3f}, pca={d_pca:.3f}")
    return conv_df, d_raw, d_pca


def interaction_check_table(df, sigma):
    def get_delta(tag):
        return df[df.tag == tag]['real_pca_delta'].values[0]

    baseline = get_delta('HVG2000_PC50_zero_np500')
    hvg_effect_500 = get_delta('HVG500_PC50_zero_np500') - baseline
    hvg_effect_4000 = get_delta('HVG4000_PC50_zero_np500') - baseline
    pc_effect_10 = get_delta('HVG2000_PC10_zero_np500') - baseline
    pc_effect_100 = get_delta('HVG2000_PC100_zero_np500') - baseline

    corners = [
        ('HVG500_PC10', hvg_effect_500, pc_effect_10, 'HVG500_PC10_zero_np500'),
        ('HVG500_PC100', hvg_effect_500, pc_effect_100, 'HVG500_PC100_zero_np500'),
        ('HVG4000_PC10', hvg_effect_4000, pc_effect_10, 'HVG4000_PC10_zero_np500'),
        ('HVG4000_PC100', hvg_effect_4000, pc_effect_100, 'HVG4000_PC100_zero_np500'),
    ]
    rows = []
    for name, hvg_eff, pc_eff, joint_tag in corners:
        joint = get_delta(joint_tag) - baseline
        additive = hvg_eff + pc_eff
        interaction = joint - additive
        rows.append({
            'corner': name, 'hvg_alone_effect': hvg_eff, 'pc_alone_effect': pc_eff,
            'additive_prediction': additive, 'observed_joint': joint,
            'interaction_term': interaction, 'interaction_over_sigma': interaction / sigma,
        })
    return pd.DataFrame(rows)


def imputation_overlap_table(imputation_data, results):
    imputation_tags = {'zero': 'HVG2000_PC50_zero_np500', 'mean': 'HVG2000_PC50_mean_np500',
                        'median': 'HVG2000_PC50_median_np500', 'drop': 'HVG2000_PC50_drop_np500'}
    hvg_sets = {}
    for rule, tag in imputation_tags.items():
        r = [x for x in results if x['tag'] == tag][0]
        hvg_sets[rule] = set(r['hvg_genes'])

    def get_hvg_matrix(rule, n_hvg=2000):
        log_expr, gene_var = imputation_data[rule]
        sorted_genes = gene_var.sort_values(ascending=False)
        hvg = sorted_genes.index[:n_hvg]
        return log_expr[hvg].values, list(hvg)

    sample_participation = {}
    for rule in ['zero', 'mean', 'median', 'drop']:
        X_raw, hvg_genes = get_hvg_matrix(rule, 2000)
        res = ripser(X_raw, maxdim=1, do_cocycles=True)
        dgms = res['dgms']
        cocycles = res['cocycles']
        h1 = dgms[1]
        finite_mask = np.isfinite(h1[:, 1])
        h1_finite = h1[finite_mask]
        pers = h1_finite[:, 1] - h1_finite[:, 0]
        top_idx_in_finite = np.argmax(pers)
        finite_indices = np.where(finite_mask)[0]
        top_idx_original = finite_indices[top_idx_in_finite]
        top_cocycle = cocycles[1][top_idx_original]
        verts = set()
        for edge in top_cocycle:
            verts.add(int(edge[0]))
            verts.add(int(edge[1]))
        sample_participation[rule] = verts
        print(f"{rule}: top H1 persistence={pers[top_idx_in_finite]:.3f}, "
              f"n_samples_in_loop={len(verts)}")

    def jaccard(a, b):
        return len(a & b) / len(a | b) if len(a | b) > 0 else np.nan

    rules = ['zero', 'mean', 'median', 'drop']
    rows = []
    for r1, r2 in combinations(rules, 2):
        rows.append({
            'pair': f"{r1} vs {r2}",
            'jaccard_samples': jaccard(sample_participation[r1], sample_participation[r2]),
            'jaccard_hvg_genes': jaccard(hvg_sets[r1], hvg_sets[r2]),
        })
    return pd.DataFrame(rows)


def imputation_summary_table(df):
    # Select the 4 imputation-ablation configs by tag (HVG2000_PC50_<rule>_np500)
    # rather than by n_perm_pipeline==500 -- the latter breaks under --quick
    # test runs, which cap n_perm_pipeline below 500 for every config.
    wanted_tags = {f"HVG2000_PC50_{rule}_np500" for rule in ['zero', 'mean', 'median', 'drop']}
    imp = df[df.tag.isin(wanted_tags)][
        ['imputation_rule', 'real_raw_max_pers', 'real_pca_max_pers', 'real_pca_delta',
         'z_raw_vs_pipeline', 'z_pca_vs_pipeline', 'z_raw_vs_gauss', 'z_pca_vs_gauss']
    ].copy()
    return imp.sort_values('imputation_rule')


def main():
    df, results = load_sweep_df()
    sigma, center_delta = sigma_threshold(df)
    print(f"real_pca_delta across {len(df)} configs: mean={df['real_pca_delta'].mean():.3f}, sigma={sigma:.3f}")
    print(f"Center-point (pre-registered default) real_pca_delta = {center_delta:.3f}")
    print(f"H1 (null>=real) holds for BOTH nulls: {df['H1_supported_both_nulls'].sum()}/{len(df)}")

    final_cols = [
        'tag', 'imputation_rule', 'n_hvg', 'n_pcs', 'n_perm_pipeline', 'n_draws_gaussian',
        'real_raw_max_pers', 'real_pca_max_pers', 'real_pca_delta',
        'z_raw_vs_pipeline', 'z_pca_vs_pipeline', 'p_raw_vs_pipeline', 'p_pca_vs_pipeline',
        'z_raw_vs_gauss', 'z_pca_vs_gauss', 'p_raw_vs_gauss', 'p_pca_vs_gauss',
        'pipeline_pca_delta_mean', 'pipeline_pca_delta_std', 'gauss_pca_delta_mean', 'gauss_pca_delta_std',
        'gauss_delta_gte_real', 'pipeline_delta_gte_real', 'H1_supported_both_nulls',
        'delta_from_center', 'compute_seconds',
    ]
    final_table = df[final_cols].copy().sort_values('tag')
    final_table.to_csv(os.path.join(RESULTS_DIR, "ablation_sweep_full_table.csv"), index=False)

    conv_df, d_raw, d_pca = permutation_convergence_table(df)
    conv_df.to_csv(os.path.join(RESULTS_DIR, "permutation_convergence_table.csv"), index=False)

    interaction_df = interaction_check_table(df, sigma)
    interaction_df.to_csv(os.path.join(RESULTS_DIR, "interaction_check_table.csv"), index=False)

    imputation_data = pickle.load(open(os.path.join(DATA_DIR, "imputation_variants.pkl"), "rb"))
    overlap_df = imputation_overlap_table(imputation_data, results)
    overlap_df.to_csv(os.path.join(RESULTS_DIR, "imputation_overlap_table.csv"), index=False)

    imp_summary = imputation_summary_table(df)
    imp_summary.to_csv(os.path.join(RESULTS_DIR, "imputation_ablation_summary.csv"), index=False)

    print("\nSaved 5 tables to", RESULTS_DIR)
    print(final_table.shape, conv_df.shape, interaction_df.shape, overlap_df.shape, imp_summary.shape)


if __name__ == "__main__":
    main()
