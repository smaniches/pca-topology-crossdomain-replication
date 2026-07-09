import numpy as np
import pandas as pd
import pickle
import warnings
warnings.filterwarnings("ignore")

# This script reproduces the GSE81089 pilot analysis end-to-end from the
# public GEO accession (no local checkpoints required): fetch raw data,
# preprocess, compute persistent homology, run both null models, and
# write final_results_table.csv.

import urllib.request, os

url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/suppl/GSE81089_FPKM_cufflinks.tsv.gz"
out = "GSE81089_FPKM.tsv.gz"
urllib.request.urlretrieve(url, out)

import pandas as pd
df = pd.read_csv("GSE81089_FPKM.tsv.gz", sep="\t", index_col=0)

import numpy as np
import re

df_clean = df.drop(index=[i for i in df.index if not str(i).startswith('ENSG')])

df_clean = df_clean.replace(-1.0, np.nan)
df_clean = df_clean.fillna(0.0)

expr = df_clean.T  # samples x genes

log_expr = np.log1p(expr)
gene_var = log_expr.var(axis=0)
nonzero = gene_var > 1e-12
log_expr = log_expr.loc[:, nonzero]
gene_var = gene_var[nonzero]

labels = pd.Series(['Tumor' if c.split('_')[0].endswith('T') else ('Normal' if c.split('_')[0].endswith('N') else 'Unknown') for c in expr.index], index=expr.index)

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

N_GENES = 2000
sorted_genes = gene_var.sort_values(ascending=False)
hvg_genes = sorted_genes.index[:N_GENES]
bottom_genes = sorted_genes.index[-N_GENES:]

X_hvg_raw = log_expr[hvg_genes].values
X_bottom_raw = log_expr[bottom_genes].values

scaler = StandardScaler()
X_hvg_scaled = scaler.fit_transform(X_hvg_raw)

pca = PCA(n_components=50, random_state=42)
X_pca = pca.fit_transform(X_hvg_scaled)

from ripser import ripser

def compute_ph(X, maxdim=2, thresh=np.inf):
    res = ripser(X, maxdim=maxdim, thresh=thresh)
    return res['dgms']

conditions = {
    "real_hvg_raw": X_hvg_raw,
    "real_bottom_raw": X_bottom_raw,
    "real_pca50": X_pca,
}

diagrams = {}
for name, X in conditions.items():
    dgms = compute_ph(X, maxdim=2)
    diagrams[name] = dgms

def h1_stats_from_dgms(dgms):
    h1 = dgms[1]
    h1 = h1[np.isfinite(h1[:,1])]
    if len(h1) == 0:
        return {'count': 0, 'total_pers': 0.0, 'max_pers': 0.0}
    pers = h1[:,1] - h1[:,0]
    return {'count': int(len(h1)), 'total_pers': float(pers.sum()), 'max_pers': float(pers.max())}

def perm_pvalue(observed, null_dist):
    return float((1 + np.sum(null_dist >= observed)) / (1 + len(null_dist)))

obs_hvg_raw = h1_stats_from_dgms(diagrams['real_hvg_raw'])
obs_hvg_pca = h1_stats_from_dgms(diagrams['real_pca50'])
obs_bottom_raw = h1_stats_from_dgms(diagrams['real_bottom_raw'])

# Pipeline-symmetric permutation null
def h1_stats(dgms):
    h1 = dgms[1]
    h1 = h1[np.isfinite(h1[:,1])]
    if len(h1) == 0:
        return {'count': 0, 'total_pers': 0.0, 'max_pers': 0.0}
    pers = h1[:,1] - h1[:,0]
    return {'count': int(len(h1)), 'total_pers': float(pers.sum()), 'max_pers': float(pers.max())}

def run_permutation_null(X_real_fixed_genes, n_perm=2000, seed=42, n_pcs=50):
    rng = np.random.RandomState(seed)
    n, d = X_real_fixed_genes.shape
    null_raw = {'count': [], 'total_pers': [], 'max_pers': []}
    null_pca = {'count': [], 'total_pers': [], 'max_pers': []}
    for i in range(n_perm):
        perm = np.empty_like(X_real_fixed_genes)
        for j in range(d):
            perm[:, j] = rng.permutation(X_real_fixed_genes[:, j])
        dgm_raw = ripser(perm, maxdim=1)['dgms']
        s_raw = h1_stats(dgm_raw)
        for k in s_raw: null_raw[k].append(s_raw[k])

        scaled = StandardScaler().fit_transform(perm)
        pcaX = PCA(n_components=n_pcs, random_state=seed).fit_transform(scaled)
        dgm_pca = ripser(pcaX, maxdim=1)['dgms']
        s_pca = h1_stats(dgm_pca)
        for k in s_pca: null_pca[k].append(s_pca[k])
    return {k: np.array(v) for k, v in null_raw.items()}, {k: np.array(v) for k, v in null_pca.items()}

N_PERM = 2000

null_raw_hvg, null_pca_hvg = run_permutation_null(X_hvg_raw, n_perm=N_PERM, seed=42)
null_raw_bottom, null_pca_bottom = run_permutation_null(X_bottom_raw, n_perm=N_PERM, seed=43)

nulls = {
    'null_raw_hvg': null_raw_hvg, 'null_pca_hvg': null_pca_hvg,
    'null_raw_bottom': null_raw_bottom, 'null_pca_bottom': null_pca_bottom
}

# Gaussian null
def run_gaussian_null(n_samples, n_genes, n_draws=500, seed=100, n_pcs=50):
    rng = np.random.RandomState(seed)
    null_raw = {'count': [], 'total_pers': [], 'max_pers': []}
    null_pca = {'count': [], 'total_pers': [], 'max_pers': []}
    for i in range(n_draws):
        X = rng.normal(size=(n_samples, n_genes))
        dgm_raw = ripser(X, maxdim=1)['dgms']
        s_raw = h1_stats_from_dgms(dgm_raw)
        for k in s_raw: null_raw[k].append(s_raw[k])
        scaled = StandardScaler().fit_transform(X)
        pcaX = PCA(n_components=n_pcs, random_state=seed).fit_transform(scaled)
        dgm_pca = ripser(pcaX, maxdim=1)['dgms']
        s_pca = h1_stats_from_dgms(dgm_pca)
        for k in s_pca: null_pca[k].append(s_pca[k])
    return {k: np.array(v) for k, v in null_raw.items()}, {k: np.array(v) for k, v in null_pca.items()}

N_DRAWS = 500
gauss_null_raw, gauss_null_pca = run_gaussian_null(n_samples=218, n_genes=2000, n_draws=N_DRAWS, seed=100)

# Build final results table
rows = []
def add_row(label, feature_set, condition, observed, null_dist, null_name):
    p = perm_pvalue(observed, null_dist)
    z = (observed - null_dist.mean())/null_dist.std() if null_dist.std()>0 else np.nan
    rows.append({
        'dataset': label, 'features': feature_set, 'condition': condition,
        'observed_max_pers': observed, 'null_mean': null_dist.mean(), 'null_std': null_dist.std(),
        'p_value': p, 'z_score': z, 'null_model': null_name, 'n_null_draws': len(null_dist)
    })

add_row('GSE81089 NSCLC', 'HVG top-2000', 'raw', obs_hvg_raw['max_pers'], nulls['null_raw_hvg']['max_pers'], 'pipeline-symmetric permutation')
add_row('GSE81089 NSCLC', 'HVG top-2000', 'PCA50', obs_hvg_pca['max_pers'], nulls['null_pca_hvg']['max_pers'], 'pipeline-symmetric permutation')
add_row('GSE81089 NSCLC', 'bottom-2000 var', 'raw', obs_bottom_raw['max_pers'], nulls['null_raw_bottom']['max_pers'], 'pipeline-symmetric permutation')
add_row('GSE81089 NSCLC', 'HVG top-2000', 'raw', obs_hvg_raw['max_pers'], gauss_null_raw['max_pers'], 'Gaussian noise (OQ-006)')
add_row('GSE81089 NSCLC', 'HVG top-2000', 'PCA50', obs_hvg_pca['max_pers'], gauss_null_pca['max_pers'], 'Gaussian noise (OQ-006)')

final_rows = []

def add_final(comparison, statistic, real_or_diff_mean, null_mean, null_std, p_value, effect_size, n_null, validation_method):
    final_rows.append({
        'comparison': comparison, 'statistic': statistic,
        'observed_or_diff': real_or_diff_mean, 'null_mean': null_mean, 'null_std': null_std,
        'p_value': p_value, 'effect_size_cohens_d': effect_size, 'n_null_draws': n_null,
        'validation_method': validation_method
    })

add_final('Real HVG-raw vs pipeline-null', 'max_H1_persistence', obs_hvg_raw['max_pers'],
          nulls['null_raw_hvg']['max_pers'].mean(), nulls['null_raw_hvg']['max_pers'].std(),
          perm_pvalue(obs_hvg_raw['max_pers'], nulls['null_raw_hvg']['max_pers']),
          (obs_hvg_raw['max_pers']-nulls['null_raw_hvg']['max_pers'].mean())/nulls['null_raw_hvg']['max_pers'].std(),
          2000, 'Permutation test (per-gene shuffle across samples), 2000 permutations')

add_final('Real PCA50 vs pipeline-null-PCA50', 'max_H1_persistence', obs_hvg_pca['max_pers'],
          nulls['null_pca_hvg']['max_pers'].mean(), nulls['null_pca_hvg']['max_pers'].std(),
          perm_pvalue(obs_hvg_pca['max_pers'], nulls['null_pca_hvg']['max_pers']),
          (obs_hvg_pca['max_pers']-nulls['null_pca_hvg']['max_pers'].mean())/nulls['null_pca_hvg']['max_pers'].std(),
          2000, 'Permutation test (per-gene shuffle, then PCA(50) on permuted data), 2000 permutations')

add_final('Real bottom-var-raw vs pipeline-null', 'max_H1_persistence', obs_bottom_raw['max_pers'],
          nulls['null_raw_bottom']['max_pers'].mean(), nulls['null_raw_bottom']['max_pers'].std(),
          perm_pvalue(obs_bottom_raw['max_pers'], nulls['null_raw_bottom']['max_pers']),
          np.nan, 2000, 'Permutation test, bottom-2000-variance genes, 2000 permutations')

add_final('Real HVG-raw vs Gaussian-null', 'max_H1_persistence', obs_hvg_raw['max_pers'],
          gauss_null_raw['max_pers'].mean(), gauss_null_raw['max_pers'].std(),
          perm_pvalue(obs_hvg_raw['max_pers'], gauss_null_raw['max_pers']),
          (obs_hvg_raw['max_pers']-gauss_null_raw['max_pers'].mean())/gauss_null_raw['max_pers'].std(),
          500, 'Monte Carlo null (500 iid Gaussian draws, dimension-matched)')

add_final('Real PCA50 vs Gaussian-null-PCA50', 'max_H1_persistence', obs_hvg_pca['max_pers'],
          gauss_null_pca['max_pers'].mean(), gauss_null_pca['max_pers'].std(),
          perm_pvalue(obs_hvg_pca['max_pers'], gauss_null_pca['max_pers']),
          (obs_hvg_pca['max_pers']-gauss_null_pca['max_pers'].mean())/gauss_null_pca['max_pers'].std(),
          500, 'Monte Carlo null (500 iid Gaussian draws), PCA(50) on each draw')

diff_gauss = gauss_null_pca['max_pers'] - gauss_null_raw['max_pers']
add_final('Gaussian null: PCA50 vs raw (paired)', 'max_H1_persistence_delta', diff_gauss.mean(), 0.0, diff_gauss.std(),
           float(__import__('scipy.stats', fromlist=['wilcoxon']).wilcoxon(gauss_null_pca['max_pers'], gauss_null_raw['max_pers'])[1]),
           diff_gauss.mean()/diff_gauss.std(), 500, 'Paired Wilcoxon signed-rank test, 500 paired draws')

diff_pipeline = nulls['null_pca_hvg']['max_pers'] - nulls['null_raw_hvg']['max_pers']
add_final('Pipeline null: PCA50 vs raw (paired)', 'max_H1_persistence_delta', diff_pipeline.mean(), 0.0, diff_pipeline.std(),
           float(__import__('scipy.stats', fromlist=['wilcoxon']).wilcoxon(nulls['null_pca_hvg']['max_pers'], nulls['null_raw_hvg']['max_pers'])[1]),
           diff_pipeline.mean()/diff_pipeline.std(), 2000, 'Paired Wilcoxon signed-rank test, 2000 paired permutations')

final_df = pd.DataFrame(final_rows)
final_df.to_csv("final_results_table.csv", index=False)
pd.set_option('display.width', 200)
pd.set_option('display.max_colwidth', 60)
print(final_df.to_string(index=False))