"""
TOPOLOGICA ablation sweep -- Step 2: run all 18 configurations.

For each config: compute the real-data H1 persistence statistic (raw-gene
space and PCA-projected space) and compare it against both null models
(pipeline-symmetric per-gene permutation null, and dimension-matched
Gaussian-noise null) at that config's draw count.

Reuses the pre-registered GSE81089 pilot pipeline (log1p -> HVG selection ->
z-score -> PCA -> Vietoris-Rips persistent homology via ripser), varying
HVG count, PC count, imputation rule, and permutation draw count per
configs.json.

REQUIRES: data/imputation_variants.pkl (produced by 00_download_and_preprocess.py)
          configs.json (produced by 01_make_configs.py)

Writes (checkpointed after every config, so the run is resumable):
    results/sweep_results.json      -- one record per config (all statistics)
    results/nulldist_<tag>.pkl      -- per-config null distributions + PH diagrams

COMPUTE NOTE: a full run at the pre-registered n_perm=2000 (used only for the
convergence-check configs) plus n_perm=500 elsewhere costs several hours on
4 CPU workers. To smoke-test the pipeline quickly, pass --quick to cap every
config's n_perm/n_draws at a small number (this changes the null distributions
and therefore the z-scores/p-values, but NOT the real-data statistics
real_raw_max_pers / real_pca_max_pers / real_pca_delta, which are
deterministic given the data and do not depend on permutation count --
verify_default.py checks those to full precision regardless of --quick).

Run: python 02_run_sweep.py [--quick N] [--workers N]
"""
import argparse
import json
import os
import pickle
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from ripser import ripser
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

_ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_ROOT, "data")
RESULTS_DIR = os.path.join(_ROOT, "results")


def load_imputation_data():
    with open(os.path.join(DATA_DIR, "imputation_variants.pkl"), "rb") as f:
        return pickle.load(f)  # {'zero':(log_expr, gene_var), 'mean':..., 'median':..., 'drop':...}


def get_hvg_matrix(imputation_data, imputation_rule, n_hvg):
    log_expr, gene_var = imputation_data[imputation_rule]
    sorted_genes = gene_var.sort_values(ascending=False)
    hvg = sorted_genes.index[:n_hvg]
    return log_expr[hvg].values, list(hvg)


def h1_stats(dgms):
    h1 = dgms[1]
    h1 = h1[np.isfinite(h1[:, 1])]
    if len(h1) == 0:
        return {'count': 0, 'total_pers': 0.0, 'max_pers': 0.0}
    pers = h1[:, 1] - h1[:, 0]
    return {'count': int(len(h1)), 'total_pers': float(pers.sum()), 'max_pers': float(pers.max()),
            'top_bar_idx': int(np.argmax(pers))}


def real_stats(X_raw, n_pcs, seed=42):
    dgm_raw = ripser(X_raw, maxdim=1)['dgms']
    s_raw = h1_stats(dgm_raw)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_raw)
    npcs_eff = min(n_pcs, X_raw.shape[0] - 1, X_raw.shape[1])
    pca = PCA(n_components=npcs_eff, random_state=seed)
    X_pca = pca.fit_transform(X_scaled)
    dgm_pca = ripser(X_pca, maxdim=1)['dgms']
    s_pca = h1_stats(dgm_pca)
    return s_raw, s_pca, dgm_raw, dgm_pca


def _one_perm_iter(args):
    X, seed_i, n_pcs = args
    rng = np.random.RandomState(seed_i)
    n, d = X.shape
    perm = np.empty_like(X)
    for j in range(d):
        perm[:, j] = rng.permutation(X[:, j])
    dgm_raw = ripser(perm, maxdim=1)['dgms']
    s_raw = h1_stats(dgm_raw)
    scaled = StandardScaler().fit_transform(perm)
    npcs_eff = min(n_pcs, perm.shape[0] - 1, perm.shape[1])
    pcaX = PCA(n_components=npcs_eff, random_state=seed_i).fit_transform(scaled)
    dgm_pca = ripser(pcaX, maxdim=1)['dgms']
    s_pca = h1_stats(dgm_pca)
    return s_raw['max_pers'], s_pca['max_pers'], s_raw['count'], s_pca['count']


def _one_gauss_iter(args):
    n_samples, n_genes, seed_i, n_pcs = args
    rng = np.random.RandomState(seed_i)
    X = rng.normal(size=(n_samples, n_genes))
    dgm_raw = ripser(X, maxdim=1)['dgms']
    s_raw = h1_stats(dgm_raw)
    scaled = StandardScaler().fit_transform(X)
    npcs_eff = min(n_pcs, n_samples - 1, n_genes)
    pcaX = PCA(n_components=npcs_eff, random_state=seed_i).fit_transform(scaled)
    dgm_pca = ripser(pcaX, maxdim=1)['dgms']
    s_pca = h1_stats(dgm_pca)
    return s_raw['max_pers'], s_pca['max_pers'], s_raw['count'], s_pca['count']


def run_pipeline_null_parallel(X, n_perm, base_seed, n_pcs, pool, n_workers):
    args = [(X, base_seed + i, n_pcs) for i in range(n_perm)]
    results = pool.map(_one_perm_iter, args, chunksize=max(1, n_perm // (n_workers * 4)))
    arr = np.array(results)  # n_perm x 4 [raw_max, pca_max, raw_count, pca_count]
    return {'raw_max_pers': arr[:, 0], 'pca_max_pers': arr[:, 1],
            'raw_count': arr[:, 2], 'pca_count': arr[:, 3]}


def run_gaussian_null_parallel(n_samples, n_genes, n_draws, base_seed, n_pcs, pool, n_workers):
    args = [(n_samples, n_genes, base_seed + i, n_pcs) for i in range(n_draws)]
    results = pool.map(_one_gauss_iter, args, chunksize=max(1, n_draws // (n_workers * 4)))
    arr = np.array(results)
    return {'raw_max_pers': arr[:, 0], 'pca_max_pers': arr[:, 1],
            'raw_count': arr[:, 2], 'pca_count': arr[:, 3]}


def zscore(observed, null_dist):
    sd = null_dist.std()
    if sd == 0:
        return np.nan
    return (observed - null_dist.mean()) / sd


def perm_pvalue(observed, null_dist):
    return float((1 + np.sum(null_dist >= observed)) / (1 + len(null_dist)))


def run_config(tag, imputation_rule, n_hvg, n_pcs, n_perm_pipeline, n_draws_gaussian,
               imputation_data, pool, n_workers, gauss_seed_base=100, perm_seed_base=42,
               n_samples=218):
    t0 = time.time()
    X_raw, hvg_genes = get_hvg_matrix(imputation_data, imputation_rule, n_hvg)
    s_raw, s_pca, dgm_raw, dgm_pca = real_stats(X_raw, n_pcs, seed=42)

    pipeline_null = run_pipeline_null_parallel(X_raw, n_perm_pipeline, perm_seed_base, n_pcs, pool, n_workers)
    gauss_null = run_gaussian_null_parallel(n_samples, n_hvg, n_draws_gaussian, gauss_seed_base, n_pcs, pool, n_workers)

    out = {
        'tag': tag, 'imputation_rule': imputation_rule, 'n_hvg': n_hvg, 'n_pcs': n_pcs,
        'n_perm_pipeline': n_perm_pipeline, 'n_draws_gaussian': n_draws_gaussian,
        'real_raw_max_pers': s_raw['max_pers'], 'real_pca_max_pers': s_pca['max_pers'],
        'real_raw_count': s_raw['count'], 'real_pca_count': s_pca['count'],
        'real_pca_delta': s_pca['max_pers'] - s_raw['max_pers'],
        'pipeline_null_raw_mean': float(pipeline_null['raw_max_pers'].mean()),
        'pipeline_null_raw_std': float(pipeline_null['raw_max_pers'].std()),
        'pipeline_null_pca_mean': float(pipeline_null['pca_max_pers'].mean()),
        'pipeline_null_pca_std': float(pipeline_null['pca_max_pers'].std()),
        'gauss_null_raw_mean': float(gauss_null['raw_max_pers'].mean()),
        'gauss_null_raw_std': float(gauss_null['raw_max_pers'].std()),
        'gauss_null_pca_mean': float(gauss_null['pca_max_pers'].mean()),
        'gauss_null_pca_std': float(gauss_null['pca_max_pers'].std()),
        'z_raw_vs_pipeline': zscore(s_raw['max_pers'], pipeline_null['raw_max_pers']),
        'z_pca_vs_pipeline': zscore(s_pca['max_pers'], pipeline_null['pca_max_pers']),
        'z_raw_vs_gauss': zscore(s_raw['max_pers'], gauss_null['raw_max_pers']),
        'z_pca_vs_gauss': zscore(s_pca['max_pers'], gauss_null['pca_max_pers']),
        'p_raw_vs_pipeline': perm_pvalue(s_raw['max_pers'], pipeline_null['raw_max_pers']),
        'p_pca_vs_pipeline': perm_pvalue(s_pca['max_pers'], pipeline_null['pca_max_pers']),
        'p_raw_vs_gauss': perm_pvalue(s_raw['max_pers'], gauss_null['raw_max_pers']),
        'p_pca_vs_gauss': perm_pvalue(s_pca['max_pers'], gauss_null['pca_max_pers']),
        'pipeline_pca_delta_mean': float((pipeline_null['pca_max_pers'] - pipeline_null['raw_max_pers']).mean()),
        'pipeline_pca_delta_std': float((pipeline_null['pca_max_pers'] - pipeline_null['raw_max_pers']).std()),
        'gauss_pca_delta_mean': float((gauss_null['pca_max_pers'] - gauss_null['raw_max_pers']).mean()),
        'gauss_pca_delta_std': float((gauss_null['pca_max_pers'] - gauss_null['raw_max_pers']).std()),
        'hvg_genes': hvg_genes,
        'compute_seconds': time.time() - t0,
    }
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(os.path.join(RESULTS_DIR, f"nulldist_{tag}.pkl"), "wb") as f:
        pickle.dump({'pipeline_null': pipeline_null, 'gauss_null': gauss_null,
                     'dgm_raw': dgm_raw, 'dgm_pca': dgm_pca}, f)
    print(f"[{tag}] done in {out['compute_seconds']:.1f}s | real_raw={out['real_raw_max_pers']:.3f} "
          f"real_pca={out['real_pca_max_pers']:.3f} delta={out['real_pca_delta']:.3f} "
          f"z_raw_pipe={out['z_raw_vs_pipeline']:.2f} z_pca_pipe={out['z_pca_vs_pipeline']:.2f}", flush=True)
    sys.stdout.flush()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", type=int, default=None,
                     help="Cap n_perm_pipeline and n_draws_gaussian at this value for every "
                          "config (smoke-test mode; real-data statistics are unaffected).")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--tags", type=str, default=None,
                     help="Comma-separated subset of config tags to run (default: all 18).")
    args = ap.parse_args()

    imputation_data = load_imputation_data()
    os.makedirs(RESULTS_DIR, exist_ok=True)

    configs = json.load(open(os.path.join(_ROOT, "configs.json")))
    if args.tags:
        wanted = set(args.tags.split(","))
        configs = [c for c in configs if c['tag'] in wanted]
    if args.quick is not None:
        for c in configs:
            c['n_perm_pipeline'] = min(c['n_perm_pipeline'], args.quick)
            c['n_draws_gaussian'] = min(c['n_draws_gaussian'], args.quick)
        print(f"--quick {args.quick}: capping permutation/draw counts "
              f"(z-scores/p-values will differ from the full run; real-data "
              f"statistics will not).", flush=True)

    results_path = os.path.join(RESULTS_DIR, "sweep_results.json")
    all_results = []
    done_tags = set()
    if os.path.exists(results_path):
        prev = json.load(open(results_path))
        all_results = prev
        done_tags = {r['tag'] for r in prev}
        print(f"Resuming: {len(done_tags)} configs already done", flush=True)

    pool = Pool(args.workers)
    try:
        for cfg in configs:
            if cfg['tag'] in done_tags:
                continue
            try:
                res = run_config(pool=pool, n_workers=args.workers,
                                  imputation_data=imputation_data, **cfg)
                all_results.append(res)
            except Exception as e:
                import traceback
                print(f"[{cfg['tag']}] FAILED: {e}\n{traceback.format_exc()}", flush=True)
                all_results.append({'tag': cfg['tag'], 'FAILED': True, 'error': str(e), **cfg})
            with open(results_path, "w") as f:
                json.dump(all_results, f, indent=2, default=str)
    finally:
        pool.close()
        pool.join()
    print("SWEEP COMPLETE", flush=True)


if __name__ == "__main__":
    main()
