"""
Shared helper functions for the TCGA-LUAD confirmatory replication cohort
(RNA-seq + methylation layers). Imported by 01_rnaseq_results_table.py and
02_methylation_results_table.py -- not intended to be run directly.

No hardcoded absolute paths: callers pass in an explicit data directory.
"""
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import laplacian
from scipy.sparse import csr_matrix
from ripser import ripser


def select_hvg(mat_genes_x_samples, n_top=2000, bottom=False):
    """mat: genes x samples dataframe (rows=genes/probes). Returns top or bottom
    n_top rows by variance across samples."""
    var = mat_genes_x_samples.var(axis=1)
    order = var.sort_values(ascending=False)
    if bottom:
        sel = order.index[-n_top:]
    else:
        sel = order.index[:n_top]
    return mat_genes_x_samples.loc[sel]


def standardize_samples_x_features(mat_genes_x_samples):
    """Input genes x samples -> output samples x genes, standardized
    zero-mean unit-variance per feature (gene/probe)."""
    X = mat_genes_x_samples.T.values  # samples x genes
    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    return (X - mu) / sd


def estimate_intrinsic_dim_mle(X, k=10):
    k = min(k, X.shape[0] - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    dists, _ = nn.kneighbors(X)
    dists = dists[:, 1:]
    log_ratios = np.log(dists[:, -1:] / np.maximum(dists[:, :-1], 1e-15))
    d_hat = (k - 1) / np.sum(log_ratios, axis=1)
    return float(np.median(d_hat))


def estimate_intrinsic_dim_twonn(X):
    nn = NearestNeighbors(n_neighbors=3).fit(X)
    dists, _ = nn.kneighbors(X)
    r1 = dists[:, 1]
    r2 = dists[:, 2]
    mu = r2 / np.maximum(r1, 1e-15)
    n = len(X)
    d_hat = n / np.sum(np.log(np.maximum(mu, 1.0 + 1e-10)))
    return float(d_hat)


def dim_gate(X, k=10):
    d_mle = estimate_intrinsic_dim_mle(X, k=k)
    d_twonn = estimate_intrinsic_dim_twonn(X)
    d_amb = X.shape[1]
    ratio = d_mle / d_amb
    regime = "JL" if ratio < 0.3 else ("transitional" if ratio < 0.7 else "high_intrinsic")
    return {"d_int_mle": d_mle, "d_int_twonn": d_twonn, "d_amb": d_amb, "ratio": ratio, "regime": regime}


def max_h1_persistence_from_D(D, maxdim=2):
    """D: precomputed distance matrix. Returns dict of max_H{dim}/n_H{dim}."""
    res = ripser(D, distance_matrix=True, maxdim=maxdim)
    dgms = res["dgms"]
    out = {}
    for dim in range(maxdim + 1):
        dgm = dgms[dim]
        finite = dgm[np.isfinite(dgm[:, 1])]
        if len(finite) == 0:
            out[f"max_H{dim}"] = 0.0
            out[f"n_H{dim}"] = 0
        else:
            pers = finite[:, 1] - finite[:, 0]
            out[f"max_H{dim}"] = float(pers.max())
            out[f"n_H{dim}"] = int(len(pers))
    out["dgms"] = dgms
    return out


def max_h1_persistence(X, maxdim=2):
    """X: samples x features point cloud (Euclidean metric)."""
    D = squareform(pdist(X, metric="euclidean"))
    return max_h1_persistence_from_D(D, maxdim=maxdim)


def compute_spectral_distance(X, k=15, n_components=None):
    """Section-5-mandated metric for the 'transitional' intrinsic-dimension
    regime: graph-Laplacian eigenmap distance in place of raw Euclidean
    distance on the PCA-reduced point cloud. Used for TCGA-LUAD methylation
    PCA(35) space, which falls in the transitional regime (see
    tcga_luad_report.md Section 5)."""
    nn = NearestNeighbors(n_neighbors=min(k + 1, X.shape[0])).fit(X)
    dists, indices = nn.kneighbors(X)
    n = len(X)
    kk = dists.shape[1] - 1
    rows, cols, vals = [], [], []
    sigma = np.median(dists[:, 1]) if dists.shape[1] > 1 else 1.0
    if sigma == 0:
        sigma = 1e-6
    for i in range(n):
        for j_idx in range(1, kk + 1):
            j = indices[i, j_idx]
            w = np.exp(-dists[i, j_idx] ** 2 / (2 * sigma ** 2))
            rows.extend([i, j])
            cols.extend([j, i])
            vals.extend([w, w])
    W = csr_matrix((vals, (rows, cols)), shape=(n, n))
    L = laplacian(W, normed=False)
    if n_components is None:
        n_components = min(n - 2, 30)
    L_dense = L.toarray().astype(float)
    eigenvalues, eigenvectors = np.linalg.eigh(L_dense)
    mask = eigenvalues > 1e-8
    eigenvalues = eigenvalues[mask][:n_components]
    eigenvectors = eigenvectors[:, mask][:, :n_components]
    scaled = eigenvectors / np.sqrt(eigenvalues)[np.newaxis, :]
    D = squareform(pdist(scaled, metric="euclidean"))
    return D


def compute_z_and_p(real_val, null_vals):
    """Standardized deviation from null (z-score) and empirical one-sided p-value."""
    null_vals = np.array(null_vals)
    null_mean = null_vals.mean()
    null_std = null_vals.std(ddof=1)
    z = (real_val - null_mean) / null_std if null_std > 0 else np.nan
    p = (np.sum(null_vals >= real_val) + 1) / (len(null_vals) + 1)
    return {"real": float(real_val), "null_mean": float(null_mean), "null_std": float(null_std),
            "z": float(z) if not np.isnan(z) else None, "p_empirical": float(p), "n_null": len(null_vals)}


def benjamini_hochberg(p_values, alpha=0.05):
    p = np.asarray(p_values, dtype=np.float64)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(m)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        val = ranked[i] * m / (i + 1)
        prev = min(prev, val)
        adjusted[order[i]] = prev
    adjusted = np.minimum(adjusted, 1.0)
    return (adjusted <= alpha), adjusted


def extract_h1(records):
    return np.array([r["max_H1"] for r in records])
