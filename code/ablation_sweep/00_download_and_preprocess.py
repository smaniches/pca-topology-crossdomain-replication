"""
TOPOLOGICA ablation sweep -- Step 0: fetch GSE81089 raw data and build the
4 imputation-rule variants of the preprocessed expression matrix used by
the sweep (run_sweep.py) and by the confound-attribution audit.

Reproduces from scratch: no local checkpoints required (though this project's
repo also has results/pilot_GSE81089/checkpoints/*.pkl from the original
pilot run -- this script does NOT read those; it always regenerates from the
public GEO accession so the sweep is independently reproducible).

Output: data/imputation_variants.pkl
    {'zero': (log_expr_df, gene_var_series), 'mean': (...), 'median': (...), 'drop': (...)}
    each log_expr_df is samples x genes (log1p-transformed, nonzero-variance-filtered).
Output: data/sample_labels.pkl  (Series: Tumor/Normal/Unknown per sample)

Run: python 00_download_and_preprocess.py
"""
import os
import urllib.request
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

_ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(_ROOT, "data")
os.makedirs(DATA_DIR, exist_ok=True)

GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81089/suppl/GSE81089_FPKM_cufflinks.tsv.gz"
RAW_PATH = os.path.join(DATA_DIR, "GSE81089_FPKM_cufflinks.tsv.gz")


def fetch_raw_data():
    if os.path.exists(RAW_PATH):
        print(f"Using cached {RAW_PATH}")
        return
    print(f"Downloading {GEO_URL} ...")
    urllib.request.urlretrieve(GEO_URL, RAW_PATH)
    print(f"Saved to {RAW_PATH}")


def load_expression():
    """Load raw FPKM table, keep only ENSG-indexed rows (genes)."""
    df = pd.read_csv(RAW_PATH, sep="\t", index_col=0)
    df_clean = df.loc[[i for i in df.index if str(i).startswith("ENSG")]]
    return df_clean


def preprocess_impute(expr_raw, rule="zero"):
    """Apply an imputation rule to the pipeline's -1.0 sentinel ("missing")
    values, matching the pre-registered pilot's convention that -1.0 in the
    cufflinks FPKM output denotes a missing/undefined value (not a true zero
    or negative expression)."""
    X = expr_raw.copy()
    mask = (X == -1.0)
    if rule == "zero":
        X = X.mask(mask, 0.0)
    elif rule == "mean":
        X = X.mask(mask, np.nan)
        X = X.apply(lambda col: col.fillna(col.mean()), axis=0)
    elif rule == "median":
        X = X.mask(mask, np.nan)
        X = X.apply(lambda col: col.fillna(col.median()), axis=0)
    elif rule == "drop":
        genes_with_sentinel = mask.any(axis=0)
        X = X.loc[:, ~genes_with_sentinel]
    else:
        raise ValueError(f"unknown imputation rule: {rule}")
    return X


def build_log_expr(X_imputed):
    log_expr = np.log1p(X_imputed)
    gene_var = log_expr.var(axis=0)
    nonzero = gene_var > 1e-12
    return log_expr.loc[:, nonzero], gene_var[nonzero]


def main():
    fetch_raw_data()
    df_clean = load_expression()
    print("ENSG rows:", df_clean.shape)

    n_sentinel = (df_clean.values == -1.0).sum()
    print(f"sentinel -1.0 count: {n_sentinel} of {df_clean.size} "
          f"({n_sentinel / df_clean.size:.4%})")
    genes_with_sentinel = (df_clean == -1.0).any(axis=1).sum()
    print("genes with >=1 sentinel:", genes_with_sentinel)

    expr_raw = df_clean.T  # samples x genes, sentinels still present
    print("expr_raw shape (samples x genes):", expr_raw.shape)

    labels = pd.Series(
        ['Tumor' if c.split('_')[0].endswith('T')
         else ('Normal' if c.split('_')[0].endswith('N') else 'Unknown')
         for c in expr_raw.index],
        index=expr_raw.index,
    )
    print(labels.value_counts().to_dict())

    data_to_save = {}
    for rule in ["zero", "mean", "median", "drop"]:
        X_imputed = preprocess_impute(expr_raw, rule)
        log_e, var_e = build_log_expr(X_imputed)
        data_to_save[rule] = (log_e, var_e)
        print(f"{rule}: {log_e.shape}")

    out_path = os.path.join(DATA_DIR, "imputation_variants.pkl")
    import pickle
    with open(out_path, "wb") as f:
        pickle.dump(data_to_save, f)
    labels.to_pickle(os.path.join(DATA_DIR, "sample_labels.pkl"))
    print(f"Saved {out_path} ({os.path.getsize(out_path)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
