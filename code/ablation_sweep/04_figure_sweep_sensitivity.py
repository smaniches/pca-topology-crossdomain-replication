"""
TOPOLOGICA -- Figure 3 generation: ablation/hyperparameter-sensitivity sweep.

Recovered reproduction script for paper/figures/fig3_ablation_sensitivity.png,
following the same request as the earlier-recovered permutation-test code:
this figure previously existed only as a static PNG with no committed
generation script. Recovered via host.lineage on the session frame that
produced it (frame_id 1bbce5c0-c78d-49e7-8fa7-c0478e2f0aaf), then rewritten
here to read exclusively from this directory's already-committed aggregate
CSVs (../../results/ablation_sweep/*.csv) rather than from the raw
sweep_results.json / imputation_variants.pkl the original interactive
session used -- those two raw files are not retained in the repository (see
Data/Code Availability in the manuscript for which intermediates are kept),
but every number plotted here is present in the aggregate tables that
02_run_sweep.py + 03_aggregate_results.py already produce and commit.

Panel (a): three sub-panels from ablation_sweep_full_table.csv +
permutation_convergence_table.csv --
  (left)   gene-count sweep at PC=50, zero-fill, n_perm=500
  (center) PC-count sweep at HVG=2000, zero-fill, n_perm=500
  (right)  permutation-count convergence (z-score with bootstrap 95% CI)

Panel (b): two sub-panels from interaction_check_table.csv +
imputation_overlap_table.csv --
  (left)  gene-count x PC-count interaction check at the four grid corners
  (right) imputation-rule HVG gene-set Jaccard overlap heatmap

Run from the repository root:
    python code/ablation_sweep/04_figure_sweep_sensitivity.py
Produces: paper/figures/fig3_ablation_sensitivity.png
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd

MASTER_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(MASTER_DIR, "results", "ablation_sweep")
FIGURES_OUT = os.path.join(MASTER_DIR, "paper", "figures")
TMP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_tmp_panels")
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(FIGURES_OUT, exist_ok=True)

META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


# ---------------------------------------------------------------------------
# Panel (a): sweep_sensitivity (3 sub-panels)
# ---------------------------------------------------------------------------
df_sweep = pd.read_csv(os.path.join(RESULTS_DIR, "ablation_sweep_full_table.csv"))
conv_df = pd.read_csv(os.path.join(RESULTS_DIR, "permutation_convergence_table.csv"))

apply_figure_style()
fig_a, axes = plt.subplots(1, 3, figsize=(11, 3.6))

hvg_sub = df_sweep[(df_sweep.n_pcs == 50) & (df_sweep.imputation_rule == "zero") & (df_sweep.n_perm_pipeline == 500)].sort_values("n_hvg")
ax = axes[0]
ax.plot(hvg_sub.n_hvg, hvg_sub.real_pca_delta, "o-", color="#2166ac", label="Real PCA-delta")
ax.plot(hvg_sub.n_hvg, hvg_sub.gauss_pca_delta_mean, "s--", color="#b2182b", label="Gaussian-null PCA-delta")
ax.plot(hvg_sub.n_hvg, hvg_sub.pipeline_pca_delta_mean, "^--", color="#f4a582", label="Pipeline-null PCA-delta")
ax.axhline(0, color=META_GREY, lw=0.8)
ax.set_xlabel("HVG gene count")
ax.set_ylabel("PCA − raw max-H1-persistence")
ax.set_title("Gene-count sweep (PC=50, zero-fill)")
ax.legend(frameon=False, fontsize=6, loc="upper left")

pc_sub = df_sweep[(df_sweep.n_hvg == 2000) & (df_sweep.imputation_rule == "zero") & (df_sweep.n_perm_pipeline == 500)].sort_values("n_pcs")
ax = axes[1]
ax.plot(pc_sub.n_pcs, pc_sub.real_pca_delta, "o-", color="#2166ac", label="Real PCA-delta")
ax.plot(pc_sub.n_pcs, pc_sub.gauss_pca_delta_mean, "s--", color="#b2182b", label="Gaussian-null PCA-delta")
ax.plot(pc_sub.n_pcs, pc_sub.pipeline_pca_delta_mean, "^--", color="#f4a582", label="Pipeline-null PCA-delta")
ax.axhline(0, color=META_GREY, lw=0.8)
ax.set_xlabel("PCA component count")
ax.set_ylabel("PCA − raw max-H1-persistence")
ax.set_title("PC-count sweep (HVG=2000, zero-fill)")

ax = axes[2]
ax.plot(conv_df.n_perm, conv_df.z_raw_mean, "o-", color="#2166ac", label="z (raw HVG)")
ax.fill_between(conv_df.n_perm, conv_df.z_raw_ci_lo, conv_df.z_raw_ci_hi, color="#2166ac", alpha=0.2)
ax.plot(conv_df.n_perm, conv_df.z_pca_mean, "s-", color="#b2182b", label="z (PCA50)")
ax.fill_between(conv_df.n_perm, conv_df.z_pca_ci_lo, conv_df.z_pca_ci_hi, color="#b2182b", alpha=0.2)
ax.set_xscale("log")
ax.set_xlabel("Permutation count (n_perm)")
ax.set_ylabel("z-score vs pipeline-null (bootstrap 95% CI)")
ax.set_title("Permutation-count convergence")
ax.legend(frameon=False, fontsize=6, loc="center left")

fig_a.tight_layout()
panel_top_path = os.path.join(TMP_DIR, "sweep_sensitivity.png")
fig_a.savefig(panel_top_path, dpi=300)
plt.close(fig_a)

# ---------------------------------------------------------------------------
# Panel (b): interaction check + imputation Jaccard overlap (2 sub-panels)
# ---------------------------------------------------------------------------
interaction_df = pd.read_csv(os.path.join(RESULTS_DIR, "interaction_check_table.csv"))
overlap_pairs = pd.read_csv(os.path.join(RESULTS_DIR, "imputation_overlap_table.csv"))

# sigma of real_pca_delta across the full 18-config sweep (matches original)
sigma = df_sweep["real_pca_delta"].std()

# Reconstruct the symmetric 4x4 gene-set Jaccard matrix from the pairwise table
rules_order = ["zero", "mean", "median", "drop"]
overlap_matrix = np.ones((4, 4))
pair_lookup = {}
for _, row in overlap_pairs.iterrows():
    a, b = [s.strip() for s in row["pair"].split(" vs ")]
    pair_lookup[(a, b)] = row["jaccard_hvg_genes"]
    pair_lookup[(b, a)] = row["jaccard_hvg_genes"]
for i, r1 in enumerate(rules_order):
    for j, r2 in enumerate(rules_order):
        if i != j:
            overlap_matrix[i, j] = pair_lookup[(r1, r2)]

apply_figure_style()
fig_b2, axes2 = plt.subplots(1, 2, figsize=(9, 3.6))

ax = axes2[0]
corner_order = ["HVG500_PC10", "HVG500_PC100", "HVG4000_PC10", "HVG4000_PC100"]
interaction_df = interaction_df.set_index("corner").loc[corner_order].reset_index()
labels_int = ["HVG500\nPC10", "HVG500\nPC100", "HVG4000\nPC10", "HVG4000\nPC100"]
x = np.arange(len(labels_int))
width = 0.35
ax.bar(x - width / 2, interaction_df["additive_prediction"], width, label="Additive prediction (LOO sum)", color="#f4a582")
ax.bar(x + width / 2, interaction_df["observed_joint"], width, label="Observed joint effect", color="#2166ac")
ax.axhline(0.5 * sigma, color=META_GREY, lw=0.8, linestyle=":")
ax.axhline(-0.5 * sigma, color=META_GREY, lw=0.8, linestyle=":")
ax.set_xticks(x)
ax.set_xticklabels(labels_int, fontsize=6)
ax.set_ylabel("Change in real PCA-delta from center point")
ax.set_title("Gene-count × PC-count interaction check")
ax.legend(frameon=False, fontsize=6, loc="upper left")

ax = axes2[1]
im = ax.imshow(overlap_matrix, cmap="YlGnBu", vmin=0.9, vmax=1.0)
ax.set_xticks(range(4)); ax.set_xticklabels(rules_order, fontsize=7)
ax.set_yticks(range(4)); ax.set_yticklabels(rules_order, fontsize=7)
for i in range(4):
    for j in range(4):
        ax.text(j, i, f"{overlap_matrix[i, j]:.2f}", ha="center", va="center", fontsize=7,
                color="white" if overlap_matrix[i, j] > 0.97 else "black")
ax.set_title("Imputation rule: HVG gene-set Jaccard overlap")
cbar = fig_b2.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize=6)

fig_b2.tight_layout()
panel_bot_path = os.path.join(TMP_DIR, "interaction_and_imputation_overlap.png")
fig_b2.savefig(panel_bot_path, dpi=300)
plt.close(fig_b2)

# ---------------------------------------------------------------------------
# Compose final Figure 3
# ---------------------------------------------------------------------------
apply_figure_style()
fig = plt.figure(figsize=(7.4, 6.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.25, top=0.93, bottom=0.03, left=0.02, right=0.98)

ax_top = fig.add_subplot(gs[0])
ax_top.imshow(mpimg.imread(panel_top_path)); ax_top.axis("off"); panel_letter(ax_top, "a")

ax_bot = fig.add_subplot(gs[1])
ax_bot.imshow(mpimg.imread(panel_bot_path)); ax_bot.axis("off"); panel_letter(ax_bot, "b")

fig.suptitle("Figure 3. Ablation/hyperparameter-sensitivity sweep (18 configurations, GSE81089): H1 fails at PC=10 and at HVG4000/PC100",
             fontsize=9.5, y=0.985)
out_path = os.path.join(FIGURES_OUT, "fig3_ablation_sensitivity.png")
fig.savefig(out_path, dpi=300, transparent=False)
plt.close(fig)
print(f"saved {out_path}")
