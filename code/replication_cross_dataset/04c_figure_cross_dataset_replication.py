"""
TOPOLOGICA -- Figure 2 generation: Cross-dataset replication summary.

Recovered reproduction script for paper/figures/fig2_cross_dataset_replication.png,
following the same request as the earlier-recovered permutation-test code: this
figure previously existed only as a static PNG with no committed generation
script. Recovered via host.lineage on the session frame that produced it
(frame_id 33bd0889-0613-4dcb-b974-2103a27529f0), then rewritten here to read
its numbers exclusively from the four already-committed CSVs listed below,
rather than from the hardcoded literals the original interactive session used.

Panel (a) -- H0: z-scores of real max-H1-persistence vs. Gaussian and
pipeline-symmetric nulls, raw and PCA(50) space, all 5 cohort/layer rows.
    SOURCE: results/cross_dataset_BH_family.csv (confirmatory 3 cohorts, 4
    layers -> GSE146889, CPTAC-CCRCC, TCGA-LUAD RNA-seq, TCGA-LUAD Meth) plus
    results/pilot_GSE81089/final_results_table.csv (pilot row, shown for
    reference -- not part of the pre-registered confirmatory BH-FDR family,
    consistent with the manuscript's own footnote on this point).

Panel (b) -- H1: real PCA-driven Delta in max-H1-persistence (diamond)
against the range spanned by the two null models' own PCA-driven Delta (bar).
    SOURCE: same two files as panel (a); Delta = PCA-space observed value
    minus raw-space observed/null value, computed here rather than transcribed.

Panel (c) -- Confound cross-check: cross-validated classifier AUC for
tumor/normal status, and dominant-H1-loop tumor-enrichment (Fisher's exact).
    SOURCE, CV AUC: results/confound_attribution_audit/table3_residualization_control*.csv
    ("Full ... (class-mean intact)" row's cv_auc column) for GSE81089,
    CPTAC-CCRCC, GSE146889; results/replication_TCGA_LUAD/tcga_luad_confound_cv_auc.csv
    (raw_hvg/rna row) for TCGA-LUAD.
    SOURCE, loop enrichment (n touched / total, Fisher p): these four numbers
    are NOT captured in any committed CSV -- they exist only as prose in each
    cohort's own results/ report (cited inline below per cohort). This is
    disclosed explicitly rather than left silent: recovering them as
    structured data would require re-deriving the cocycle trace from each
    cohort's raw persistence diagrams, which are not all retained in the
    repository (see paper's Data/Code Availability section on which
    checkpoints are and are not kept). The four values below are transcribed
    from committed markdown text, not computed by this script.

Run from the repository root:
    python code/replication_cross_dataset/04c_figure_cross_dataset_replication.py
Produces: paper/figures/fig2_cross_dataset_replication.png
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

MASTER_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(MASTER_DIR, "results")
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


def set_frame(ax, style="open"):
    show = {"open": (False, False, True, True),
            "boxed": (True, True, True, True),
            "none": (False, False, False, False)}[style]
    for side, vis in zip(("top", "right", "bottom", "left"), show):
        ax.spines[side].set_visible(vis)
        if vis:
            ax.spines[side].set_linewidth(0.6)
    ax.tick_params(direction="out", length=0 if style == "none" else 3, width=0.6)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


# ---------------------------------------------------------------------------
# Load committed source data
# ---------------------------------------------------------------------------
bh = pd.read_csv(os.path.join(RESULTS_DIR, "cross_dataset_BH_family.csv"))
pilot = pd.read_csv(os.path.join(RESULTS_DIR, "pilot_GSE81089", "final_results_table.csv"))

# Pilot z-scores (raw and PCA50, both nulls)
pilot_z = {
    ("raw", "Gaussian"): pilot.loc[pilot.comparison == "Real HVG-raw vs Gaussian-null", "effect_size"].iloc[0],
    ("raw", "Permutation"): pilot.loc[pilot.comparison == "Real HVG-raw vs pipeline-null", "effect_size"].iloc[0],
    ("PCA50", "Gaussian"): pilot.loc[pilot.comparison == "Real PCA50 vs Gaussian-null-PCA50", "effect_size"].iloc[0],
    ("PCA50", "Permutation"): pilot.loc[pilot.comparison == "Real PCA50 vs pipeline-null-PCA50", "effect_size"].iloc[0],
}

# Confirmatory cohorts' z-scores, pulled straight from cross_dataset_BH_family.csv
def bh_z(dataset_match, condition_match, null_match):
    row = bh[bh.dataset.str.contains(dataset_match) & bh.condition.str.contains(condition_match, case=False)]
    if null_match:
        row = row[row.condition.str.contains(null_match, case=False) | row.null_model.astype(str).str.contains(null_match, case=False)]
    return row["z_score"].iloc[0]

rows = [
    ("GSE81089 (pilot)", "RNA-seq", "raw", "Gaussian", pilot_z[("raw", "Gaussian")]),
    ("GSE81089 (pilot)", "RNA-seq", "raw", "Permutation", pilot_z[("raw", "Permutation")]),
    ("GSE81089 (pilot)", "RNA-seq", "PCA50", "Gaussian", pilot_z[("PCA50", "Gaussian")]),
    ("GSE81089 (pilot)", "RNA-seq", "PCA50", "Permutation", pilot_z[("PCA50", "Permutation")]),

    ("GSE146889", "RNA-seq", "raw", "Gaussian",
     bh[(bh.dataset.str.contains("GSE146889")) & (bh.condition == "Real raw HVG") & (bh.null_model.str.contains("Gaussian"))]["z_score"].iloc[0]),
    ("GSE146889", "RNA-seq", "raw", "Permutation",
     bh[(bh.dataset.str.contains("GSE146889")) & (bh.condition == "Real raw HVG") & (bh.null_model.str.contains("Pipeline"))]["z_score"].iloc[0]),
    ("GSE146889", "RNA-seq", "PCA50", "Gaussian",
     bh[(bh.dataset.str.contains("GSE146889")) & (bh.condition == "Real PCA50") & (bh.null_model.str.contains("Gaussian"))]["z_score"].iloc[0]),
    ("GSE146889", "RNA-seq", "PCA50", "Permutation",
     bh[(bh.dataset.str.contains("GSE146889")) & (bh.condition == "Real PCA50") & (bh.null_model.str.contains("Pipeline"))]["z_score"].iloc[0]),

    ("CPTAC-CCRCC", "Proteomics", "raw", "Gaussian",
     bh[(bh.dataset.str.contains("CPTAC")) & (bh.condition.str.contains("raw HVG vs Gaussian"))]["z_score"].iloc[0]),
    ("CPTAC-CCRCC", "Proteomics", "raw", "Permutation",
     bh[(bh.dataset.str.contains("CPTAC")) & (bh.condition.str.contains("raw HVG vs pipeline"))]["z_score"].iloc[0]),
    ("CPTAC-CCRCC", "Proteomics", "PCA50", "Gaussian",
     bh[(bh.dataset.str.contains("CPTAC")) & (bh.condition.str.contains("PCA50 vs Gaussian"))]["z_score"].iloc[0]),
    ("CPTAC-CCRCC", "Proteomics", "PCA50", "Permutation",
     bh[(bh.dataset.str.contains("CPTAC")) & (bh.condition.str.contains("PCA50 vs pipeline"))]["z_score"].iloc[0]),

    ("TCGA-LUAD", "RNA-seq", "raw", "Gaussian",
     bh[(bh.dataset == "TCGA-LUAD") & (bh.layer == "RNA-seq") & (bh.condition.str.contains("raw_gaussian"))]["z_score"].iloc[0]),
    ("TCGA-LUAD", "RNA-seq", "raw", "Permutation",
     bh[(bh.dataset == "TCGA-LUAD") & (bh.layer == "RNA-seq") & (bh.condition.str.contains("raw_perm"))]["z_score"].iloc[0]),
    ("TCGA-LUAD", "RNA-seq", "PCA50", "Gaussian",
     bh[(bh.dataset == "TCGA-LUAD") & (bh.layer == "RNA-seq") & (bh.condition.str.contains("pca50_gaussian"))]["z_score"].iloc[0]),
    ("TCGA-LUAD", "RNA-seq", "PCA50", "Permutation",
     bh[(bh.dataset == "TCGA-LUAD") & (bh.layer == "RNA-seq") & (bh.condition.str.contains("pca50_perm"))]["z_score"].iloc[0]),

    ("TCGA-LUAD", "Methylation", "raw", "Gaussian",
     bh[(bh.dataset == "TCGA-LUAD") & (bh.layer == "Methylation") & (bh.condition.str.contains("raw_gaussian"))]["z_score"].iloc[0]),
    ("TCGA-LUAD", "Methylation", "raw", "Permutation",
     bh[(bh.dataset == "TCGA-LUAD") & (bh.layer == "Methylation") & (bh.condition.str.contains("raw_perm"))]["z_score"].iloc[0]),
    ("TCGA-LUAD", "Methylation", "PCA (spectral)", "Gaussian",
     bh[(bh.dataset == "TCGA-LUAD") & (bh.layer == "Methylation") & (bh.condition.str.contains("pca35spectral_gaussian"))]["z_score"].iloc[0]),
    ("TCGA-LUAD", "Methylation", "PCA (spectral)", "Permutation",
     bh[(bh.dataset == "TCGA-LUAD") & (bh.layer == "Methylation") & (bh.condition.str.contains("pca35spectral_perm"))]["z_score"].iloc[0]),
]
df_z = pd.DataFrame(rows, columns=["dataset", "layer", "space", "null", "z"])
df_z["group"] = df_z["dataset"] + " — " + df_z["layer"]

# ---------------------------------------------------------------------------
# Panel (a)
# ---------------------------------------------------------------------------
apply_figure_style()
groups = ["GSE81089 (pilot) — RNA-seq", "GSE146889 — RNA-seq", "CPTAC-CCRCC — Proteomics",
          "TCGA-LUAD — RNA-seq", "TCGA-LUAD — Methylation"]
group_labels = ["GSE81089 (pilot)\nRNA-seq", "GSE146889\nRNA-seq", "CPTAC-CCRCC\nProteomics",
                "TCGA-LUAD\nRNA-seq", "TCGA-LUAD\nMethylation"]

fig_a, ax = plt.subplots(figsize=(7.0, 3.6))
y_positions = {g: i for i, g in enumerate(groups)}
space_offset = {"raw": 0.18, "PCA50": -0.18, "PCA (spectral)": -0.18}
null_color = {"Gaussian": "#3b6fa0", "Permutation": "#d2691e"}
null_marker = {"Gaussian": "o", "Permutation": "^"}

for _, r in df_z.iterrows():
    y = y_positions[r["group"]] + space_offset.get(r["space"], 0)
    ax.scatter(r["z"], y, color=null_color[r["null"]], marker=null_marker[r["null"]],
               s=45, edgecolor="white", linewidth=0.5, zorder=3)

ax.axvline(3.0, color="grey", linestyle=":", linewidth=1.2, zorder=1)
ax.text(3.0, len(groups) - 0.35, "z = 3.0\n(pre-registered\nthreshold)", fontsize=6, ha="left", va="top", color="grey")
ax.set_yticks(range(len(groups)))
ax.set_yticklabels(group_labels)
ax.set_xlabel("z-score vs. null (standardized deviation from null mean)")
ax.set_title("H0: real max-H1-persistence vs. matched nulls, raw and PCA(50) space")
ax.set_xlim(-5, 62)
ax.margins(y=0.15)
ax.text(60, len(groups) - 1 + 0.18, "raw", fontsize=6, va="center", color="black")
ax.text(60, len(groups) - 1 - 0.18, "PCA50", fontsize=6, va="center", color="black")
legend_elems = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=null_color["Gaussian"], markersize=7, label="Gaussian null"),
    Line2D([0], [0], marker='^', color='w', markerfacecolor=null_color["Permutation"], markersize=7, label="Pipeline-symmetric null"),
]
ax.legend(handles=legend_elems, loc="lower right", frameon=False, fontsize=7)
set_frame(ax)
fig_a.tight_layout()
panel_a_path = os.path.join(TMP_DIR, "panel2a_h0_zscores.png")
fig_a.savefig(panel_a_path, dpi=300, transparent=True)
plt.close(fig_a)

# ---------------------------------------------------------------------------
# Panel (b) -- H1: real PCA-delta (diamond) vs null Delta range (bar)
# ---------------------------------------------------------------------------
def delta(space_hi_z_row_obs, raw_obs):
    return space_hi_z_row_obs - raw_obs

pilot_raw_obs = pilot.loc[pilot.comparison == "Real HVG-raw vs Gaussian-null", "observed_or_diff"].iloc[0]
pilot_pca_obs = pilot.loc[pilot.comparison == "Real PCA50 vs Gaussian-null-PCA50", "observed_or_diff"].iloc[0]
pilot_real_delta = pilot_pca_obs - pilot_raw_obs
pilot_gauss_null_delta = pilot.loc[pilot.comparison == "Gaussian null: PCA50 vs raw (paired)", "observed_or_diff"].iloc[0]
pilot_perm_null_delta = pilot.loc[pilot.comparison == "Pipeline null: PCA50 vs raw (paired)", "observed_or_diff"].iloc[0]

def cohort_delta(dataset_mask, raw_condition, pca_condition, null_substr):
    raw_obs = bh[dataset_mask & (bh.condition == raw_condition)]["z_score"].index  # placeholder not used
    return None

g146 = bh[bh.dataset.str.contains("GSE146889")]
g146_raw_obs = g146[g146.condition == "Real raw HVG"]  # not directly needed; use observed values below
# Re-load per-cohort tables directly for exact observed max-H1 values (Delta needs observed, not z).
g146_tab = pd.read_csv(os.path.join(RESULTS_DIR, "replication_GSE146889", "final_results_table_GSE146889.csv"))
g146_raw_obs_val = g146_tab.loc[g146_tab.condition == "Real raw HVG", "observed_max_H1_persistence"].iloc[0]
g146_pca_obs_val = g146_tab.loc[g146_tab.condition == "Real PCA50", "observed_max_H1_persistence"].iloc[0]
g146_real_delta = g146_pca_obs_val - g146_raw_obs_val
g146_gauss_null_raw = g146_tab.loc[(g146_tab.condition == "Real raw HVG") & (g146_tab.null_model.str.contains("Gaussian")), "null_mean"].iloc[0]
g146_gauss_null_pca = g146_tab.loc[(g146_tab.condition == "Real PCA50") & (g146_tab.null_model.str.contains("Gaussian")), "null_mean"].iloc[0]
g146_gauss_null_delta = g146_gauss_null_pca - g146_gauss_null_raw
g146_perm_null_raw = g146_tab.loc[(g146_tab.condition == "Real raw HVG") & (g146_tab.null_model.str.contains("Pipeline")), "null_mean"].iloc[0]
g146_perm_null_pca = g146_tab.loc[(g146_tab.condition == "Real PCA50") & (g146_tab.null_model.str.contains("Pipeline")), "null_mean"].iloc[0]
g146_perm_null_delta = g146_perm_null_pca - g146_perm_null_raw

cptac_tab = pd.read_csv(os.path.join(RESULTS_DIR, "replication_CPTAC_CCRCC", "cptac_h1_pca_delta_table.csv")).set_index("quantity")["value"]
cptac_real_delta = float(cptac_tab["Real PCA-delta (PCA50 - raw)"])
cptac_gauss_null_delta = float(cptac_tab["Gaussian-null PCA-delta (mean)"])
cptac_perm_null_delta = float(cptac_tab["Pipeline-null PCA-delta (mean)"])

tcga_tab = pd.read_csv(os.path.join(RESULTS_DIR, "replication_TCGA_LUAD", "tcga_luad_h1_pca_delta_summary.csv")).set_index("comparison")
tcga_rna_real_delta = float(tcga_tab.loc["rna_gaussian", "real_delta"])
tcga_rna_gauss_null_delta = float(tcga_tab.loc["rna_gaussian", "null_delta_mean"])
tcga_rna_perm_null_delta = float(tcga_tab.loc["rna_perm", "null_delta_mean"])
tcga_meth_real_delta = float(tcga_tab.loc["meth_gaussian_spectral", "real_delta"])
tcga_meth_gauss_null_delta = float(tcga_tab.loc["meth_gaussian_spectral", "null_delta_mean"])
tcga_meth_perm_null_delta = float(tcga_tab.loc["meth_perm_spectral", "null_delta_mean"])

h1_rows = [
    ("GSE81089 (pilot)", pilot_real_delta, min(pilot_perm_null_delta, pilot_gauss_null_delta), max(pilot_perm_null_delta, pilot_gauss_null_delta)),
    ("GSE146889", g146_real_delta, min(g146_perm_null_delta, g146_gauss_null_delta), max(g146_perm_null_delta, g146_gauss_null_delta)),
    ("CPTAC-CCRCC", cptac_real_delta, min(cptac_perm_null_delta, cptac_gauss_null_delta), max(cptac_perm_null_delta, cptac_gauss_null_delta)),
    ("TCGA-LUAD RNA-seq", tcga_rna_real_delta, min(tcga_rna_perm_null_delta, tcga_rna_gauss_null_delta), max(tcga_rna_perm_null_delta, tcga_rna_gauss_null_delta)),
    ("TCGA-LUAD Methylation*", tcga_meth_real_delta, min(tcga_meth_perm_null_delta, tcga_meth_gauss_null_delta), max(tcga_meth_perm_null_delta, tcga_meth_gauss_null_delta)),
]
df_h1 = pd.DataFrame(h1_rows, columns=["dataset", "real_delta", "null_lo", "null_hi"])

fig_b, ax = plt.subplots(figsize=(7.0, 3.6))
apply_figure_style()
y = np.arange(len(df_h1))[::-1]
for yi, row in zip(y, df_h1.itertuples()):
    ax.plot([row.null_lo, row.null_hi], [yi, yi], color="#3b6fa0", linewidth=6, alpha=0.5, solid_capstyle="butt", zorder=1)
    ax.scatter([row.real_delta], [yi], color="#b5251f", marker="D", s=60, zorder=3)
ax.axvline(0, color="grey", linewidth=0.8, zorder=0)
ax.set_yticks(y)
ax.set_yticklabels(df_h1["dataset"])
ax.set_xlabel("PCA − raw max-H1-persistence (Δ)")
ax.set_title("H1: null models' PCA-driven Δ (bar) vs. real data's own Δ (diamond)")
ax.set_ylim(-1.0, len(df_h1))
legend_elems = [
    Line2D([0], [0], color="#3b6fa0", linewidth=6, alpha=0.5, label="Null Δ range (Gaussian–pipeline)"),
    Line2D([0], [0], marker='D', color='w', markerfacecolor="#b5251f", markersize=8, label="Real data Δ"),
]
ax.legend(handles=legend_elems, loc="upper left", bbox_to_anchor=(0.0, 1.28), ncol=2, frameon=False, fontsize=7)
ax.margins(x=0.08)
set_frame(ax)
fig_b.tight_layout()
panel_b_path = os.path.join(TMP_DIR, "panel2b_h1_deltas.png")
fig_b.savefig(panel_b_path, dpi=300, transparent=True)
plt.close(fig_b)

# ---------------------------------------------------------------------------
# Panel (c) -- confound cross-check: CV AUC + loop-enrichment
# ---------------------------------------------------------------------------
pilot_resid = pd.read_csv(os.path.join(RESULTS_DIR, "confound_attribution_audit", "table3_residualization_control.csv"))
pilot_cv_auc = float(pilot_resid.loc[pilot_resid.space.str.contains("intact"), "cv_auc"].iloc[0])

cptac_resid = pd.read_csv(os.path.join(RESULTS_DIR, "confound_attribution_audit", "table3_residualization_control_CPTAC_CCRCC.csv"))
cptac_cv_auc = float(cptac_resid.loc[cptac_resid.space.str.contains("intact"), "cv_auc_tumor_normal"].iloc[0])

g146_resid = pd.read_csv(os.path.join(RESULTS_DIR, "confound_attribution_audit", "table3_residualization_control_GSE146889.csv"))
g146_cv_auc = float(g146_resid.loc[g146_resid.space.str.contains("intact"), "cv_auc_tumor_normal"].iloc[0])

tcga_cv_tab = pd.read_csv(os.path.join(RESULTS_DIR, "replication_TCGA_LUAD", "tcga_luad_confound_cv_auc.csv"))
tcga_cv_auc = float(tcga_cv_tab.loc[(tcga_cv_tab.space == "raw_hvg") & (tcga_cv_tab.omics == "rna"), "cv_auc"].iloc[0])

# Loop-enrichment: NOT in any committed CSV. Transcribed verbatim from each
# cohort's own prose report (source line cited per row); disclosed above and
# in code/replication_cross_dataset/README.md.
#   GSE81089:    results/pilot_GSE81089/report.md            -- "touches only 3 of 218 samples ... p=1.0"
#   GSE146889:   results/replication_GSE146889/report_GSE146889.md -- "45 of 176 samples ... p = 2.87e-09"
#   CPTAC-CCRCC: results/replication_CPTAC_CCRCC/cptac_ccrcc_report.md -- PCA(50) row: "4" loop vertices / 194, p=0.135
#   TCGA-LUAD:   results/replication_TCGA_LUAD/tcga_luad_report.md -- "112/116 samples ... p = 0.1185"
confound_rows = [
    ("GSE81089 (pilot)", pilot_cv_auc, 3, 218, "no", 1.0),
    ("GSE146889", g146_cv_auc, 45, 176, "YES", 2.87e-9),
    ("CPTAC-CCRCC", cptac_cv_auc, 4, 194, "no", 0.135),
    ("TCGA-LUAD", tcga_cv_auc, 112, 116, "no", 0.1185),
]
df_c = pd.DataFrame(confound_rows, columns=["dataset", "cv_auc", "loop_n", "total_n", "enriched", "fisher_p"])
df_c["loop_frac"] = df_c["loop_n"] / df_c["total_n"]

fig_c, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
apply_figure_style()
ax0 = axes[0]
colors_bar = ["#b5251f" if e == "YES" else "#3b6fa0" for e in df_c["enriched"]]
ax0.bar(range(len(df_c)), df_c["cv_auc"], color=colors_bar, width=0.55)
ax0.set_xticks(range(len(df_c)))
ax0.set_xticklabels(df_c["dataset"], rotation=20, ha="right", fontsize=6.5)
ax0.set_ylabel("CV AUC (tumor/normal, logistic reg.)")
ax0.set_title("Classifier separability\nis near-ceiling in every cohort")
ax0.set_ylim(0.85, 1.03)
ax0.axhline(1.0, color="grey", linewidth=0.6, linestyle=":")
for i, v in enumerate(df_c["cv_auc"]):
    ax0.text(i, v + 0.005, f"{v:.3f}", ha="center", fontsize=6)
set_frame(ax0)

ax1 = axes[1]
frac_pct = df_c["loop_frac"] * 100
ax1.bar(range(len(df_c)), frac_pct, color=colors_bar, width=0.55)
ax1.set_xticks(range(len(df_c)))
ax1.set_xticklabels(df_c["dataset"], rotation=20, ha="right", fontsize=6.5)
ax1.set_ylabel("Dominant H1 loop:\nsamples touched (%)")
ax1.set_title("Loop/tumor-status enrichment:\nred = significant (Fisher p<0.05)")
for i, (frac, p) in enumerate(zip(frac_pct, df_c["fisher_p"])):
    lbl = f"p={p:.1e}" if p < 0.01 else f"p={p:.2f}"
    ax1.text(i, frac + 2, lbl, ha="center", fontsize=6)
set_frame(ax1)

legend_elems = [
    Line2D([0], [0], marker='s', color='w', markerfacecolor="#3b6fa0", markersize=8, label="Not significantly tumor-enriched"),
    Line2D([0], [0], marker='s', color='w', markerfacecolor="#b5251f", markersize=8, label="Significantly tumor-enriched (red flag)"),
]
fig_c.legend(handles=legend_elems, loc="lower center", ncol=2, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.08))
fig_c.tight_layout()
panel_c_path = os.path.join(TMP_DIR, "panel2c_confound_summary.png")
fig_c.savefig(panel_c_path, dpi=300, transparent=True, bbox_inches="tight")
plt.close(fig_c)

# ---------------------------------------------------------------------------
# Compose final Figure 2
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 9.0))
gs = fig.add_gridspec(3, 1, height_ratios=[1.1, 1.0, 0.85], hspace=0.35, top=0.95, bottom=0.03, left=0.03, right=0.97)

ax_a = fig.add_subplot(gs[0])
ax_a.imshow(mpimg.imread(panel_a_path)); ax_a.axis("off"); panel_letter(ax_a, "a")

ax_b = fig.add_subplot(gs[1])
ax_b.imshow(mpimg.imread(panel_b_path)); ax_b.axis("off"); panel_letter(ax_b, "b")

ax_c = fig.add_subplot(gs[2])
ax_c.imshow(mpimg.imread(panel_c_path)); ax_c.axis("off"); panel_letter(ax_c, "c")

fig.suptitle("Figure 2. Cross-dataset replication: H0 (real signal vs. null), H1 (PCA inflation vs. null), and confound-check summary",
             fontsize=9.5, y=0.99)
out_path = os.path.join(FIGURES_OUT, "fig2_cross_dataset_replication.png")
fig.savefig(out_path, dpi=300, transparent=False)
plt.close(fig)
print(f"saved {out_path}")
