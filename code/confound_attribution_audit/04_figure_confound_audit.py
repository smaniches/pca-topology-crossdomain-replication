"""
TOPOLOGICA -- Figure 4 generation: confound-attribution audit across all four
cohorts.

Recovered reproduction script for paper/figures/fig4_confound_audit.png,
following the same request as the earlier-recovered permutation-test code:
this figure previously existed only as a static PNG with no committed
generation script. Recovered via host.lineage on the session frame that
produced it (frame_id 8e6e2a10-1cdc-4b6d-b902-b11f3a5c7b1f), then rewritten
here to read every plotted value from this directory's already-committed
table1-4 CSVs (per-cohort) instead of the hardcoded literal lists the
original interactive session used. Every number below was cross-checked
against its source CSV before this script was written; see the loader
functions for the exact row/column pulled.

Panel (a): within-class decomposition -- z-score (vs. Gaussian null) for
mixed / tumor-only / normal-only subsets, each cohort.
    SOURCE: table1_within_class_decomposition[_<COHORT>].csv, z_vs_gaussian
    (pilot/TCGA-LUAD column naming) or z_gaussian (CPTAC/GSE146889 naming).

Panel (b): within-stratum control -- the WEAKEST confound-quartile's
z-score (vs. Gaussian null) in each cohort.
    SOURCE: table2_within_stratum_control[_<COHORT>].csv, minimum of
    z_score / z_gaussian / z_gauss across the four quartile rows.

Panel (c): residualization in PCA(50) space (the pre-registered primary
test space) -- z vs. the more conservative permutation null, post
class-mean residualization.
    SOURCE: table3_residualization_control[_<COHORT>].csv,
    "...residualized HVG->PCA50" row's z_permutation column (pilot,
    CPTAC-CCRCC, GSE146889 rows are already in PCA50 space by construction;
    TCGA-LUAD uses its explicit "...residualized HVG->PCA50" row since that
    cohort's table also carries raw-HVG-space rows).

Panel (d): block-bootstrap 95% CI on the tumor-only signal-beyond-null
(z-scale).
    SOURCE: table4_block_bootstrap_ci[_<COHORT>].csv. Row selection differs
    slightly by cohort because each cohort's audit reported this statistic
    under a different subset label -- see loader below; this mirrors the
    manuscript's own per-cohort table (Table CI column) rather than
    reinterpreting the underlying audits.

Run from the repository root:
    python code/confound_attribution_audit/04_figure_confound_audit.py
Produces: paper/figures/fig4_confound_audit.png
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MASTER_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(MASTER_DIR, "results", "confound_attribution_audit")
FIGURES_OUT = os.path.join(MASTER_DIR, "paper", "figures")
os.makedirs(FIGURES_OUT, exist_ok=True)

META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
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


apply_figure_style()

datasets_order = ["GSE81089\n(pilot)", "CPTAC-CCRCC", "GSE146889", "TCGA-LUAD\n(RNA-seq)"]

# ---------------------------------------------------------------------------
# Panel A: within-class decomposition (z vs. Gaussian null)
# ---------------------------------------------------------------------------
t1_pilot = pd.read_csv(os.path.join(RESULTS_DIR, "table1_within_class_decomposition.csv"))
t1_cptac = pd.read_csv(os.path.join(RESULTS_DIR, "table1_within_class_decomposition_CPTAC_CCRCC.csv"))
t1_g146 = pd.read_csv(os.path.join(RESULTS_DIR, "table1_within_class_decomposition_GSE146889.csv"))
t1_tcga = pd.read_csv(os.path.join(RESULTS_DIR, "table1_within_class_decomposition_TCGA_LUAD.csv"))

mixed_z = [
    t1_pilot.loc[t1_pilot.condition == "mixed", "z_vs_gaussian"].iloc[0],
    t1_cptac.loc[t1_cptac.condition == "Mixed (tumor+normal)", "z_gaussian"].iloc[0],
    t1_g146.loc[t1_g146.condition == "Mixed (tumor+normal)", "z_gaussian"].iloc[0],
    t1_tcga.loc[t1_tcga.condition == "mixed", "z_vs_gaussian"].iloc[0],
]
tumoronly_z = [
    t1_pilot.loc[t1_pilot.condition == "tumor", "z_vs_gaussian"].iloc[0],
    t1_cptac.loc[t1_cptac.condition == "Tumor-only", "z_gaussian"].iloc[0],
    t1_g146.loc[t1_g146.condition == "Tumor-only", "z_gaussian"].iloc[0],
    t1_tcga.loc[t1_tcga.condition == "tumor", "z_vs_gaussian"].iloc[0],
]
normalonly_z = [
    t1_pilot.loc[t1_pilot.condition == "normal", "z_vs_gaussian"].iloc[0],
    t1_cptac.loc[t1_cptac.condition == "Normal-only", "z_gaussian"].iloc[0],
    t1_g146.loc[t1_g146.condition == "Normal-only", "z_gaussian"].iloc[0],
    t1_tcga.loc[t1_tcga.condition == "normal", "z_vs_gaussian"].iloc[0],
]

# ---------------------------------------------------------------------------
# Panel B: within-stratum control -- weakest quartile's z (vs. Gaussian null)
# ---------------------------------------------------------------------------
t2_pilot = pd.read_csv(os.path.join(RESULTS_DIR, "table2_within_stratum_control.csv"))
t2_cptac = pd.read_csv(os.path.join(RESULTS_DIR, "table2_within_stratum_control_CPTAC_CCRCC.csv"))
t2_g146 = pd.read_csv(os.path.join(RESULTS_DIR, "table2_within_stratum_control_GSE146889.csv"))
t2_tcga = pd.read_csv(os.path.join(RESULTS_DIR, "table2_within_stratum_control_TCGA_LUAD.csv"))

min_quartile_z_gauss = [
    t2_pilot["z_score"].min(),
    t2_cptac["z_gaussian"].min(),
    t2_g146["z_gauss"].min(),
    t2_tcga["z_score"].min(),
]

# ---------------------------------------------------------------------------
# Panel C: residualization in PCA(50) space -- z vs. permutation null
# ---------------------------------------------------------------------------
t3_pilot = pd.read_csv(os.path.join(RESULTS_DIR, "table3_residualization_control.csv"))
t3_cptac = pd.read_csv(os.path.join(RESULTS_DIR, "table3_residualization_control_CPTAC_CCRCC.csv"))
t3_g146 = pd.read_csv(os.path.join(RESULTS_DIR, "table3_residualization_control_GSE146889.csv"))
t3_tcga = pd.read_csv(os.path.join(RESULTS_DIR, "table3_residualization_control_TCGA_LUAD.csv"))

resid_z_perm = [
    t3_pilot.loc[t3_pilot.space.str.contains("residualized"), "z_vs_permutation"].iloc[0],
    t3_cptac.loc[t3_cptac.space.str.contains("residualized"), "z_permutation"].iloc[0],
    t3_g146.loc[t3_g146.space.str.contains("residualized"), "z_permutation"].iloc[0],
    t3_tcga.loc[t3_tcga.space == "Class-mean-residualized HVG->PCA50", "z_vs_permutation"].iloc[0],
]

# ---------------------------------------------------------------------------
# Panel D: block-bootstrap 95% CI on tumor-only signal-beyond-null (z-scale)
# ---------------------------------------------------------------------------
t4_pilot = pd.read_csv(os.path.join(RESULTS_DIR, "table4_block_bootstrap_ci.csv"))
t4_cptac = pd.read_csv(os.path.join(RESULTS_DIR, "table4_block_bootstrap_ci_CPTAC_CCRCC.csv"))
t4_g146 = pd.read_csv(os.path.join(RESULTS_DIR, "table4_block_bootstrap_ci_GSE146889.csv"))
t4_tcga = pd.read_csv(os.path.join(RESULTS_DIR, "table4_block_bootstrap_ci_TCGA_LUAD.csv"))

pilot_row = t4_pilot.iloc[0]
ci_lower_pilot, ci_upper_pilot, point_pilot = pilot_row["ci_low_z"], pilot_row["ci_high_z"], pilot_row["point_estimate_z"]

cptac_row = t4_cptac.loc[t4_cptac.statistic.str.startswith("z (Delta")].iloc[0]
ci_lower_cptac, ci_upper_cptac, point_cptac = cptac_row["ci_lower"], cptac_row["ci_upper"], cptac_row["point_estimate"]

g146_row = t4_g146.loc[t4_g146.subset == "tumor-only"].iloc[0]
ci_lower_g146, ci_upper_g146, point_g146 = g146_row["z_ci95_lo"], g146_row["z_ci95_hi"], g146_row["point_z"]

tcga_row = t4_tcga.iloc[0]
ci_lower_tcga, ci_upper_tcga, point_tcga = tcga_row["ci_low_z"], tcga_row["ci_high_z"], tcga_row["point_estimate_z"]

ci_lower = [ci_lower_pilot, ci_lower_cptac, ci_lower_g146, ci_lower_tcga]
ci_upper = [ci_upper_pilot, ci_upper_cptac, ci_upper_g146, ci_upper_tcga]
point = [point_pilot, point_cptac, point_g146, point_tcga]

# ---------------------------------------------------------------------------
# Compose figure
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
x = np.arange(4)

ax = axes[0, 0]
w = 0.25
ax.bar(x - w, mixed_z, width=w, label="Mixed", color="#4C72B0")
ax.bar(x, tumoronly_z, width=w, label="Tumor-only", color="#DD8452")
ax.bar(x + w, normalonly_z, width=w, label="Normal-only", color="#55A868")
ax.axhline(3.0, color="grey", linestyle="--", linewidth=1)
ax.set_xticks(x)
ax.set_xticklabels(datasets_order, fontsize=7)
ax.set_ylabel("z-score vs. Gaussian null")
ax.set_title("Within-class decomposition: does a single class alone carry the signal?")
ax.legend(frameon=False, fontsize=7, loc="upper right")
ax.margins(y=0.08)

ax = axes[0, 1]
bar_colors = ["#4C72B0" if v >= 3.0 else "#C44E52" for v in min_quartile_z_gauss]
ax.bar(x, min_quartile_z_gauss, color=bar_colors)
ax.axhline(3.0, color="grey", linestyle="--", linewidth=1)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(datasets_order, fontsize=7)
ax.set_ylabel("Weakest quartile's z-score (vs. Gaussian null)")
ax.set_title("Within-stratum control: weakest confound-quartile's signal")
ax.margins(y=0.08)

ax = axes[1, 0]
bar_colors2 = ["#4C72B0" if v >= 3.0 else "#C44E52" for v in resid_z_perm]
ax.bar(x, resid_z_perm, color=bar_colors2)
ax.axhline(3.0, color="grey", linestyle="--", linewidth=1)
ax.set_xticks(x)
ax.set_xticklabels(datasets_order, fontsize=7)
ax.set_ylabel("z vs. permutation null, post-residualization")
ax.set_title("Residualization in PCA(50) space (pre-registered test space)")
ax.margins(y=0.1)
for xi, v in zip(x, resid_z_perm):
    ax.text(xi, v + (0.6 if v >= 0 else -1.0), f"{v:.2f}", ha="center", fontsize=7)

ax = axes[1, 1]
for xi, lo, hi, pt in zip(x, ci_lower, ci_upper, point):
    color = "#4C72B0" if lo > 0 else "#C44E52"
    ax.plot([xi, xi], [lo, hi], color=color, linewidth=2.5)
    ax.plot(xi, pt, 'o', color=color, markersize=5)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_xticks(x)
ax.set_xticklabels(datasets_order, fontsize=7)
ax.set_ylabel("Bootstrap 95% CI, z (signal beyond null)")
ax.set_title("Block-bootstrap CI: is the tumor-only excess stable under resampling?")
ax.margins(y=0.1)

fig.suptitle("Confound-attribution audit generalizes cleanly in 2 of 4 cohorts; weaker/failing in the other 2", fontsize=10, y=1.00)
fig.tight_layout()
out_path = os.path.join(FIGURES_OUT, "fig4_confound_audit.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"saved {out_path}")
