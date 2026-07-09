import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from ripser import ripser

META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")

    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
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


def focal_palette(labels, focal, focal_color, other="muted", base_colors=None):
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    focal_set = {focal} if isinstance(focal, str) else set(focal)
    n = len(labels)
    if not focal_set & set(labels):
        raise ValueError(f"focal {focal!r} not found in labels")
    if base_colors is None:
        base_colors = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["#444444"])
    base_colors = [base_colors[i % len(base_colors)] for i in range(n)]
    if other == "grey":
        rest = ["#BCBCBC"] * n
    elif other == "ordinal":
        nf = max(1, n - len(focal_set))
        ramp = [mcolors.to_hex((v, v, v)) for v in
                ([0.55] if nf == 1 else [0.80 - 0.35 * i / (nf - 1) for i in range(nf)])]
        rest, k = [], 0
        for l in labels:
            rest.append(ramp[min(k, nf - 1)]); k += (l not in focal_set)
    else:
        def mute(c):
            r, g, b = mcolors.to_rgb(c)
            m = (r + g + b) / 3
            return mcolors.to_hex((0.3 * r + 0.7 * m, 0.3 * g + 0.7 * m, 0.3 * b + 0.7 * m))
        rest = [mute(c) for c in base_colors]
    return [focal_color if l in focal_set else rest[i] for i, l in enumerate(labels)]


apply_figure_style()

import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_CHECKPOINTS = os.path.join(_ROOT, 'results', 'pilot_GSE81089', 'checkpoints')

with open(os.path.join(_CHECKPOINTS, 'diagrams.pkl'), 'rb') as f:
    diagrams = pickle.load(f)

def get_h1(dgms):
    h1 = dgms[1]
    return h1[np.isfinite(h1[:,1])]

h1_real_raw = get_h1(diagrams['real_hvg_raw'])
h1_real_pca = get_h1(diagrams['real_pca50'])

rng_viz = np.random.RandomState(7)
pooled_raw = []
pooled_pca = []
for i in range(15):
    X = rng_viz.normal(size=(218, 2000))
    d = ripser(X, maxdim=1)['dgms'][1]
    d = d[np.isfinite(d[:,1])]
    pooled_raw.append(d)
    Xs = StandardScaler().fit_transform(X)
    Xp = PCA(n_components=50, random_state=7).fit_transform(Xs)
    d2 = ripser(Xp, maxdim=1)['dgms'][1]
    d2 = d2[np.isfinite(d2[:,1])]
    pooled_pca.append(d2)
pooled_raw = np.vstack(pooled_raw)
pooled_pca = np.vstack(pooled_pca)

colors = focal_palette(['Real NSCLC', 'Gaussian null'], focal='Real NSCLC', focal_color='#C0392B')
real_color, null_color = colors

fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))

def plot_pd(ax, h1, color, label, alpha=0.5, s=14):
    if len(h1) == 0:
        return
    ax.scatter(h1[:,0], h1[:,1], s=s, color=color, alpha=alpha, label=label, edgecolors='none')

ax = axes[0]
plot_pd(ax, pooled_raw, null_color, 'Gaussian noise null (15 draws pooled)')
plot_pd(ax, h1_real_raw, real_color, 'Real NSCLC (raw HVG)', alpha=0.9, s=26)
lims = [0, max(h1_real_raw[:,1].max(), pooled_raw[:,1].max())*1.1]
ax.plot(lims, lims, '--', color='0.6', linewidth=0.8, zorder=0)
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel('Birth'); ax.set_ylabel('Death')
ax.set_title('Raw feature space (2000 HVGs)')
ax.legend(frameon=False, loc='lower right', fontsize=6)
ax.margins(0.04)

ax = axes[1]
plot_pd(ax, pooled_pca, null_color, 'Gaussian noise null, PCA50 (15 draws pooled)')
plot_pd(ax, h1_real_pca, real_color, 'Real NSCLC, PCA50', alpha=0.9, s=26)
lims = [0, max(h1_real_pca[:,1].max(), pooled_pca[:,1].max())*1.1]
ax.plot(lims, lims, '--', color='0.6', linewidth=0.8, zorder=0)
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel('Birth'); ax.set_ylabel('Death')
ax.set_title('PCA-reduced space (50 PCs)')
ax.legend(frameon=False, loc='lower right', fontsize=6)
ax.margins(0.04)

fig.suptitle('GSE81089 NSCLC bulk RNA-seq (n=218 samples): H1 persistence, real vs Gaussian-noise null (OQ-006)', fontsize=9)
fig.tight_layout()
fig.savefig('persistence_diagrams_real_vs_gaussian.png', dpi=300)