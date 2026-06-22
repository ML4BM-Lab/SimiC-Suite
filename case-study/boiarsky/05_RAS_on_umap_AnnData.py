# To run:
# nohup python3 -u scripts/02_tf_activity_umap.py > logs/02_tf_activity_umap.log 2>&1 &

# =============================================================================
# 02_tf_activity_umap.py
#
# Overlays SimiC TF activity scores (wAUC) onto the processed UMAP.
# Reads the already-computed h5ad and the wAUC matrix, joins them, and
# generates a grid of UMAP plots — one panel per TF.
# =============================================================================

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from adjustText import adjust_text
from pathlib import Path

HERE = Path(__file__).parent.resolve()

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

H5AD_PATH = HERE.parent / "output" / "GSE193531_processed.h5ad"

AUC_CSV = (
    HERE.parent.parent
    / "Boiarsky_run" / "outputSimic" / "matrices" / "NBM-SMM-MM"
    / "NBM-SMM-MM_L1_1e-06_L2_0.0001_wAUC_matrices_filtered_BIC_collected.csv"
)

FIGURES_DIR = HERE.parent / "output" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# TFs to visualize.  Set to None to plot all 100.
# Biologically relevant plasma-cell TFs to start with:
# TF_NAMES = [
#     "IRF4", "XBP1", "PRDM1", "MAF", "IKZF1",
#     "MYC", "RUNX1", "BACH1", "JUNB", "STAT3",
# ]
TF_NAMES = [
    "IRF2", "MAF", "POU2F2", "BCL11A", "FOXO3",
    "NR4A1", "E2F4", "SMAD5", "BACH1", "FOXN3",
]

# Leiden resolution to show on the left anchor panel (must exist in h5ad)
LEIDEN_KEY = "leiden_res1.5"

# Number of columns in the TF grid
N_COLS = 3

# colormap for TF activity scores
CMAP = "plasma"

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
print("Loading h5ad ...")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

print("Loading wAUC matrix ...")
auc = pd.read_csv(AUC_CSV, index_col=0)
print(f"  {auc.shape[0]} cells x {auc.shape[1]} TFs")

# Validate index overlap
missing = adata.obs.index.difference(auc.index)
if len(missing) > 0:
    print(f"WARNING: {len(missing)} cells in h5ad are NOT in the AUC matrix — they will be NaN.")
else:
    print("  All cell indices match.")

# Resolve TF list
all_tfs = auc.columns.tolist()
if TF_NAMES is None:
    tfs = all_tfs
else:
    tfs = [t for t in TF_NAMES if t in all_tfs]
    missing_tfs = [t for t in TF_NAMES if t not in all_tfs]
    if missing_tfs:
        print(f"WARNING: TFs not found in AUC matrix and will be skipped: {missing_tfs}")
    print(f"  Plotting {len(tfs)} TFs: {tfs}")

# Join activity scores into adata.obs (reindex to match adata order)
auc_sub = auc.loc[:, tfs].reindex(adata.obs.index)
for tf in tfs:
    adata.obs[f"{tf}_AS"] = auc_sub[tf].values

# -----------------------------------------------------------------------------
# Plot: anchor cluster UMAP + one panel per TF
# -----------------------------------------------------------------------------
n_panels = 1 + len(tfs)          # leiden cluster panel + one per TF
n_cols   = N_COLS
n_rows   = int(np.ceil(n_panels / n_cols))

fig, axes = plt.subplots(n_rows, n_cols,
                          figsize=(n_cols * 4, n_rows * 4))
axes = axes.flatten()

sc.settings.set_figure_params(dpi=120)

# Panel 0: Leiden clusters — color legend on the right + repelled labels on data
if LEIDEN_KEY in adata.obs.columns:
    sc.pl.umap(adata, color=LEIDEN_KEY, legend_loc="right margin",
               legend_fontsize=6, title=f"Leiden ({LEIDEN_KEY})",
               ax=axes[0], show=False)

    # Compute per-cluster median UMAP coordinate as label anchor
    umap_df = pd.DataFrame(
        adata.obsm["X_umap"], index=adata.obs.index, columns=["u1", "u2"]
    )
    umap_df["cluster"] = adata.obs[LEIDEN_KEY].values
    centroids = umap_df.groupby("cluster")[["u1", "u2"]].median()

    texts = [
        axes[0].text(row.u1, row.u2, str(cl),
                     fontsize=6, fontweight="bold", color="black",
                     bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.5, lw=0))
        for cl, row in centroids.iterrows()
    ]
    adjust_text(texts, ax=axes[0], arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
else:
    sc.pl.umap(adata, color="disease_stage",
               title="disease stage", ax=axes[0], show=False)

# One panel per TF
for i, tf in enumerate(tfs, start=1):
    sc.pl.umap(adata, color=f"{tf}_AS", cmap=CMAP,
               title=f"{tf} activity", ax=axes[i], show=False)

# Hide any unused axes
for j in range(n_panels, len(axes)):
    axes[j].set_visible(False)

plt.suptitle("TF activity scores on UMAP (SimiC wAUC)", y=1.01, fontsize=14)
plt.tight_layout()

out_path = FIGURES_DIR / "tf_activity_umap.pdf"
fig.savefig(out_path, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved -> {out_path}")
print("Done.")


# =============================================================================
# Per-TF figure: disease stage UMAP + TF activity UMAP (one page per TF)
# =============================================================================
from matplotlib.backends.backend_pdf import PdfPages

DISEASE_COLORS = ["#3B7EA1", "#E66101", "#B2182B"]

# Ensure disease_stage is categorical and pin custom colors
if "disease_stage" in adata.obs.columns:
    adata.obs["disease_stage"] = adata.obs["disease_stage"].astype("category")
    n_cats = len(adata.obs["disease_stage"].cat.categories)
    adata.uns["disease_stage_colors"] = DISEASE_COLORS[:n_cats]

out_per_tf = FIGURES_DIR / "tf_activity_per_tf.pdf"
print(f"\nGenerating per-TF figures -> {out_per_tf}")

with PdfPages(out_per_tf) as pdf:
    for tf in tfs:
        col = f"{tf}_AS"

        fig, (ax_stage, ax_tf) = plt.subplots(1, 2, figsize=(12, 5))

        # Left panel: disease stage, legend on right margin, no on-data labels
        sc.pl.umap(
            adata,
            color="disease_stage",
            legend_loc="right margin",
            legend_fontsize=9,
            title="Disease Stage",
            ax=ax_stage,
            show=False,
        )

        # Right panel: TF activity score
        sc.pl.umap(
            adata,
            color=col,
            cmap=CMAP,
            title=f"{tf} Activity",
            ax=ax_tf,
            show=False,
        )

        fig.suptitle(tf, fontsize=14, fontweight="bold", y=1.02)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        print(f"  {tf} done")

print(f"Saved per-TF figures -> {out_per_tf}")


# =============================================================================
# Focus figure: 2 x 2 UMAP grid for IRF2/QC review
# =============================================================================
IRF2_COL = "IRF2_AS"
NEO_COL = "normal_or_neoplastic"
RIBO_COL = "pct_counts_ribo"
focus_out = FIGURES_DIR / "IRF2_umaps.pdf"

if IRF2_COL not in adata.obs.columns:
    if "IRF2" not in auc.columns:
        raise KeyError("IRF2 was not found in the wAUC matrix.")
    adata.obs[IRF2_COL] = auc["IRF2"].reindex(adata.obs.index).values

if RIBO_COL not in adata.obs.columns:
    if "umi_counts" not in adata.layers:
        raise KeyError(
            "'pct_counts_ribo' is missing and adata.layers['umi_counts'] is not available "
            "to calculate it."
        )
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["ribo"],
        layer="umi_counts",
        percent_top=None,
        log1p=False,
        inplace=True,
    )

required_cols = ["disease_stage", NEO_COL, IRF2_COL, RIBO_COL]
missing_cols = [col for col in required_cols if col not in adata.obs.columns]
if missing_cols:
    raise KeyError(f"Missing required columns for focus UMAP: {missing_cols}")

adata.obs[NEO_COL] = adata.obs[NEO_COL].astype("category")


def normal_neo_color(label):
    label = str(label).lower()
    if "normal" in label:
        return "#2ca25f"
    if "neoplastic" in label:
        return "#de2d26"
    return "#808080"


adata.uns[f"{NEO_COL}_colors"] = [
    normal_neo_color(cat) for cat in adata.obs[NEO_COL].cat.categories
]

print(f"\nGenerating focus UMAP figure -> {focus_out}")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

sc.pl.umap(
    adata,
    color="disease_stage",
    legend_loc="right margin",
    legend_fontsize=9,
    title="Disease Stage",
    ax=axes[0],
    show=False,
)

sc.pl.umap(
    adata,
    color=NEO_COL,
    legend_loc="right margin",
    legend_fontsize=9,
    title="normal_or_neoplastic (Boiarsky et al. classification)",
    ax=axes[1],
    show=False,
)

sc.pl.umap(
    adata,
    color=IRF2_COL,
    cmap=CMAP,
    title="IRF2 Activity Score",
    ax=axes[2],
    show=False,
)

sc.pl.umap(
    adata,
    color=RIBO_COL,
    color_map="RdYlBu_r",
    title="% Ribosomal",
    ax=axes[3],
    show=False,
)

plt.tight_layout()
fig.savefig(focus_out, bbox_inches="tight")
plt.close(fig)
print(f"Saved focus UMAP figure -> {focus_out}")
