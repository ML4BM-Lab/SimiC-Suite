# To run w nohup:
# nohup python3 -u scripts/01_data_prep.py > logs/01_data_prep.log 2>&1 &
# 2895197
# nohup python3 -u scripts/01_data_prep.py > logs/01_data_prep_clustering.log 2>&1 &

# =============================================================================
# 01_data_prep.py
#
# Boiarsky et al. 2022 (NBM, SMM, MM): data preparation.
#
# Replicates the preprocessing pipeline from the authors' published notebook:
#   0_reproduce_results_from_raw_data.ipynb
#
# Pipeline order:
#   1. Download GSE193531 from GEO
#   2. Generate IG-locus gene list from GRCh38 GTF (mirrors Ig_genes.ipynb)
#   3. Build AnnData with all samples
#   4. Drop MGUS cells
#   5. Export raw counts for SIMIC (NBM / SMM / MM only)
#   6. Normalize + process exactly as the paper (log1p, HVG, hemoglobin
#      regression, scale, PCA, UMAP)
#
# Inputs (downloaded automatically if absent):
#   GEO:     GSE193531_umi-count-matrix.csv.gz
#            GSE193531_cell-level-metadata.csv.gz
#   Ensembl: Homo_sapiens.GRCh38.<release>.gtf.gz  (for IG gene list)
#
# Outputs (under boiarsky/output/):
#   SIMIC inputs:   singlecell_matrix.mtx, cell_ids.txt,
#                   genes_ids.txt, metadata.csv
#   Processed:      GSE193531_processed.h5ad
#   QC plots:       figures/
#
# Usage:
#   cd SimiC/boiarsky
#   python 01_data_prep.py
# =============================================================================

import re
import urllib.request
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import scipy.sparse
import scipy.io
import matplotlib.pyplot as plt
from pathlib import Path


try:
    import datatable as dt
    HAS_DT = True
except ImportError:
    HAS_DT = False
    print("datatable not found; falling back to pandas (slower).")

try:
    from gtfparse import read_gtf
    HAS_GTFPARSE = True
except ImportError:
    HAS_GTFPARSE = False
    print("gtfparse not found. Install with: pip install gtfparse")
    print("Falling back to symbol-prefix matching for IG genes.")

# -----------------------------------------------------------------------------
# 0. Configuration
# -----------------------------------------------------------------------------
ACCESSION  = "GSE193531"
DATA_DIR   = Path("data") / ACCESSION
OUTPUT_DIR = Path("output")
FIGURES_DIR = OUTPUT_DIR / "figures"

# Set to True to skip steps 1-6 and load the already-processed h5ad directly.
# Useful when you only need to (re)run clustering and UMAP plots.
REUSE_PROCESSED = False

# Set to True to skip everything except step 8 (composition plots).
# Requires a saved h5ad that already contains leiden_res* columns.
PLOTTING_ONLY = True

for d in [DATA_DIR, OUTPUT_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# GEO supplementary files
GEO_FTP = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/{ACCESSION}/suppl"
GEO_FILES = {
    "matrix":   f"{ACCESSION}_umi-count-matrix.csv.gz",
    "metadata": f"{ACCESSION}_cell-level-metadata.csv.gz",
}

# GTF for IG gene list (Ensembl GRCh38 release 109 — matches hg38 coordinates
# used in the authors' Ig_genes.ipynb; any GRCh38 release works because the
# IG locus coordinates are hardcoded)
ENSEMBL_RELEASE = "109"
GTF_FNAME  = f"Homo_sapiens.GRCh38.{ENSEMBL_RELEASE}.gtf.gz"
GTF_URL    = (
    f"https://ftp.ensembl.org/pub/release-{ENSEMBL_RELEASE}/gtf/homo_sapiens/"
    f"{GTF_FNAME}"
)
GTF_PATH   = DATA_DIR / GTF_FNAME
IG_GENES_FILE = DATA_DIR / "ig_locus_genes.txt"


out_h5ad = OUTPUT_DIR / f"{ACCESSION}_processed.h5ad"

if PLOTTING_ONLY or REUSE_PROCESSED:
    # -------------------------------------------------------------------------
    # Fast path: load already-processed h5ad and jump straight to clustering
    # (REUSE_PROCESSED) or plotting (PLOTTING_ONLY).
    # The file must exist (run once with both flags False first).
    # -------------------------------------------------------------------------
    label = "PLOTTING_ONLY" if PLOTTING_ONLY else "REUSE_PROCESSED"
    print(f"{label}=True: loading saved h5ad, skipping steps 1-6.")
    cd138_adata = sc.read_h5ad(out_h5ad)
    print(f"  Loaded: {cd138_adata.n_obs} cells x {cd138_adata.n_vars} genes")

else:
    # =========================================================================
    # STEP 1 — Download GEO supplementary files
    # =========================================================================
    print("=" * 60)
    print("STEP 1: Downloading GEO supplementary files")
    print("=" * 60)
    for fname in GEO_FILES.values():
        dest = DATA_DIR / fname
        if dest.exists():
            print(f"  {fname}: already present.")
        else:
            url = f"{GEO_FTP}/{fname}"
            print(f"  Downloading {fname} ...")
            urllib.request.urlretrieve(url, dest)
            print(f"  -> {dest}")

    mat_gz  = DATA_DIR / GEO_FILES["matrix"]
    meta_gz = DATA_DIR / GEO_FILES["metadata"]

    # =========================================================================
    # STEP 2 — Generate IG-locus gene list (mirrors authors' Ig_genes.ipynb)
    #
    # IGH locus: chr14  105,586,437 – 106,879,844
    # IGL locus: chr22   22,026,076 –  22,922,913
    # IGK locus: chr2    88,857,361 –  90,235,368
    #
    # Any gene whose start OR end overlaps the locus is included.
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP 2: Building IG-locus gene list")
    print("=" * 60)

    IG_LOCI = {
        "IGH": ("14",  105_586_437, 106_879_844),
        "IGL": ("22",   22_026_076,  22_922_913),
        "IGK": ("2",    88_857_361,  90_235_368),
    }

    if IG_GENES_FILE.exists():
        ig_genes = pd.read_csv(IG_GENES_FILE, header=None).iloc[:, 0].tolist()
        print(f"  ig_locus_genes.txt already present ({len(ig_genes)} genes).")

    elif HAS_GTFPARSE:
        if not GTF_PATH.exists():
            print(f"  Downloading GTF ({GTF_FNAME}) — this may take a few minutes...")
            urllib.request.urlretrieve(GTF_URL, GTF_PATH)
            print(f"  -> {GTF_PATH}")

        print("  Parsing GTF...")
        gene_loci = read_gtf(str(GTF_PATH))
        gene_loci = gene_loci[gene_loci["feature"] == "gene"]

        ig_genes = []
        for locus_name, (chrom, locus_start, locus_end) in IG_LOCI.items():
            hits = gene_loci.loc[
                (gene_loci["seqname"] == chrom) &
                (
                    ((locus_start <= gene_loci["start"]) & (gene_loci["start"] <= locus_end)) |
                    ((locus_start <= gene_loci["end"])   & (gene_loci["end"]   <= locus_end))
                ),
                "gene_name",
            ].drop_duplicates()
            print(f"  {locus_name}: {len(hits)} genes")
            ig_genes += hits.tolist()

        ig_genes = list(dict.fromkeys(ig_genes))   # deduplicate, preserve order
        print(f"  Total IG genes: {len(ig_genes)}")

        pd.Series(ig_genes).to_csv(IG_GENES_FILE, index=False, header=False)
        print(f"  Saved -> {IG_GENES_FILE}")

    else:
        print("  WARNING: gtfparse unavailable; IG genes will be matched by symbol")
        print("  prefix (^IG[HKL]|JCHAIN). Install gtfparse for the exact list.")
        ig_genes = None   # resolved later against the gene universe

    # =========================================================================
    # STEP 3 — Load data and build AnnData
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP 3: Loading UMI matrix and metadata")
    print("=" * 60)

    print("  Loading count matrix...")
    if HAS_DT:
        umi = dt.fread(str(mat_gz)).to_pandas()
    else:
        umi = pd.read_csv(mat_gz)

    first_col = umi.columns[0]
    umi.index = umi[first_col].tolist()
    umi.drop(columns=first_col, inplace=True)
    umi = umi.astype("float64").T          # now cells x genes
    print(f"  Matrix: {umi.shape[0]} cells x {umi.shape[1]} genes")

    print("  Loading cell-level metadata...")
    metadata = pd.read_csv(meta_gz, index_col=0)
    print(f"  Metadata: {metadata.shape[0]} cells x {metadata.shape[1]} columns")

    cd138_adata = anndata.AnnData(X=umi, obs=metadata)

    cd138_adata.obs["sample_ID"] = pd.Categorical(
        cd138_adata.obs["sample_ID"],
        categories=[
            "NBM-1","NBM-2","NBM-3","NBM-4","NBM-6","NBM-7","NBM-8","NBM-10","NBM-11",
            "MGUS-1","MGUS-2","MGUS-3","MGUS-4","MGUS-5","MGUS-6",
            "SMM-1","SMM-2","SMM-3","SMM-4","SMM-5","SMM-6","SMM-7",
            "SMM-8","SMM-9","SMM-10","SMM-11","SMM-12",
            "MM-1","MM-2","MM-3","MM-4","MM-5","MM-6","MM-7","MM-8",
        ],
        ordered=True,
    )
    cd138_adata.obs["disease_stage"] = pd.Categorical(
        cd138_adata.obs["disease_stage"],
        categories=["NBM", "MGUS", "SMM", "MM"],
        ordered=True,
    )

    cd138_adata.obs.loc[
        cd138_adata.obs["normal_or_neoplastic"] == "nan", "normal_or_neoplastic"
    ] = None
    if hasattr(cd138_adata.obs["normal_or_neoplastic"], "cat"):
        cd138_adata.obs["normal_or_neoplastic"] = (
            cd138_adata.obs["normal_or_neoplastic"].cat.remove_unused_categories()
        )

    cd138_adata.layers["umi_counts"] = scipy.sparse.csr_matrix(cd138_adata.X)

    print(f"  AnnData (all stages): {cd138_adata.n_obs} cells x {cd138_adata.n_vars} genes")
    print(f"  Stage breakdown:\n{cd138_adata.obs['disease_stage'].value_counts().sort_index()}")

    # =========================================================================
    # STEP 4 — Drop MGUS cells
    #
    # MGUS samples have ~73% normal plasma cells on average (mean tumor purity
    # ~27%), so "MGUS" as a SimiC label would conflate donor-normal PCs with
    # MGUS tumor cells. We exclude MGUS to keep phenotype labels clean for SimiC.
    # Downstream normalization and processing then runs on NBM / SMM / MM only.
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP 4: Removing MGUS cells")
    print("=" * 60)

    n_before = cd138_adata.n_obs
    cd138_adata = cd138_adata[cd138_adata.obs["disease_stage"] != "MGUS"].copy()
    n_after  = cd138_adata.n_obs
    print(f"  Removed {n_before - n_after} MGUS cells.")
    print(f"  Remaining: {n_after} cells (NBM / SMM / MM)")
    print(f"  Stage breakdown:\n{cd138_adata.obs['disease_stage'].value_counts().sort_index()}")

    # =========================================================================
    # STEP 5 — Export raw counts for SIMIC
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP 5: Exporting data for SIMIC")
    print("=" * 60)

    pd.DataFrame({"Cells": cd138_adata.obs_names}).to_csv(
        OUTPUT_DIR / "cell_ids.txt", index=False, header=False
    )
    pd.DataFrame({"Genes": cd138_adata.var_names}).to_csv(
        OUTPUT_DIR / "genes_ids.txt", index=False, header=False
    )
    cd138_adata.obs.to_csv(OUTPUT_DIR / "metadata.csv")

    m_raw = cd138_adata.layers["umi_counts"]
    scipy.io.mmwrite(str(OUTPUT_DIR / "singlecell_matrix.mtx"), m_raw)

    print("  Exported:")
    print(f"    cell_ids.txt        ({cd138_adata.n_obs} cells)")
    print(f"    genes_ids.txt       ({cd138_adata.n_vars} genes)")
    print(f"    metadata.csv")
    print(f"    singlecell_matrix.mtx")

    # =========================================================================
    # STEP 6 — Normalize and process (mirrors paper / authors' notebook exactly)
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP 6: Normalizing and processing (paper pipeline)")
    print("=" * 60)

    # 6a. Normalize: total count = 10,000, exclude genes >20% in any cell
    sc.pp.normalize_total(
        cd138_adata, target_sum=1e4,
        exclude_highly_expressed=True, max_fraction=0.2
    )
    sc.pp.log1p(cd138_adata)

    cd138_adata.layers["lognorm"] = cd138_adata.X.copy()
    cd138_adata.raw = cd138_adata
    print("  Normalization done.")

    # 6b. Highly variable genes
    sc.pp.highly_variable_genes(
        cd138_adata, layer="lognorm",
        min_mean=0.0125, max_mean=4, min_disp=0.5
    )

    if ig_genes is None:
        ig_genes = [g for g in cd138_adata.var_names
                    if re.match(r"^(IG[HKL]|JCHAIN)", g)]
        print(f"  IG genes matched by prefix: {len(ig_genes)}")

    sex_genes = ["XIST", "RPS4Y1"]
    exclude = set(ig_genes + sex_genes)
    cd138_adata.var.loc[cd138_adata.var_names.isin(exclude), "highly_variable"] = False
    n_hvg = cd138_adata.var["highly_variable"].sum()
    print(f"  HVGs after IG/sex exclusion: {n_hvg}")

    # 6c. Regress out hemoglobin contamination score
    hb_pattern = re.compile(r"^HB.*")
    hb_genes = [g for g in cd138_adata.var_names if hb_pattern.match(g)]
    for non_hb in ["HBEGF", "HBS1L", "HBP1"]:
        if non_hb in hb_genes:
            hb_genes.remove(non_hb)

    lognorm_hb = cd138_adata[:, hb_genes].layers["lognorm"]
    # Paper specifies mean (not sum) of log-normalized hemoglobin gene expression
    cd138_adata.obs["hemoglobin_score"] = np.array(lognorm_hb.mean(axis=1)).ravel()

    print("  Regressing out hemoglobin score (may take a few minutes)...")
    sc.pp.regress_out(cd138_adata, keys="hemoglobin_score")

    # 6d. Scale, PCA, neighbors, UMAP
    sc.pp.scale(cd138_adata, max_value=10)
    sc.tl.pca(cd138_adata, svd_solver="arpack")
    sc.pp.neighbors(cd138_adata, n_neighbors=15, n_pcs=14)
    sc.tl.umap(cd138_adata)

    cd138_adata.uns["disease_stage_colors"] = [
        "cornflowerblue", "#DC7209", "#880B0B"   # NBM, SMM, MM (MGUS removed)
    ]
    print("  Scale / PCA / neighbors / UMAP done.")


if not PLOTTING_ONLY:
    # =========================================================================
    # STEP 6e — Leiden clustering
    # Paper: resolution=1.5; we also run 0.5 and 1.0 for comparison.
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP 6e: Leiden clustering")
    print("=" * 60)

    for res in [0.5, 1.0, 1.5]:
        key = f"leiden_res{res}"
        sc.tl.leiden(cd138_adata, resolution=res, key_added=key)
        n_clusters = cd138_adata.obs[key].nunique()
        print(f"  res={res}: {n_clusters} clusters")

    # =========================================================================
    # STEP 6f — Basic UMAP plots
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP 6f: UMAP plots")
    print("=" * 60)

    sc.settings.figdir = str(FIGURES_DIR)
    sc.settings.set_figure_params(figsize=[6, 6], dpi=120, format="pdf")
    sc.pl.umap(cd138_adata, color="disease_stage",
               title="disease stage", save="_by_stage.pdf", show=False)
    sc.pl.umap(cd138_adata, color="sample_ID",
               title="sample", save="_by_sample.pdf", show=False)
    for res in [0.5, 1.0, 1.5]:
        key = f"leiden_res{res}"
        sc.pl.umap(cd138_adata, color=key,
                   title=f"Leiden (res={res})", legend_loc="on data",
                   save=f"_leiden_res{res}.pdf", show=False)
    print(f"  UMAP plots saved to {FIGURES_DIR}/")

    # =========================================================================
    # STEP 7 — Save processed AnnData
    # =========================================================================
    print("\n" + "=" * 60)
    print("STEP 7: Saving processed AnnData")
    print("=" * 60)

    cd138_adata.write_h5ad(out_h5ad)
    print(f"  Saved -> {out_h5ad}")


# =============================================================================
# STEP 8 — QC + Composition figures
# =============================================================================
print("\n" + "=" * 60)
print("STEP 8: QC + Composition figures")
print("=" * 60)

sc.settings.figdir = str(FIGURES_DIR)
sc.settings.set_figure_params(dpi=120, format="pdf")

# --- 8a. Compute QC metrics from raw UMI counts (if not already present) ---
if "pct_counts_mt" not in cd138_adata.obs.columns:
    cd138_adata.var["mt"]   = cd138_adata.var_names.str.startswith("MT-")
    cd138_adata.var["ribo"] = cd138_adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        cd138_adata, qc_vars=["mt", "ribo"],
        layer="umi_counts", percent_top=None, log1p=False, inplace=True,
    )
    print(f"  MT   range: {cd138_adata.obs['pct_counts_mt'].min():.1f}% – "
          f"{cd138_adata.obs['pct_counts_mt'].max():.1f}%")
    print(f"  Ribo range: {cd138_adata.obs['pct_counts_ribo'].min():.1f}% – "
          f"{cd138_adata.obs['pct_counts_ribo'].max():.1f}%")
else:
    print("  QC metrics already present, skipping recalculation.")

# --- 8b. UMAP colored by QC metrics ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.umap(cd138_adata, color="pct_counts_mt",
           color_map="RdYlGn_r", title="% mitochondrial",
           ax=axes[0], show=False)
sc.pl.umap(cd138_adata, color="pct_counts_ribo",
           color_map="RdYlBu_r", title="% ribosomal",
           ax=axes[1], show=False)
plt.tight_layout()
qc_path = FIGURES_DIR / "umap_qc_metrics.pdf"
fig.savefig(qc_path, bbox_inches="tight")
plt.close(fig)
print(f"  QC UMAP saved -> {qc_path.name}")

# Fixed colors for disease stages (NBM / SMM / MM)
stage_colors = {
    stage: col
    for stage, col in zip(
        ["NBM", "SMM", "MM"],
        cd138_adata.uns.get("disease_stage_colors",
                            ["cornflowerblue", "#DC7209", "#880B0B"]),
    )
}

for res in [0.5, 1.0, 1.5]:
    key = f"leiden_res{res}"
    n_clusters = cd138_adata.obs[key].nunique()

    fig, axes = plt.subplots(1, 3, figsize=(7 + n_clusters * 0.35, 6))

    # --- Panel 1: UMAP colored by Leiden clusters ---
    sc.pl.umap(cd138_adata, color=key, legend_loc="on data",
               title=f"Leiden (res={res})", ax=axes[0], show=False)

    # --- Panel 2: sample_ID proportions per cluster ---
    props_sample = (
        cd138_adata.obs
        .groupby([key, "sample_ID"], observed=True)
        .size()
        .unstack(fill_value=0)
        .pipe(lambda df: df.div(df.sum(axis=1), axis=0))
    )
    props_sample.plot(kind="bar", stacked=True, ax=axes[1],
                      legend=True, width=0.8)
    axes[1].set_xlabel("Cluster")
    axes[1].set_ylabel("Proportion")
    axes[1].set_title("Sample composition")
    axes[1].tick_params(axis="x", rotation=90)
    axes[1].legend(bbox_to_anchor=(1.01, 1), loc="upper left",
                   fontsize=6, ncol=1)

    # --- Panel 3: disease_stage proportions per cluster ---
    props_stage = (
        cd138_adata.obs
        .groupby([key, "disease_stage"], observed=True)
        .size()
        .unstack(fill_value=0)
        .pipe(lambda df: df.div(df.sum(axis=1), axis=0))
    )
    stage_cols_ordered = [s for s in ["NBM", "SMM", "MM"]
                          if s in props_stage.columns]
    props_stage[stage_cols_ordered].plot(
        kind="bar", stacked=True, ax=axes[2],
        color=[stage_colors[s] for s in stage_cols_ordered],
        legend=True, width=0.8,
    )
    axes[2].set_xlabel("Cluster")
    axes[2].set_ylabel("Proportion")
    axes[2].set_title("Disease stage composition")
    axes[2].tick_params(axis="x", rotation=90)
    axes[2].legend(bbox_to_anchor=(1.01, 1), loc="upper left")

    plt.tight_layout()
    out_path = FIGURES_DIR / f"composition_leiden_res{res}.pdf"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  res={res}: {n_clusters} clusters -> {out_path.name}")

print("\nDone.")
