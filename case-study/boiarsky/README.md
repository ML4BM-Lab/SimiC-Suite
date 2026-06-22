# SimiC-Suite case study: Boiarsky 2022 NBM, SMM, MM plasma-cell scRNA-seq

This folder applies the SimiC-Suite workflow to the dataset reported in:

> Boiarsky, R., Haradhvala, N.J., Alberge, J.B. et al.
> *Single cell characterization of myeloma and its precursor conditions
> reveals transcriptional signatures of early tumorigenesis.*
> Nat Commun 13, 7040 (2022). https://doi.org/10.1038/s41467-022-33944-z

It mirrors the structure of `../case_study/` (the GSE145977 worked example
distributed with the SimiC-Suite Application Note) so the two analyses are
directly comparable.

Cohort (after dropping MGUS): 29 samples, ~19,100 CD138+ plasma cells,
three disease stages.

- NBM (normal bone marrow donors), n = 9 samples, ~9,000 cells
- SMM, n = 12 samples, ~8,400 cells
- MM, n = 8 samples, ~1,700 cells

The full published cohort also includes 6 MGUS samples (~850 cells). We
exclude MGUS in this run because the Boiarsky paper reports ~73% normal
plasma cells in MGUS samples on average (mean tumor purity ~27%). A
stage-level SimiC label of "MGUS" would therefore conflate MGUS-derived
normal PCs with MGUS tumor cells. Excluding MGUS also mirrors the
3-stage split used in the existing `../case_study/` (GSE145977) run, so
the two analyses are directly comparable.

Compared with the existing case study, this dataset reflects a different
cohort (Dana-Farber / Ghobrial lab) and uses the QC and preprocessing
protocol defined in the Boiarsky paper Methods rather than the
SimiC-Suite-internal protocol used for GSE145977.

## Workflow

```
01_data_prep.R                    R / Seurat
   |
   v writes data_GSE193531/
                  GSE193531_umi-count-matrix.csv.gz
                  GSE193531_cell-level-metadata.csv.gz
             seurat_objects/
                  GSE193531_with_qc.rds
                  GSE193531_filtered.rds
                  GSE193531_final.rds
             output/
                  singlecell_matrix.mtx
                  cell_ids.txt
                  genes_ids.txt
                  metadata.csv
                  QC_RNA_by_stage.pdf
                  QC_RNA_by_stage_postsubset.pdf
                  UMAP_by_stage.pdf

02_simic_preprocessing.ipynb      Python / Jupyter
   |
   v writes NBM-SMM-MM/
                  magic_output/magic_imputed.{pickle,h5ad}
                  inputFiles/expression_matrix.pickle
                  inputFiles/TF_final.csv
                  inputFiles/stage_final_annotation.csv

03_run_simic.py                   Python
   |
   v writes NBM-SMM-MM/outputSimic/matrices/NBM-SMM-MM/
                  *_simic_matrices.pickle
                  *_simic_matrices_filtered_BIC.pickle
                  *_wAUC_matrices_filtered_BIC.pickle
                  *_wAUC_matrices_filtered_BIC_collected.csv

04 (interactive)                  R / SimiCviz, in your own session
```

Step 4 is run interactively rather than as a fixed script, so figures and
TFs of interest can be explored against the paper's signatures.

## Step-by-step

### Step 1. Data preparation in R (`01_data_prep.R`)

- Downloads the GSE193531 processed-data CSVs via `GEOquery`:
  `GSE193531_umi-count-matrix.csv.gz` (cells x genes UMI counts) and
  `GSE193531_cell-level-metadata.csv.gz` (per-cell sample / stage / sex /
  batch annotation). The raw FASTQs are dbGaP-restricted
  (phs001323.v3.p1) and are not used.
- Builds a Seurat object directly from the count matrix. The published
  GEO matrix already reflects the paper's QC filter on Cell Ranger v2.0.1
  output aligned to hg38, but the script re-applies that filter to be
  transparent:
  - `percent.mt < 15`
  - `nFeature_RNA > 200`
  - `nCount_RNA < 50000`
  - `nFeature_RNA < 4000`
- Normalises with `LogNormalize`, `scale.factor = 1e4`, matching the
  paper's formula `e_g,c = log(1e4 / N_c * n_g,c + 1)`.
- Removes the gene sets the paper excludes prior to downstream analyses:
  - IGH / IGL / IGK genes (gene-symbol prefix `IG[HKL]`, plus `JCHAIN`).
  - Sex genes `XIST` and `RPS4Y1` (the two genes with the greatest
    absolute fold change between male and female samples in the
    Boiarsky cohort).
- Runs a Scanpy-style sanity-check UMAP with 14 PCs and Leiden
  resolution 1.5 (matching the paper's clustering parameters) and saves
  a UMAP coloured by disease stage.
- Exports the filtered raw count matrix and cell metadata in the format
  SimiCPipeline expects.

The script does not run scDblFinder. The Boiarsky paper handles doublets
and contaminating non-CD138+ cells through coarse Leiden clustering with
manual inspection, which the public GEO matrix already reflects.

Runtime: about 5-10 minutes (mostly the GEO download).

### Step 2. SimiC preprocessing in Python (`02_simic_preprocessing.ipynb`)

- Loads the matrix exported in step 1 as an AnnData object.
- Keeps the three stages NBM, SMM, MM (MGUS already removed in step 1).
- Runs MAGIC imputation through the `MagicPipeline` wrapper.
- Initialises an `ExperimentSetup`, restricting the TF candidate pool to
  the JASPAR2024 CORE vertebrate human TF list intersected with the count
  matrix, and selects the 100 most variable TFs and the 1,000 most
  variable target genes by Median Absolute Deviation (MAD).
- Saves the expression matrix, final TF list, and stage annotation
  required by SimiCPipeline.

The JASPAR TF file path defaults to `../case_study/data/JASPAR2024_Human_TFs.csv`
so the file does not need to be duplicated. Edit `TF_LIST_PATH` if your
copy lives elsewhere.

Runtime: about 5-10 minutes. MAGIC imputation dominates.

### Step 3. SimiCPipeline run (`03_run_simic.py`)

- Loads the expression matrix, TF list, and stage annotation written in
  step 2.
- Runs cross-validation over a grid of L1 and L2 penalties:
  - `LIST_L1 = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]`
  - `LIST_L2 = [1, 1e-1, 1e-2, 1e-3, 1e-4]`
- Selects the (lambda1, lambda2) pair maximising adjusted R-squared on
  held-out evaluation cells.
- Runs the final regression with the selected pair and produces
  stage-specific weighted incidence matrices for NBM, SMM, MM.
- Filters weights using a BIC criterion (variance threshold 0.9).
- Calculates AUC matrices for the filtered weights and writes a
  collected-AUC CSV.

Runtime: expect roughly comparable to the GSE145977 case study run, which
was 6h 23m for cross-validation and 14m for AUC calculation. Three
conditions, ~19k cells total. Per-stage cell counts: NBM ~9k, SMM ~8.4k,
MM ~1.7k. Treat the per-stage AUC distribution for MM with some caution
given its lower cell count (8 samples, ~1.7k cells), even though MM tumor
purity is very high (~98-100%) in this cohort.

### Step 4. Interactive visualization (your own session)

Reuse `../case_study/04_simic_visualization.R` as a starting point.
Update:

- `PROJECT_DIR` / `RUN_NAME` to `NBM-SMM-MM` (under the boiarsky/ folder,
  so the path differs from the existing case study even though the
  project sub-name is the same).
- Stage palette: same 3 stages as the existing case study.
- The `annotation_order` to `c("NBM", "SMM", "MM")` (matches existing).

The paper signatures most worth cross-checking against SimiC regulons:

- The "lost-in-abnormal" normal plasma cell signature (CD27 among top
  members).
- IFN-inducible signature (interferon response, shared across plasma
  cells, T cells, and monocytes in their data).
- Translocation-marker genes: `CCND1`, `MMSET/WHSC1`, `FGFR3`, `MAF`,
  `MAFB`, `CCND2`, `ITGB7`.

Be cautious when interpreting MM-specific regulons: the cohort has 8 MM
samples and roughly 1,700 MM cells, so any stage-specific finding should
be reported with its underlying cell count.

## Prerequisites

R (>= 4.3) packages:

- GEOquery, Seurat (>= 5), ggplot2, cowplot, Matrix, data.table
- SimiCviz (for step 4)

Python (>= 3.10) packages:

- simicpipeline (>= 0.1.0), anndata, pandas, numpy, scipy, scikit-learn
- magic-impute, scprep

No GTF reference is required. The Boiarsky case study identifies IG genes
by symbol prefix; the GSE145977 case study used a GRCh38 GTF to pick up
locus coordinates, which is unnecessary here because the GEO matrix is
indexed by HGNC gene symbol.

## Reproducing the figures

Boiarsky et al. do not run SimiC, so this case study is not a literal
reproduction of any specific figure. The intent is to apply the same
SimiC-Suite workflow to a cohort with an extensively documented paper, so
SimiC's stage-specific regulons can be cross-checked against the paper's
NMF signatures (Fig. 3 and Supplementary Table 8) and the abnormal-vs-
normal DEG list (Supplementary Tables 6 and 7).

## Notes

- The cohort metadata includes sex, age, sample preparation batch, and
  fresh/frozen status (in addition to stage). Stage is the only variable
  used as the SimiC phenotype in this run. The other covariates are
  retained in `metadata.csv` for use in downstream sanity checks.
- One patient was biopsied twice in the original cohort, at SMM (SMM-1)
  and after progression to MM (MM-8). They are treated as two
  independent samples here, the same way the paper treats them.
- All paths in the scripts are relative to `boiarsky/`. Edit the
  `setwd()` line at the top of the R script and the `os.chdir()` line
  in the Python script to point to your local copy of this folder.

## Citation

If you use this case study, please cite both:

- Boiarsky et al., Nat Commun 13, 7040 (2022) for the data.
- The SimiC-Suite Application Note for the workflow.
