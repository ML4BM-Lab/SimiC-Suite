# SimiC-Suite

SimiC-Suite is a reproducible framework for phenotype-specific gene regulatory network (GRN) analysis from single-cell RNA-seq data.

The suite comprises two complementary software packages:

| Package | Language | Description |
|----------|----------|-------------|
| SimiCPipeline | Python | Phenotype-specific GRN inference and Regulon Activity Score (RAS) calculation |
| SimiCviz | R | Tool-agnostic RAS calculation, assessment and visualization of GRN outputs |

Together, these packages support reproducible GRN inference, RAS calculation, regulatory dissimilarity analysis, and downstream visualization across Python and R workflows.

---

## Overview

<p align="center">
  <img src="figures/simic_suite_graphical_abstract.png" alt="SimiC-Suite overview" width="1000">
</p>

---

## Components

### SimiCPipeline (Python)

📦 Repository: https://github.com/ML4BM-Lab/SimiCPipeline

SimiCPipeline is a modular Python framework for phenotype-specific GRN inference based on the previously published SimiC method (Peng et al., 2022).

Key functionality includes:

- Data preprocessing and imputation
- Automated hyperparameter tuning
- Phenotype-specific GRN inference
- Network filtering and quality assessment
- Regulon Activity Score (RAS) calculation
- Integration with the scverse ecosystem
- Containerized execution via Docker

---

### SimiCviz (R / Bioconductor)

📦 Repository: https://github.com/ML4BM-Lab/SimiCviz

SimiCviz is a tool-agnostic R package for the assessment and visualization of GRN outputs.

In addition to SimiCPipeline outputs, SimiCviz accepts outputs from other GRN inference methods that can be represented as TF-by-target weight matrices, including frameworks such as SCENIC and Pando.

Key functionality includes:

- Parallelized Regulon Activity Score (RAS) calculation
- Regulatory dissimilarity analysis
- Empirical cumulative distribution function (ECDF) visualization of RAS distributions
- TF-target weight exploration
- Regulatory activity heatmaps
- UMAP visualization
- Import and export utilities
- Integration with Bioconductor workflows

---

## Typical Workflows

### End-to-End SimiC Workflow

text Single-cell RNA-seq data           │           ▼ SimiCPipeline (Python)   • Preprocessing   • GRN inference   • RAS calculation           │           ▼ Phenotype-specific GRNs and RAS matrices           │           ▼ SimiCviz (R)   • Assessment   • Visualization   • Reporting           │           ▼ Biological insights and publication-ready outputs 

### Tool-Agnostic Workflow

text SCENIC / Pando / other GRN inference methods               │               ▼ TF–target weight matrices and expression data               │               ▼ SimiCviz (R)   • RAS calculation   • Network assessment   • Comparative analysis   • Visualization               │               ▼ RAS profiles, figures, and downstream analyses 

---

## Documentation

Please refer to the individual repositories for installation instructions, tutorials, and documentation.

- SimiCPipeline: https://github.com/ML4BM-Lab/SimiCPipeline
- SimiCviz: https://github.com/ML4BM-Lab/SimiCviz

---

## Citation

If you use SimiC-Suite in your research, please cite:

### SimiC-Suite

Citation will be added upon publication.

---

## License

Please refer to the license information provided in the individual repositories.
