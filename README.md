# Code & Software Submission

The package layout is:

```
├── README.md                 <- you are reading it
├── code/                     <- production scripts
│   ├── 1.Seurat.R            <- QC, doublet removal, harmony integration, clustering
│   ├── 2.DEGene.R            <- per-cluster DE markers (FindAllMarkers)
│   ├── DEGene_lib.R          <- helper functions for script 2
│   ├── 3.scDiffer.R          <- between-group differential expression
│   ├── scDiffer_lib.R        <- helper functions for script 3
│   ├── scDiffer_config.yaml  <- group definitions & thresholds for script 3
│   ├── 4.Cellstate.R         <- gene-set activity scoring / cell-state annotation
│   └── Cellstate_lib.R       <- helper functions for script 4
└── demo/                     <- example dataset + runnable end-to-end demo
    ├── run_demo.sh           <- one-shot demo runner (4 steps)
    ├── data/                 <- small synthetic 10X dataset (2 samples x 1000 cells)
    │   ├── CK/               <- 10X folder of the control sample
    │   ├── As/               <- 10X folder of the treated sample
    │   └── doublet.xls       <- predicted doublet barcodes to be filtered
    └── results/              <- populated when demo is run
        ├── 1.Seurat/
        ├── 2.DEGene/
        ├── 3.scDiffer/
        └── 4.Cellstate/
```

------------------------------------------------------------------------

## 1. System requirements

**Operating system:**

* CentOS 7.9.2009

**Software:**

* **R ≥ 4.3.2**

**R packages:**

All analysis is implemented in **R 4.3.2** with the following R packages
(`code/README.md` lists them in detail):

| Package         | Version (tested)   | Role                                                 |
|-----------------|--------------------|------------------------------------------------------|
| Seurat          | 5.1.0              | QC, integration, clustering, DE                      |
| harmony         | 1.2.0              | Multi-sample integration                             |
| dplyr           | 1.1.3              | Data wrangling                                       |
| ggplot2         | 3.5.2              | Plotting                                             |
| patchwork       | 1.2.0              | Multi-panel plots                                    |
| cowplot         | 1.1.1              | Multi-panel plots                                    |
| AUCell          | 1.22.0             | Gene set activity scoring                            |
| BiocParallel    | 1.34.2             | Parallel computing backend for Bioconductor workflows|
| GSEABase        | 1.62.0             | Gene set data structure                              |
| GSVA            | 1.48.2             | Gene set activity scoring                            |
| UCell           | 2.4.0              | Gene set activity scoring                            |
| extrafont       | 0.19               | Font management for publication figures              |
| extrafontdb     | 1.0                | Font database backend for extrafont                  |
| Rttf2pt1        | 1.3.12             | TTF font converter (extrafont dependency)            |
| Cairo           | 1.6.2              | High-quality anti-aliased vector/raster figure export|
| showtext        | 0.9.7              | Flexible font rendering, support for CJK characters  |

------------------------------------------------------------------------

## 2. Installation guide

```bash
# 1. Install R 4.3.2 (skip if already present; conda works too)
#    See https://cran.r-project.org/

# 2. Install the R packages listed in §1. We recommend Bioconductor:
R -e 'install.packages("BiocManager")'
R -e 'BiocManager::install(c(
        "Seurat", "harmony", "dplyr", "ggplot2", "patchwork", "cowplot",
        "AUCell", "BiocParallel", "GSEABase", "GSVA", "UCell",
        "extrafont", "extrafontdb", "Rttf2pt1", "Cairo", "showtext"))'
```

**Typical install time.** On a "normal" desktop (8-core, 16 GB, fresh
conda env):

| Step                        | Time            |
|-----------------------------|-----------------|
| R + BiocManager install     | 5–10 min        |
| `BiocManager::install(...)` | 30–45 min       |
| **Total**                   | **≈ 35–55 min** |

------------------------------------------------------------------------

## 3. Demo

### 3.1 Dataset of demo

The `demo/data/` directory ships a **simulated** 10X dataset:

* 2 conditions × 1 replicate × 1000 cells = **2000 cells**
* **80,220 genes** per sample
* 2 ground-truth groups: **CK** (control) and **As** (treated)
* The dataset is provided in 10X v2 format (gzipped):

  * `demo/data/CK/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}`
  * `demo/data/As/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}`
  * `demo/data/doublet.xls` — 20 barcodes flagged as doublets, removed
    by script 1 during QC.

### 3.2 Run on data

The demo runner executes the production scripts in place, so first copy
them next to the data, then run:

```bash
cd demo
cp ../code/*.R ../code/*.yaml .
bash run_demo.sh
```

The script runs four R scripts in sequence and writes their outputs
to `demo/results/`:

| Step                                             | Script              | R env   | Output dir                 |
|--------------------------------------------------|---------------------|---------|----------------------------|
| Seurat filter + harmony integration + clustering | `1.Seurat.R`        | R 4.3.2 | `demo/results/1.Seurat/`   |
| Find DE genes per cluster                        | `2.DEGene.R`        | R 4.3.2 | `demo/results/2.DEGene/`   |
| Group-level differential expression (CK vs As)   | `3.scDiffer.R`      | R 4.3.2 | `demo/results/3.scDiffer/` |
| Gene set activity scoring (e.g. AddModuleScore)  | `4.Cellstate.R`     | R 4.3.2 | `demo/results/4.Cellstate/`|

### 3.3 Expected output

```
demo/results/
├── 1.Seurat/
│   ├── UMAP.pdf              <- UMAP coloured by sample / cluster
│   ├── TSNE.pdf              <- t-SNE coloured by sample / cluster
├── 2.DEGene/
│   ├── DeGene.list.xls       <- per-cluster positive marker table
├── 3.scDiffer/
│   ├── AllGene.avg_exp.Ep2.xls        <- grouped average-expression table of ALL genes for each cluster  
│   ├── DifferMarker.CK-vs-As.all.xls  <- DE table for all tested genes
│   ├── DifferMarker.CK-vs-As.diff.xls <- DE table filtered at |log2FC| > 0.36 & P < 0.05
│   ├── Volcano.CK-vs-As_Ep2.pdf       <- volcano plot per cluster
│   └── Heatmap.CK-vs-As.pdf           <- top-5 DE genes heatmap per cluster
└── 4.Cellstate/
    ├── ref.yaml              <- gene sets used for scoring (one entry per state)
    └── AddModuleScore/
        └── Cells.CellSelect.list.xls <- per-cell gene-set scores & Positive/Negative calls
```

------------------------------------------------------------------------

## 4. Instructions for use

### 4.1 QC, integration and clustering on your own 10X data

`code/1.Seurat.R` is the entry script:

1. Place your 10X folders under the working directory, e.g.
   `my_data/Control/`, `my_data/Treated/`. Each subfolder must contain
   `barcodes.tsv.gz`, `genes.tsv.gz`, `matrix.mtx.gz` (10X v2 format).
2. Edit `data_name <- c("CK","As")` at the top of the script to your
   sample names.
3. Adjust the QC thresholds in the `subset()` call
   (`nFeature_RNA < 98000`, `nCount_RNA < 40000`) to your dataset.
4. Prepare a `doublet.xls` file listing predicted doublet barcodes
   (one per line); these cells are removed before integration.
5. Run `Rscript code/1.Seurat.R`.

Outputs: `TSNE.pdf`, `UMAP.pdf`, and the integrated Seurat object
`obj.Rda` used by all downstream scripts.

### 4.2 Per-cluster marker genes

`code/2.DEGene.R` loads `obj.Rda`, runs `FindAllMarkers`
(`only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25`) and writes
`DeGene.list.xls` plus `markers.Rda`. Adjust these thresholds at the top
of the script if needed, then run `Rscript code/2.DEGene.R`.

### 4.3 Between-group differential expression

`code/3.scDiffer.R` is driven by `code/scDiffer_config.yaml`:

* `group.by_sample` — map each sample to its condition;
* `differ` — list of group pairs to contrast, e.g. `[CK, As]`;
* `FindMarkers` — significance thresholds: `log2fc: 0.36`, `min_pct: 0.1`,
  `pvalue: 0.05`, `method: MAST`;
* `plots` — `top: 5` genes per cluster in the heatmap plots.

Edit the YAML, then run `Rscript code/3.scDiffer.R`. Outputs include the
full and filtered DE tables (`DifferMarker.<pair>.all.xls` /
`.diff.xls`), average-expression matrices, volcano plots and a heatmap.

### 4.4 Gene-set activity scoring / cell-state annotation

`code/4.Cellstate.R` scores each cell against user-supplied gene sets:

```bash
Rscript code/4.Cellstate.R <method> <outdir>
```

* `<method>` — one of `mean`, `AddModuleScore` (Seurat default),
  `AUCell`, or `GSVA` (ssgsea, Poisson kcdf, suited to UMI counts);
* `<outdir>` — results are written to this directory;
* gene sets are read from `ref.yaml` (a YAML file mapping each
  state name to a list of gene IDs; genes absent from the object are
  dropped automatically).

Per-cell scores are attached to the object, each cell is assigned to the
highest-scoring state (`CellState`), and `Cells.CellSelect.list.xls`
reports all scores plus Positive/Negative calls at the 75th-percentile
cutoff of each gene set.
