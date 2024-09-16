
```r
---
title: "PBMC ATAC-seq Analysis"
author: "Vidya Sai Rashmitha Reddy Yellareddy"
output: .R
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

# Introduction

This analysis performs chromatin accessibility analysis on PBMC ATAC-seq data using `Signac` and `Seurat` packages. We will go through the steps of data loading, quality control, dimensional reduction, and clustering.

# 1. Load Required Libraries

```{r libraries}
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(tidyverse)
library(hdf5r)
```

# 2. Download Data

We use the system commands to download the dataset files from 10x Genomics.

```{r download-data, eval=FALSE}
# Downloading files using curl
system("curl -O https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
system("curl -O https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv")
system("curl -O https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz")
system("curl -O https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi")
```

# 3. Load Data

We load the downloaded data into R.

```{r load-data}
# Load fragment file
frag.file <- read.delim('/Users/rashmithayellareddy/atac_v1_pbmc_10k_fragments.tsv.gz', header = FALSE, nrows = 10)
head(frag.file)

# Load the filtered peak matrix
counts <- Read10X_h5('/Users/rashmithayellareddy/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5')
counts[1:10, 1:10]

# Create Chromatin Assay
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "/Users/rashmithayellareddy/atac_v1_pbmc_10k_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)
str(chrom_assay)

# Load metadata
metadata <- read.csv('/Users/rashmithayellareddy/atac_v1_pbmc_10k_singlecell.csv', header = TRUE, row.names = 1)
View(metadata)
```

# 4. Create a Seurat Object

We now create a Seurat object using the chromatin assay and metadata.

```{r create-seurat-object}
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = metadata,
  assay = 'ATAC'
)
str(pbmc)
```

# 5. Gene Annotation

Add gene annotations to the Seurat object from the Ensembl database.

```{r gene-annotation}
# Extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# Change to UCSC style
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))

# Add gene information to the object
Annotation(pbmc) <- annotations
pbmc@assays$ATAC@annotation
```

# 6. Quality Control

We compute the quality control metrics, such as nucleosome signal, TSS enrichment, blacklist ratio, and the percentage of reads in peaks.

```{r quality-control}
# Compute nucleosome signal
pbmc <- NucleosomeSignal(pbmc)

# Compute TSS enrichment
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# Add blacklist ratio and pct reads in peaks
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

View(pbmc@meta.data)
```

# 7. Visualizing QC

We use density scatter plots and violin plots to visualize quality control metrics.

```{r visualize-qc, fig.width=10, fig.height=6}
# Density Scatter plots
a1 <- DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
a2 <- DensityScatter(pbmc, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# Combine the plots
a1 | a2

# Violin plot for QC metrics
VlnPlot(object = pbmc, 
        features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal', 'blacklist_ratio', 'pct_reads_in_peaks'),
        pt.size = 0.1,
        ncol = 6)
```

# 8. Filtering Poor Quality Cells

We filter out low-quality cells based on QC metrics.

```{r filter-cells}
pbmc <- subset(x = pbmc,
               subset = nCount_ATAC > 3000 &
                 nCount_ATAC < 30000 &
                 pct_reads_in_peaks > 15 & 
                 blacklist_ratio < 0.05 &
                 nucleosome_signal < 4 &
                 TSS.enrichment > 3)
```

# 9. Normalization and Dimensional Reduction

We normalize the data using TF-IDF, perform dimensional reduction, and compute depth correlation.

```{r normalization}
pbmc <- RunTFIDF(pbmc)  # Normalization
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')  # Select top features
pbmc <- RunSVD(pbmc)  # Dimensionality reduction
DepthCor(pbmc)  # Compute depth correlation
```

# 10. Clustering and UMAP

Finally, we perform non-linear dimensional reduction and clustering using UMAP, followed by visualization with `DimPlot`.

```{r clustering}
# UMAP and clustering
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, algorithm = 3)

# UMAP plot
DimPlot(object = pbmc, label = TRUE) + NoLegend()
```

# Conclusion
In this analysis, we successfully processed and analyzed scATAC-seq data from 10x Genomics PBMC dataset using **Signac** and **Seurat**. Key steps in the workflow included:

- Downloading and pre-processing raw data, including fragments, metadata, and peak-barcode matrix files.
- Creating a Seurat object with chromatin accessibility data and incorporating gene annotations from **EnsDb.Hsapiens.v75** to link the chromatin data to gene regions.
- Performing quality control by calculating the nucleosome signal, TSS enrichment, blacklist ratio, and fraction of reads in peaks, followed by visualizing these metrics using **VlnPlot** and **DensityScatter** plots.
- Filtering out poor-quality cells based on stringent criteria such as nucleosome signal, TSS enrichment score, and the fraction of reads in peaks.
- Normalizing the data and performing linear and non-linear dimensionality reductions using TF-IDF, SVD, and UMAP.
- Identifying clusters in the PBMC dataset using UMAP for visualization and the Louvain algorithm for clustering.

The clustering results provided insights into the chromatin accessibility profiles of different PBMC subpopulations, which could be further explored in downstream analyses, such as motif enrichment or integration with scRNA-seq data for multi-omic analysis. This workflow demonstrates a robust approach to analyzing single-cell chromatin accessibility data, helping to uncover regulatory mechanisms and cellular heterogeneity within PBMCs.


```
