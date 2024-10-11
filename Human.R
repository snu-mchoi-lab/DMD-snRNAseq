##############################
# 10x Single Cell RNA-Seq Analysis
# Author: Eun Young Jeon
# Date: 2020.02.13
# This script performs basic single-cell RNA-seq analysis
##############################

# Required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(magrittr)
library(cowplot)
library(harmony)
library(rliger)
library(SeuratWrappers)
library(velocyto.R)
library(Rcpp)

# Install additional packages if necessary
remotes::install_github('satijalab/seurat-wrappers')
devtools::install_github("velocyto-team/velocyto.R")
BiocManager::install(Rcpp)

# 1. Set Working Directory
setwd("/mnt/workspace/10x_temp_EJ")

# 2. Load 10x Data for Multiple Samples
data_list <- list(
  C1 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD001/filtered_feature_bc_matrix/"),
  C2 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD002/filtered_feature_bc_matrix/"),
  BMD1 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD003/filtered_feature_bc_matrix/"),
  DMD1 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD005/filtered_feature_bc_matrix/"),
  DMD2 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD007/outs/filtered_feature_bc_matrix/"),
  BMD2 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD008/outs/filtered_feature_bc_matrix/"),
  BMD3 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD009/outs/filtered_feature_bc_matrix/"),
  DMD3 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD010/outs/filtered_feature_bc_matrix/"),
  C3 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD011/MD011/outs/filtered_feature_bc_matrix/"),
  C4 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD012/MD012/outs/filtered_feature_bc_matrix/"),
  C5 = Read10X(data.dir = "/mnt/workspace/10x_temp_EJ/MD013/outs/filtered_feature_bc_matrix/")
)

# 3. Create Seurat Objects
seurat_objects <- lapply(names(data_list), function(name) {
  CreateSeuratObject(counts = data_list[[name]], project = name, min.cells = 3, min.features = 200)
})

# 4. Add Metadata (Stimulus Information)
stim_labels <- c("CTRL", "CTRL", "BMD", "DMD", "DMD", "BMD", "BMD", "DMD", "CTRL", "CTRL", "CTRL")
for (i in seq_along(seurat_objects)) {
  seurat_objects[[i]]$stim <- stim_labels[i]
}

# 5. Calculate Percent Mitochondrial and Ribosomal Content
seurat_objects <- lapply(seurat_objects, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  return(obj)
})

# 6. Subset Data Based on Quality Control Metrics
seurat_objects <- lapply(seurat_objects, function(obj) {
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
})

# 7. Normalize and Identify Variable Features
seurat_objects <- lapply(seurat_objects, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  return(obj)
})

# 8. Integrate Data Across Samples
DMD_anchors <- FindIntegrationAnchors(object.list = seurat_objects, dims = 1:20)
DMD_combined <- IntegrateData(anchorset = DMD_anchors, dims = 1:20)

# Set the default assay to "integrated"
DefaultAssay(DMD_combined) <- "integrated"

# 9. Scale Data
DMD_combined <- ScaleData(DMD_combined, split.by = "orig.ident", do.center = FALSE)

# 10. Run ALS and Quantile Normalization
DMD_combined <- RunOptimizeALS(DMD_combined, k = 20, lambda = 5, split.by = "orig.ident")
DMD_combined <- RunQuantileNorm(DMD_combined, split.by = "orig.ident")

# 11. Find Neighbors and Clusters
DMD_combined <- FindNeighbors(DMD_combined, reduction = "iNMF", dims = 1:20, graph.name = "test")
DMD_combined <- FindClusters(DMD_combined, resolution = 0.03, graph.name = "test")

# 12. Run UMAP for Visualization
DMD_combined <- RunUMAP(DMD_combined, dims = 1:ncol(DMD_combined[["iNMF"]]), reduction = "iNMF")

# 13. Plot UMAP
DimPlot(DMD_combined, reduction = "umap", pt.size = 0.1)
