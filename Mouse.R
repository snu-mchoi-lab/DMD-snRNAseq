######################
# 10x Single Cell Expression Data Analysis
# Author: Eun Young Jeon
# Date: 2022.01.11
######################

# This script performs basic single cell RNA sequencing analysis.
# The analysis includes:
# 1. Installing required packages
# 2. Setting working directory
# 3. Loading datasets
# 4. Merging datasets
# 5. Preprocessing (QC, scaling, PCA)
# 6. Clustering and UMAP
# 7. Marker identification and visualization
# 8. Defining subsets of interest

# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(magrittr)
library(cowplot)
library(harmony)

# 1. Set Working Directory
setwd("/Users/eunyoungjeon/Desktop/FGL/DMD_mouse/Mouse_work/10X")

# 2. Load 10x data
WT6.data <- Read10X(data.dir = "/mnt/workspace/Project_EJ/DMD/WT_6wk/outs/filtered_feature_bc_matrix")
WT28.data <- Read10X(data.dir = "/mnt/workspace/Project_EJ/DMD/WT_28wk/outs/filtered_feature_bc_matrix")
DMD6.data <- Read10X(data.dir = "/mnt/workspace/Project_EJ/DMD/D011-01/outs/filtered_feature_bc_matrix")
DMD6s.data <- Read10X(data.dir = "/mnt/workspace/Project_EJ/DMD/DMD_6wk_steroid_cellranger_output/outs/filtered_feature_bc_matrix")
DMD28.data <- Read10X(data.dir = "/mnt/workspace/Project_EJ/DMD/D003/D003-01/outs/filtered_feature_bc_matrix")
DMD28s.data <- Read10X(data.dir = "/mnt/workspace/Project_EJ/DMD/D003/D003-04/outs/filtered_feature_bc_matrix")

# 3. Create Seurat objects for each dataset
WT6 <- CreateSeuratObject(counts = WT6.data, project = "WT6", min.cells = 3, min.features = 200)
WT28 <- CreateSeuratObject(counts = WT28.data, project = "WT28", min.cells = 3, min.features = 200)
DMD6 <- CreateSeuratObject(counts = DMD6.data, project = "DMD6", min.cells = 3, min.features = 200)
DMD6s <- CreateSeuratObject(counts = DMD6s.data, project = "DMD6s", min.cells = 3, min.features = 200)
DMD28 <- CreateSeuratObject(counts = DMD28.data, project = "DMD28", min.cells = 3, min.features = 200)
DMD28s <- CreateSeuratObject(counts = DMD28s.data, project = "DMD28s", min.cells = 3, min.features = 200)

# Assign conditions to each dataset
WT6$stim <- "WT"
WT28$stim <- "WT"
DMD6$stim <- "DMD"
DMD6s$stim <- "Steroid"
DMD28$stim <- "DMD"
DMD28s$stim <- "Steroid"

# 4. Merge datasets into one Seurat object
alldata <- merge(WT6, y = c(WT28, DMD6, DMD6s, DMD28, DMD28s), 
                 add.cell.ids = c("WT_6wk", "WT_28wk", "DMD6", "DMD6s", "DMD28", "DMD28s"), 
                 project = "DMD_mouse_steroid_n6")

# 5. Preprocessing
# 5.1 QC and mitochondrial content filtering
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^mt-")
alldata[["percent.rp"]] <- PercentageFeatureSet(alldata, pattern = "^Rp[sl]")

VlnPlot(alldata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)

# Exclude cells with high mitochondrial counts and extreme gene counts
alldata <- subset(alldata, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

# 5.2 Normalize the data
alldata <- NormalizeData(alldata, normalization.method = "LogNormalize", scale.factor = 10000)

# 5.3 Identify highly variable features
alldata <- FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 2000)

# 5.4 Scale the data
all.genes <- rownames(alldata)
alldata <- ScaleData(alldata, features = all.genes)

# 6. Perform PCA
alldata <- RunPCA(alldata, features = VariableFeatures(object = alldata))
ElbowPlot(alldata)

# 7. Cluster cells using Harmony and UMAP
alldata <- RunHarmony(alldata, "orig.ident")
alldata <- FindNeighbors(alldata, dims = 1:8, reduction = 'harmony')
alldata <- FindClusters(alldata, resolution = 0.07, reduction = 'harmony')
alldata <- RunUMAP(alldata, dims = 1:7, reduction = 'harmony')

# Visualization
DimPlot(alldata, reduction = "umap", pt.size = 0.3)
DimPlot(alldata, reduction = "umap", group.by = "orig.ident", pt.size = 0.3)

# 8. Identify markers for clusters
markers <- FindAllMarkers(alldata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
