# Install Bioconductor and required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AUCell", "RcisTarget", "GENIE3", "SingleCellExperiment"))

# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Install SCopeLoomR
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

# Install lazy loading packages
packages_lazyinstall <- c("zoo", "mixtools", "rbokeh", "DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne", "doMC", "doRNG")
install.packages(packages_lazyinstall)

# Load necessary libraries
library(zoo)
library(mixtools)
library(rbokeh)
library(DT)
library(NMF)
library(ComplexHeatmap)
library(R2HTML)
library(Rtsne)
library(doMC)
library(doRNG)
library(SeuratDisk)
library(SCopeLoomR)
library(hdf5r)
library(feather)
library(SCENIC)

# Install Seurat v3.1.4 for compatibility
remotes::install_version(package = "Seurat", version = "3.1.4")

# Load data and subset
DMD_mouse_GS_n9 <- readRDS("DMD_mouse_GS_n9.rds")
DMD_mouse_GS_n9_satellite <- subset(DMD_mouse_GS_n9, subset = seurat_clusters == "5")
DimPlot(DMD_mouse_GS_n9_satellite)

# Subset data for SP-stimulated cells and plot
mdx_n9 <- subset(DMD_mouse_GS_n9, subset = stim == "SP")
DimPlot(mdx_n9)

# Convert Seurat object to loom format and extract expression matrix and cell info
mdx_n9_loom <- as.loom(mdx_n9, filename = "mdx_n9.loom", verbose = FALSE)
exprMat <- get_dgem(mdx_n9_loom)
cellInfo <- get_cell_annotation(mdx_n9_loom)

# Initialize SCENIC settings
org <- "mgi"
dbDir <- "cisTarget_databases"
myDatasetTitle <- "SCENIC on DMD Mouse"
scenicOptions <- initializeScenic(org = org, dbDir = dbDir, datasetTitle = myDatasetTitle, nCores = 10)

# Gene filtering
genesKept <- geneFiltering(exprMat, scenicOptions = scenicOptions, minCountsPerGene = 3 * .01 * ncol(exprMat), minSamples = ncol(exprMat) * .01)
exprMat_filtered <- exprMat[genesKept, ]

# Run correlation and GENIE3 for network inference
runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions)

# Run SCENIC steps
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log = log2(exprMat + 1), skipTsne = TRUE, skipHeatmap = TRUE)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

# Export to loom and save SCENIC options
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file = "int/scenicOptions.Rds")

# Regulon specificity score (RSS) calculation
rss <- calcRSS(AUC = getAUC(loadInt(scenicOptions, "aucell_regulonAUC")), cellAnnotation = cellInfo$seurat_clusters)
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

# Heatmap visualization
Heatmap(rss, name = "Regulon activity")

# Export final results
saveRDS(cellInfo, file = "int/cellInfo_mdx.rds")
