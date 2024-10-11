# Install necessary packages (only once, no need to repeat)
BiocManager::install(c("slingshot", "tradeSeq", "fgsea", "gridExtra", "scater"))
install.packages("/mnt/workspace/data/podo/Projects/Project_EJ/DMD/1.1.0.tar.gz", repos = NULL, type = "source")

# Load libraries
library(scater)
library(Seurat)
suppressPackageStartupMessages({
  library(slingshot)
  library(SingleCellExperiment)
  library(RColorBrewer)
  library(scales)
  library(viridis)
  library(UpSetR)
  library(pheatmap)
  library(msigdbr)
  library(fgsea)
  library(knitr)
  library(ggplot2)
  library(gridExtra)
  library(tradeSeq)
  library(DelayedMatrixStats)
})

# Read and subset data
DMD_human <- readRDS("DMD_human.rds")
DMD_human_sate_muscle <- subset(DMD_human, subset = seurat_clusters %in% c("0", "1", "3"))
DimPlot(DMD_human_sate_muscle)

DMD_human_n3_sate_muscle <- subset(DMD_human_sate_muscle, subset = orig.ident %in% c("Ctrl1", "BMD1", "DMD1"))
table(DMD_human_n3_sate_muscle$orig.ident, DMD_human_n3_sate_muscle$seurat_clusters)

# Convert to SingleCellExperiment and run Slingshot
sce_human_n3 <- as.SingleCellExperiment(DMD_human_n3_sate_muscle, assay = "RNA")
sce_human_n3 <- slingshot(sce_human_n3, reducedDim = 'UMAP', clusterLabels = colData(sce_human_n3)$seurat_clusters, start.clus = '5', approx_points = 150)
plotUMAP(sce_human_n3, colour_by = "slingPseudotime_1")

# Fit generalized additive model (GAM)
sce_human_n3 <- fitGAM(sce_human_n3, conditions = factor(colData(sce_human_n3)$stim), nknots = 5)

# Association test and significant gene identification
assocRes <- associationTest(sce_human_n3, lineages = TRUE, l2fc = log2(2))
mockGenes <- rownames(assocRes)[which(p.adjust(assocRes$pvalue_lineage1_conditionCtrl1, "fdr") <= 0.05)]
DMDGenes <- rownames(assocRes)[which(p.adjust(assocRes$pvalue_lineage1_conditionDMD1, "fdr") <= 0.05)]

# Visualization using UpSet plot
UpSetR::upset(fromList(list(mock = mockGenes, DMD = DMDGenes)))

# Heatmap of smoothed predictions
yhatSmooth <- predictSmooth(sce_human_n3, gene = mockGenes, nPoints = 50, tidy = FALSE)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE)
cl <- sort(cutree(heatSmooth$tree_row, k = 6))
table(cl)

# FGSEA analysis using GO terms
geneSets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
geneSets <- geneSets[geneSets$gene_symbol %in% names(sce_human_n3), ]
m_list <- split(geneSets$gene_symbol, geneSets$gs_name)
stats <- assocRes$waldStat_lineage1_conditionWT
stats_naomit <- na.omit(stats)
eaRes <- fgsea(pathways = m_list, stats = stats_naomit, nperm = 5e4, minSize = 10)

# Differential expression analysis
condRes <- conditionTest(sce_human_n3, l2fc = log2(2), global = FALSE, pairwise = TRUE)
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
oo <- order(condRes$waldStat, decreasing = TRUE)

# Plot smoothers for specific genes
plotSmoothers(sce_human_n3, assays(sce_human_n3)$counts, gene = "Ezh2") + ggtitle("Ezh2")
plotSmoothers(sce_human_n3, assays(sce_human_n3)$counts, gene = "Dcn") + ggtitle("Dcn")

# FGSEA result visualization
kable(head(eaRes[order(eaRes$pval), 1:3], n = 20))
write.csv(eaRes, file = "240107_Tradeseq_eaRes_wt_DMD_GPSP.csv")
