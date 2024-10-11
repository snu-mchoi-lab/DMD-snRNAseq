# Load data and initial plot
DMD_human <- readRDS("DMD_human.rds")
DimPlot(DMD_human)

# Install necessary packages
BiocManager::install(version = 'devel')
BiocManager::install("GSVA", force = TRUE)
devtools::install_github("ncborcherding/escape")

# Load required library
library(escape)

# Retrieve Hallmark gene sets
GS.hallmark <- getGeneSets(library = "H")

# Define custom gene sets
gene.sets <- list(
  Necrosis = c("BAK1", "BAX", "BIRC2", "BIRC3", "CASP1", "CASP3", "CASP4", "CASP5", "CASP8", 
               "CDC37", "CFLAR", "CHMP2A", "CHMP2B", "CHMP3", "CHMP4A", "CHMP4B", "CHMP4C", 
               "CHMP6", "CHMP7", "CYCS", "ELANE", "FADD", "FAS", "FASLG", "FLOT1", "FLOT2", 
               "GSDMD", "GSDME", "GZMB", "HMGB1", "HSP90AA1", "IL18", "IL1A", "IL1B", "IRF1", 
               "IRF2", "ITCH", "MLKL", "OGT", "PDCD6IP", "PELI1", "PRKN", "RIPK1", "RIPK3", 
               "RPS27A", "SDCBP", "STUB1", "TNFRSF10A", "TNFRSF10B", "TNFSF10", "TP53", "TP63", 
               "TRADD", "TRAF2", "UBA52", "UBB", "UBC", "UBE2L3", "XIAP"),
  Myeloid = c("SPI1", "FCER1G", "CSF1R"),
  Tcells = c("CD3E", "CD3D", "CD3G", "CD7", "CD8A")
)

# Run ssGSEA with Hallmark gene sets
scRep_example <- runEscape(DMD_human, 
                           method = "ssGSEA",
                           gene.sets = GS.hallmark, 
                           groups = 5000, 
                           min.size = 0,
                           new.assay.name = "escape.ssGSEA")

# Run ssGSEA with custom gene sets
scRep_example2 <- runEscape(DMD_human, 
                            method = "ssGSEA",
                            gene.sets = gene.sets, 
                            groups = 5000, 
                            min.size = 0,
                            new.assay.name = "escape.ssGSEA")

# Normalize the results
scRep_example <- performNormalization(scRep_example, 
                                      assay = "escape.ssGSEA", 
                                      gene.sets = GS.hallmark, 
                                      scale.factor = scRep_example$nFeature_RNA)

scRep_example2 <- performNormalization(scRep_example2, 
                                       assay = "escape.ssGSEA", 
                                       gene.sets = gene.sets)

# Create heatmap enrichment plots
heatmapEnrichment(scRep_example, 
                  group.by = "ident",
                  assay = "escape.ssGSEA",
                  gene.set.use = GS.hallmark,
                  scale = TRUE,
                  cluster.rows = TRUE,
                  cluster.columns = TRUE)

heatmapEnrichment(scRep_example, 
                  group.by = "ident",
                  assay = "escape.ssGSEA",
                  gene.set = "HALLMARK-APOPTOSIS")

heatmapEnrichment(scRep_example2, 
                  group.by = "ident",
                  assay = "escape.ssGSEA",
                  gene.set = "Necrosis", 
                  facet.by = "stim")

# Create geyser enrichment plots
geyserEnrichment(scRep_example2, 
                 group.by = "ident",
                 assay = "escape.ssGSEA",
                 gene.set = "Necrosis", 
                 facet.by = "stim")

geyserEnrichment(scRep_example, 
                 assay = "escape.ssGSEA",
                 gene.set = "HALLMARK-APOPTOSIS")

geyserEnrichment(scRep_example, 
                 assay = "escape.ssGSEA",
                 gene.set = "HALLMARK-APOPTOSIS", 
                 facet.by = "stim")

geyserEnrichment(scRep_example, 
                 assay = "escape.ssGSEA",
                 gene.set = "HALLMARK-G2M-CHECKPOINT", 
                 facet.by = "stim")

geyserEnrichment(scRep_example, 
                 assay = "escape.ssGSEA",
                 gene.set = "HALLMARK-MYOGENESIS", 
                 facet.by = "stim")

geyserEnrichment(scRep_example, 
                 assay = "escape.ssGSEA",
                 gene.set = "HALLMARK-INTERFERON-ALPHA-RESPONSE", 
                 facet.by = "stim")

geyserEnrichment(scRep_example, 
                 assay = "escape.ssGSEA",
                 gene.set = "HALLMARK-TGF-BETA-SIGNALING", 
                 facet.by = "stim")

geyserEnrichment(scRep_example, 
                 assay = "escape.ssGSEA",
                 gene.set = "HALLMARK-INFLAMMATORY-RESPONSE", 
                 facet.by = "stim")

geyserEnrichment(scRep_example2, 
                 assay = "escape.ssGSEA",
                 gene.set = "Necrosis", 
                 facet.by = "stim")
