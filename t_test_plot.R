# Install and load the required packages
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
library(ggpubr)
library(Seurat)
library(purrr)
library(cowplot)

# Define custom levels for the factor
my_levels <- c("CTRL", "BMD", "DMD")
DMD_human_satellite_pro$stim <- factor(DMD_human_satellite_pro$stim, levels = my_levels)

# Function to generate violin plots and save them as PDFs
vp_case1 <- function(gene_signature, file_name, test_sign, y_max) {
  
  plot_case1 <- function(signature) {
    VlnPlot(DMD_human_satellite_pro, features = signature,
            pt.size = 0.01, 
            group.by = "stim", 
            cols = c("#a9a9a9", "#767676", "#333333"),
            y.max = y_max) + 
      stat_compare_means(comparisons = test_sign, label = "p.signif", method = "t.test")
  }
  
  # Generate plots for each gene signature
  plots <- purrr::map(gene_signature, plot_case1)
  cowplot::plot_grid(plotlist = plots)
  
  # Save the plot as a PDF
  ggsave(paste0(file_name, "_r.pdf"), width = 12, height = 8)
}

# Example gene signature and comparison list
gene_sig <- c("DOCK10", "DOCK8", "PTPRJ", "CAMK1D", "PARP8", "TNS3")
comparisons <- list(c("CTRL", "DMD"), c("CTRL", "BMD"), c("BMD", "DMD"))

# Generate and save the plots
vp_case1(gene_signature = gene_sig, file_name = "Figure2_DMD_human_myeloid", test_sign = comparisons, y_max = 7)


