# Update CellChat objects
cellchat_G <- updateCellChat(cellchat)
cellchat_mdx <- updateCellChat(cellchat_mdx)

# Merge CellChat objects
object.list <- list(mdx = cellchat_mdx, G = cellchat_G)
cellchat_G_mdx <- mergeCellChat(object.list, add.names = names(object.list))

# Compare interactions
gg1 <- compareInteractions(cellchat_G_mdx, show.legend = FALSE, group = c(1, 2))
gg2 <- compareInteractions(cellchat_G_mdx, show.legend = FALSE, group = c(1, 2), measure = "weight")
gg1 + gg2

# Visualize differential interactions
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat_G_mdx, weight.scale = TRUE)
netVisual_diffInteraction(cellchat_G_mdx, weight.scale = TRUE, measure = "weight")

# Heatmap of interactions
gg1 <- netVisual_heatmap(cellchat_G_mdx)
gg2 <- netVisual_heatmap(cellchat_G_mdx, measure = "weight")
gg1 + gg2

# Visualize circle plots of interaction counts
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
par(mfrow = c(1, 2), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_circle(
    object.list[[i]]@net$count, 
    weight.scale = TRUE, 
    label.edge = FALSE, 
    edge.weight.max = weight.max[2], 
    edge.width.max = 12, 
    title.name = paste0("Number of interactions - ", names(object.list)[i])
  )
}

# Identify and visualize signaling networks based on functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

# Compare outgoing signaling with heatmaps
library(ComplexHeatmap)
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i + 1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# Visualize specific pathway (RESISTIN) across datasets
pathways.show <- c("RESISTIN")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

par(mfrow = c(1, 2), xpd = TRUE)
for (i in seq_along(object.list)) {
  netVisual_aggregate(
    object.list[[i]], 
    signaling = pathways.show, 
    layout = "circle", 
    edge.weight.max = weight.max[1], 
    edge.width.max = 10, 
    signaling.name = paste(pathways.show, names(object.list)[i])
  )
}
