library(CellChat)
library(Seurat)
library(ggplot2)                  
library(patchwork)
library(igraph)

setwd("~/Documents/cell_Chat_Comparison/")
#Rds file contains the clustered Seurat object for single nuclear data. Run cellchat analysis separately for both control and mutant. 
all_cell <- readRDS("1.clustering_all_cell.rds")

fname = "control"
control <- subset(all_cell, sample == "CTR")


ctr.data <- GetAssayData(control, assay = "SCT", slot = "data") # normalized data matrix
ctr.label <- data.frame("labels" = control@meta.data$cell.type, row.names = colnames(ctr.data))
ctr.meta <- data.frame("labels" = ctr.label$labels, row.names = rownames(ctr.label)) # create a dataframe of the cell labels
ctr.cc <- createCellChat(object = ctr.data, meta = ctr.meta, group.by = "labels")


#Run cell chat for mutant only first. 
dir.create(fname)
CellChatDB <- CellChatDB.mouse 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
ctr.cc@DB <- CellChatDB.use
ctr.cc <- subsetData(ctr.cc)
future::plan("multisession", workers = 4) # do parallel
ctr.cc <- identifyOverExpressedGenes(ctr.cc)
ctr.cc <- identifyOverExpressedInteractions(ctr.cc)
ctr.cc <- computeCommunProb(ctr.cc)
ctr.cc <- filterCommunication(ctr.cc, min.cells = 5)
ctr.cc <- computeCommunProbPathway(ctr.cc)
ctr.cc <- aggregateNet(ctr.cc)

groupSize <- as.numeric(table(ctr.cc@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(ctr.cc@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(ctr.cc@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- ctr.cc@net$weight
par(mfrow = c(2,3), xpd=TRUE)
dir.create(paste0(fname, "/Sep_Interactions"))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  pdf(paste0(fname, "/Sep_Interactions/",rownames(mat)[i], ".pdf"))
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  dev.off()
}


pathways.show.all <- ctr.cc@netP$pathways

dir.create(paste0(fname, "/Results"))
dir.create(paste0(fname, "/Results/LR"))

cell.ident <- levels(ctr.cc@idents)
vertex.receiver = c(1,2,5,6)

for(i in 1:length(pathways.show.all)){
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  setwd(paste0("~/Documents/cell_Chat_Comparison/", fname, "/Results"))
  
  netVisual(ctr.cc, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  
  pairLR <- extractEnrichedLR(ctr.cc, signaling = pathways.show.all[i], geneLR.return = FALSE)
  
  setwd("~/Documents/cell_Chat_Comparison/")
  dir.create(paste0(fname, "/Results/", pathways.show.all[i]))
  
  for(lr in 1:nrow(pairLR)){
    LR.show <- pairLR[lr,] # show one ligand-receptor pair
    # Hierarchy plot
    pdf(paste0(fname, "/Results/", pathways.show.all[i], "/", LR.show, ".pdf"))
    par(mfrow = c(1,3), xpd=TRUE)
    netVisual_individual(ctr.cc, signaling = pathways.show.all[i],  pairLR.use = LR.show, vertex.receiver = vertex.receiver,
                         layout = "hierarchy")
    dev.off()
  }
  gg <- netAnalysis_contribution(ctr.cc, signaling = pathways.show.all[i])
  ggsave(filename=paste0(fname, "/Results/LR/", pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 7, height = 7, units = 'in', dpi = 300)
}

for(i in 1:length(cell.ident)){
  all_c = c(1:6)
  pdf(paste0(fname, "/", cell.ident[i], "to_all.pdf"))
  j = all_c[-i]
  p <- netVisual_bubble(ctr.cc, sources.use = 1, targets.use = j, remove.isolate = FALSE)
  print(p)
  dev.off()
}

ctr.cc <- netAnalysis_computeCentrality(ctr.cc, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

dir.create(paste0(fname, "/CentralityPlots"))
for(i in 1:length(pathways.show.all)){
  pdf(paste0(fname, "/CentralityPlots/", pathways.show.all[i], ".pdf"))
  netAnalysis_signalingRole_network(ctr.cc, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10)
  dev.off()
}

dir.create(paste0(fname, "/Incoming_Outgoing_Strength"))
for(i in 1:length(pathways.show.all)){
  pdf(paste0(fname, "/Incoming_Outgoing_Strength/", pathways.show.all[i], ".pdf"))
  gg2 <- netAnalysis_signalingRole_scatter(ctr.cc, signaling = pathways.show.all[i], title = pathways.show.all[i])
  print(gg2)
  dev.off()
}
pdf(paste0(fname, "/Signalling_pattern.pdf"), width = 10, height = 10)
ht1 <- netAnalysis_signalingRole_heatmap(ctr.cc, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(ctr.cc, pattern = "incoming")
ht1 + ht2
dev.off()


library(NMF)
library(ggalluvial)
selectK(ctr.cc, pattern = "outgoing")
nPatterns = 3
pdf(paste0(fname, "/outgoing_patterns.pdf"), width = 10, height = 5)
ctr.cc <- identifyCommunicationPatterns(ctr.cc, pattern = "outgoing", k = nPatterns)
dev.off()

pdf(paste0(fname, "/Outgoing_pattern_river.pdf"))
netAnalysis_river(ctr.cc, pattern = "outgoing")
dev.off()

pdf(paste0(fname, "/outgoing_pattern_dot.pdf"))
netAnalysis_dot(ctr.cc, pattern = "outgoing")
dev.off()

##Incomingnow
selectK(ctr.cc, pattern = "incoming")
nPatterns = 4
pdf(paste0(fname, "/incoming_patterns.pdf", width = 10))
ctr.cc <- identifyCommunicationPatterns(ctr.cc, pattern = "incoming", k = nPatterns)
dev.off()

# river plot
pdf(paste0(fname, "/incoming_river.pdf"))
netAnalysis_river(ctr.cc, pattern = "incoming")
dev.off()
pdf(paste0(fname, '/incoming_dot.pdf'))
netAnalysis_dot(ctr.cc, pattern = "incoming")
dev.off()


ctr.cc <- computeNetSimilarity(ctr.cc, type = "functional")

library(reticulate)
use_python("/Users/sandeshacharya/Documents/CellChat/miniconda/envs/r-reticulate/bin/python")
ctr.cc <- netEmbedding(ctr.cc, type = "functional")
ctr.cc <- netClustering(ctr.cc, type = "functional", do.parallel = FALSE)
dir.create(paste0(fname, "/UMAP"))

pdf(paste0(fname, "/UMAP/Umap_Functional.pdf"))
netVisual_embedding(ctr.cc, type = "functional", label.size = 3.5)
dev.off()
pdf(paste0(fname, "/UMAP/Umap_Functional_sep.pdf"))
netVisual_embeddingZoomIn(ctr.cc, type = "functional", nCol = 2)
dev.off()


ctr.cc <- computeNetSimilarity(ctr.cc, type = "structural")
ctr.cc <- netEmbedding(ctr.cc, type = "structural")
ctr.cc <- netClustering(ctr.cc, type = "structural", do.parallel = FALSE)
pdf(paste0(fname, "/UMAP/umap_structural.pdf"))
netVisual_embedding(ctr.cc, type = "structural", label.size = 3.5)
dev.off()
pdf(paste0(fname, "/UMAP/umap_structural_sep.pdf"))
netVisual_embeddingZoomIn(ctr.cc, type = "structural", nCol = 2)
dev.off()

saveRDS(ctr.cc, file = "ctr.cc_all_all.rds")
