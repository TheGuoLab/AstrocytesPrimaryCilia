setwd("~/Documents/CB_data_for_SCENIC_to_Sandesh/scenic_counts/")
library(Seurat)

#Processed data (for Cortex or Cerebellum) contains a seurat object with cluster information. 
combined <- readRDS("../processed_data.rds")
# name of Seurat object used
cell.type.set <- unique(combined$seurat_clusters) #cell type labels

# Identify genes with initial evidence of being differentially expressed between the two conditions per cell type
rna.marker.list <- list()
mutation = "3.SA"
mut_name <- "SA"
  
diffexp = c()

#Note: the DEGs for each mutation compared with the control for each cluster is already identified and saved in separate files. 
for(i in 0:3){
  filename = paste0("DEGs/CB_AST_", mutation, "_vs_CTR_CLuster..", i, ".._DEGs.csv")
  print(c("Currently Reading file: ", str_replace(filename, "DEGs/", "")))
  data <- fread(filename)
  data <- data[data$p_val_adj <= 0.1, ]
  rna.marker.list[[i+1]] = data
  }
# get entire RNA perturb/differential genes
rna.perturb.set <- sapply(rna.marker.list, function(x)  x$V1 )
rna.perturb.set <- unlist(rna.perturb.set)
names(rna.perturb.set) <- NULL
rna.perturb.set <- unique(rna.perturb.set)

rna.marker.num <- sapply(rna.marker.list, function(x) nrow(x))


Idents(combined) <- "seurat_clusters"
rna.m <- (combined[["SCT"]]@data[rna.perturb.set,])
rownames(rna.m ) <- paste0("SCT-", rownames(rna.m )  )
combined[["both"]] <- CreateAssayObject(data = rna.m )

pertb_scores = c()

for (i in 0:3){
  i.cell <- WhichCells(combined , expression = seurat_clusters == i)
  
  if(( nrow(rna.marker.list[[i+1]])  ) >1 ){
    combined.i  <- subset(combined,  cells = i.cell )
    
    DefaultAssay(combined.i) <- "both"
    combined.i <- DietSeurat(combined.i, assays = "both")
    combined.i.0 <- subset(combined.i, subset = group.name ==mutation)
    combined.i.2 <- subset(combined.i, subset = group.name =="1.CTR")
    message(paste0("Finished subsetting ",i))

    global.perturb.rna <- NULL
    if( length(rna.marker.list[[ i+1  ]]$V1) > 0 ){
      perturb.feature.rna <- paste0( "SCT-",rna.marker.list[[ i+1 ]]$V1)
      global.perturb.rna <-  rna.marker.list[[ i+1 ]]$avg_log2FC
      global.perturb.rna[global.perturb.rna==Inf] <- max(global.perturb.rna[!global.perturb.rna==Inf])
      names(  global.perturb.rna ) <- perturb.feature.rna ### !!!
    }
    global.perturb.vec <- c(global.perturb.rna)
    global.perturb.vec <- abs(as.matrix(global.perturb.vec))
    
    # normalize its length
    global.perturb.vec <- global.perturb.vec/ sqrt(sum(global.perturb.vec**2))
    perturb.feature <- rownames( global.perturb.vec )
    
    subcells <- combined.i[["both"]]@data[perturb.feature.rna, ]
    scoress <- t(subcells) %*% global.perturb.vec
    scoress <- scoress/ sqrt(sum(scoress**2))
    scoress <- data.frame(scoress)
    scoress$samples <- colnames(subcells)
    pertb_scores <- rbindlist(list(pertb_scores, data.frame(scoress)))
    
  } else {
    scoress <- data.frame(scoress = 0, samples = i.cell )
    pertb_scores <- rbindlist(list(pertb_scores, data.frame(scoress)))
    
  }
}

pertb_scores <- pertb_scores[match(rownames(combined@meta.data), pertb_scores$samples), ]
sum(pertb_scores$samples == rownames(combined@meta.data))


eval(parse(text = paste(paste0("combined$", mut_name), "pertb_scores$scoress", sep = "<-")))
dir.create("DensityPlots")
svg(paste0("DensityPlots/",mut_name, "Density.svg"))
p <- plot_density(combined, features = mut_name) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.7,0.1),
        legend.direction="horizontal",
        line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

print(p)
  dev.off()

