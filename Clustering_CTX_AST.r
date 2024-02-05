rm(list=ls())
######################Library packages###########################
library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)
library(sctransform)
library(harmony)

######################Create Seurat object########################
setwd('D:/Lizheng/astrocytes_GUO_lab_scRNA_seq/Analysis/2023Analysis for NN revision/For GEO')
all_cell <- readRDS('all_cell_matrix.rds')
meta <- read.csv('all_meta.csv', row.names=1, header=T)

all_cell <- CreateSeuratObject(all_cell,project = "all_cell", min.cells = 3, min.features = 200, meta.data = meta)

##########################QC####################################
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
all_cell[["percent.mt"]] <- PercentageFeatureSet(all_cell, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(all_cell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"),pt.size = 0, ncol = 3)

plot1 <- FeatureScatter(all_cell, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_cell, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

dim(all_cell) 
all_cell <- subset(all_cell, subset = nFeature_RNA >200 &
                     nFeature_RNA < 6000 &
                     percent.mt <25)
dim(all_cell)

VlnPlot(all_cell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol = 4))
rm(list=c('plot1', 'plot2'))


###########Normalization and runPCA########################
all_cell <- SCTransform(all_cell, vars.to.regress = 'percent.mt', verbose = FALSE)

all_cell<-RunPCA(all_cell,assay="SCT",verbose = FALSE)

###########remove batch effect using Harmony#################
table(all_cell$sample_unique_name)
all_cell=RunHarmony(all_cell,group.by.vars="sample_unique_name",assay.use="SCT", plot_convergence = TRUE,max.iter.harmony =50 )
#######Determine the best PCs for use
pct <- all_cell [["harmony"]]@stdev / sum( all_cell [["harmony"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs


#Use the best pcs for UMAP and FindNeighbors
bestpc=1:pcs
all_cell<- all_cell %>% RunUMAP(reduction = "harmony", dims = bestpc) %>% 
  FindNeighbors(reduction = "harmony", dims = bestpc)

all_cell<-FindClusters(all_cell,resolution = 0.2, future.seed=TRUE)

DimPlot(all_cell, reduction = "umap", label = TRUE,repel = T, pt.size = .3)

saveRDS(all_cell,'Astrocytes_Seurat.rds')
