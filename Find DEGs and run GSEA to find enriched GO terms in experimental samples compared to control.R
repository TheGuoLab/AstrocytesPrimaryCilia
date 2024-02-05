##################Find DEGs without threshold compared to Control#################################
rm(list=ls())
setwd('xxx')
scRNA<- readRDS('Astrocytes_Seurat.rds')

Idents(scRNA)<- 'group.name'
table(scRNA$group.name)
clusters<- c(0:4)
samples<- c('Arl13b_cKO', 'Ift88_cKO','Smo_cKO','SmoM2-Arl13b_cKO')
Control<- 'Control'

setwd('yyy') #path to save the DEG files
for(i in 1:length(clusters)){
  for(j in 1:length(samples)){
    temp<- FindMarkers(subset(scRNA,seurat_clusters==clusters[i]), ident.1= samples[j], ident.2 = Control, 
                       logfc.threshold = 0, verbose =T , min.pct)
    temp$comparison <- paste0(samples[j],'_vs_Control_', clusters[i])
    write.csv(temp,file = paste0(samples[j],'_vs_Control_', clusters[i],'_DEGs.csv'))
  } 
}


##################Run GSEA using DEG files#################################

library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(stats)


msigdbr_species()
msigdbr_collections()

genesets = msigdbr(species = "Mus musculus", category = "C5") 
View(msigdbr_collections()) 
unique(genesets$gs_subcat)  
genesets <- subset(genesets, gs_subcat=="GO:BP", select = c("gs_name", "gene_symbol"))

length(unique(genesets$gs_name))

setwd('yyy') # path for the DEG files
files<- dir(pattern = 'DEGs.csv')
files

for(a in 1: length(files)){
  dt<- read.csv(files[a], header = T, row.names = 1)
  x<- dt
  x = x[order(x$avg_log2FC, decreasing = T),]
  genelist <- structure(x$avg_log2FC, names = rownames(x))
  genelist
  res <- GSEA(genelist, TERM2GENE = genesets, eps = 0,pvalueCutoff = 1)
  setwd('path to save GSEA results') #path to save GSEA results
  saveRDS(res, file = paste0('DEGs_',files[a],'_','_GOBP_GSEA_res.rds'))
  res<- as.data.frame(res)
  write.csv(res, file = paste0('DEGs_',files[a],'_','_GOBP_GSEA_res.csv'))
  setwd('yyy') # path for the DEG files
}


