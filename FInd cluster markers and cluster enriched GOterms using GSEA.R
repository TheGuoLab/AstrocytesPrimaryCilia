###################Find cluster marker genes with no threshold##################
scRNA<- readRDS('Astrocytes_Seurat.rds')
scRNA.markers <- FindAllMarkers(scRNA, only.pos = FALSE, min.pct = 0, logfc.threshold = 0)

setwd('xxxx')
saveRDS(scRNA.markers,'AST_all_cluster.maerker_no.logfc.threshold.rds')


###################Find cluster enriched GO terms##################
#library packages
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(stats)
msigdbr_species()
msigdbr_collections()


Idents(scRNA)<- 'seurat_clusters'


#GO:BP GSEA
genesets = msigdbr(species = "Mus musculus", category = "C5") 
#View(msigdbr_collections()) 
unique(genesets$gs_subcat)  
genesets <- subset(genesets, gs_subcat=="GO:BP", select = c("gs_name", "gene_symbol"))
length(unique(genesets$gs_name)) #7656 terms

#GSEA analysis
for(i in c(0:4)){
  x<- subset(scRNA.markers, cluster==i)
  x = x[order(x$avg_log2FC, decreasing = T),]
  row.names(x)<- x$gene
  genelist <- structure(x$avg_log2FC, names = rownames(x))
  res <- GSEA(genelist, TERM2GENE = genesets, eps = 0,minGSSize = 0,maxGSSize = 1000, pvalueCutoff = 1)
  saveRDS(res, file = paste0('Cluster_',i,'_enriched_GOBP_GSEA_res.rds'))
  res<- as.data.frame(res)
  write.csv(res, file = paste0('Cluster_',i,'_enriched_GOBP_GSEA_res.csv'))
}