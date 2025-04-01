setwd('*/')
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
MM_t <- readRDS('./MM_t_celltype_inte.RDS')
mycol_t <- c('CD4 TN'='#A6CEE3','CD4 TCM'='#3385BB','CD4 TEM'='#F9AA64',
             'Th17'='#b4b3d7','Treg'='#6DBD57','CD8 TN'='#a98588','CD8 TM'='#D8CD99',
             'CD8 tox'='#D4A7AA','IFN T'="#9f2920",
             'gdT'="#64A4CC",'NKT'="#569993",'NKdim'="#6968af",'NKbright'='#b3c970')
DefaultAssay(MM_t) <-'RNA'
Idents(MM_t)<-'celltype'
DimPlot(MM_t,reduction = 'umap',label = F,cols =mycol_t,repel = T,raster = T)
table(MM_t$celltype)

MM_t.sce <- as.SingleCellExperiment(MM_t)

reducedDims(MM_t.sce) <- SimpleList(UMAP=Embeddings(MM_t,reduction='umap') )
print(MM_t.sce)

library(slingshot)
library(SingleCellExperiment)
sim <- slingshot(MM_t.sce,clusterLabels = 'celltype',reducedDim ='UMAP')

colData(sim)$celltype <- as.character(colData(sim)$celltype)
plot(reducedDims(sim)$UMAP,col =  mycol_t[colData(sim)$celltype], pch=16, asp = 1,
     cex=0.6,
     axes=F,ann=F)
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
title(main = 'Slingshot')



###
mm_myeloid <- readRDS('*/mm_myeloid_celltype.RDS')
mycol <- c('cMono (CD14+)'="#A6CEE3",'cMono (CD14+IFI6+)'="#3385BB",
           'cMono (CD14+THBS1+)'="#B2DF8A",'cMono (CD14+THBD+)'="#33A02C",
           'Intermediate Mono'='#ddbbb9','Nonclassical Mono'="#d57878",
           'Macrophages'="#e34a1a",
           'prolMono'= "#d46710",'preMono'= "#bc8756",'GMP'="#bc6556",
           'cDC1'="#D8CDE4",'cDC2'="#966bc4",'cDC5'='grey','pDC'="#858500" )
DefaultAssay(mm_myeloid) <-'RNA'
Idents(mm_myeloid)<-'celltype'
DimPlot(mm_myeloid,reduction = 'umap',label = F,cols =mycol,repel = T,raster = T)
table(mm_myeloid$celltype)
##sce
mm_myeloid.sce <- as.SingleCellExperiment(mm_myeloid)

reducedDims(mm_myeloid.sce) <- SimpleList(UMAP=Embeddings(mm_myeloid,reduction='umap') )
print(mm_myeloid.sce)
###singleshot
library(slingshot)
library(SingleCellExperiment)
sim <- slingshot(mm_myeloid.sce,clusterLabels = 'celltype',reducedDim ='UMAP')

colData(sim)$celltype <- as.character(colData(sim)$celltype)
plot(reducedDims(sim)$UMAP,col =  mycol[colData(sim)$celltype], pch=16, asp = 1,
     cex=0.6,#点的大小
     axes=F,ann=F)
lines(SlingshotDataSet(sim), lwd=2, type = 'lineages', col = 'black', show.constraints = TRUE)
title(main = 'Slingshot')
