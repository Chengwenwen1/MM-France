setwd('*/')
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
filename <- list.files(getwd())#,pattern = '_out$')
list <-list()
for (i in 1:length(filename)){
  list[[i]] <- Read10X(filename[i])#paste('./',filename[i],'/outs/filtered_feature_bc_matrix/',sep = '')
  list[[i]] <- CreateSeuratObject(list[[i]],project = filename[i],min.cells = 5,min.features = 50)
  list[[i]] <- RenameCells(list[[i]],add.cell.id=filename[i])
}
list
normal.bm <- merge(list[[1]],list[2:length(list)])
table(normal.bm$orig.ident)
normal.bm$Patient <- normal.bm$orig.ident
normal.bm$orig.ident <-'normal BM'

normal.bm[['percent.mt']] <- PercentageFeatureSet(normal.bm,pattern = '^MT-')
#normal.bm[['percent.rp']] <- PercentageFeatureSet(normal.bm,pattern = '^RP[SL]')
VlnPlot(normal.bm,features = c('nCount_RNA','nFeature_RNA','percent.mt'),pt.size = 0,group.by = 'Patient',ncol = 3,
        cols = resp.col,split.by = 'Patient')&labs(x='')

normal.bm <- subset(normal.bm,percent.mt < 8 & 500 < nFeature_RNA & nCount_RNA < 50000 ) 
table(normal.bm$orig.ident)
#Normalize #scale ##PCA
normal.bm <- NormalizeData(normal.bm,normalization.method = "LogNormalize")
normal.bm <- FindVariableFeatures(normal.bm,selection.method = 'vst',nfeatures = 2000, verbose = FALSE)
normal.bm <- ScaleData(normal.bm,features = VariableFeatures(normal.bm))
normal.bm <- RunPCA(normal.bm,features = VariableFeatures(normal.bm))
LabelPoints(VariableFeaturePlot(normal.bm),points = head(VariableFeatures(normal.bm),10),repel = T)
DimPlot(normal.bm,reduction = 'pca',group.by = 'Patient')
ElbowPlot(normal.bm,ndims = 50)
normal.bm <- FindNeighbors(normal.bm,dims = 1:20)
normal.bm <- FindClusters(normal.bm,resolution = 1)
table(normal.bm$seurat_clusters)
normal.bm <- RunUMAP(normal.bm,dims = 1:20) ##umap/tsne
############batch effect
cart.col <- colorRampPalette(brewer.pal(n=3,name = 'Dark2'))(7)
patie.col <- colorRampPalette(brewer.pal(n=9,name = 'Set1'))(46)
resp.col <- c("#35978F","#BF812D","#85658D")
DimPlot(normal.bm, reduction = 'umap',label = T,cols = patie.col,group.by = 'Patient')
mycol_all1 <- colorRampPalette(c(brewer.pal(n=12,name='Paired'),brewer.pal(n=8,name='Accent')))(80)
p1 <- DimPlot(normal.bm, reduction = 'umap',label = T,cols = mycol_all,group.by = 'seurat_clusters')+NoLegend()

Idents(normal.bm) <-'seurat_clusters'
p2<- DotPlot(normal.bm,features = big.marker %>% unique(),cluster.idents = T,group.by = 'seurat_clusters')+
  theme(axis.text.x = element_text(size=10,angle = 60,vjust = 1, hjust=1)) + NoLegend() 
cowplot::plot_grid(p1,p2,rel_widths = c(0.7,1))
##
normal.bm <- subset(normal.bm,seurat_clusters %in% setdiff(0:30,27:29))
saveRDS(normal.bm,file = '*/normal_BM.rds')
##------------------------------------------------------

##Step---integration for normal & MM cells--------------------------------------
#read normal & MM data
normal.bm <- readRDS('*/normal_BM.rds')
MM_obj <- readRDS('*/MM_obj_celltype.RDS')
#
MM_obj$orig.ident; normal.bm$orig.ident
mm_list <- list(MM_obj,normal.bm)
mm_list <- lapply(mm_list, function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = mm_list)
mm_mye <- FindIntegrationAnchors(object.list = mm_list, anchor.features = features,
                                 dims = 1:30,k.anchor = 5,k.filter = 20)
mm_normal <- IntegrateData(anchorset = mm_mye, dims = 1:30)
DefaultAssay(mm_normal) <- "integrated"
mm_normal <- ScaleData(mm_normal,features = VariableFeatures(mm_normal))
mm_normal <- RunPCA(mm_normal,features = VariableFeatures(mm_normal))
LabelPoints(VariableFeaturePlot(mm_normal),points = head(VariableFeatures(mm_normal),10),repel = T)
DimPlot(mm_normal,reduction = 'pca',group.by = 'orig.ident')
ElbowPlot(mm_normal,ndims = 50)
mm_normal <- FindNeighbors(mm_normal,dims = 1:20)
mm_normal <- FindClusters(mm_normal,resolution = 2)
table(mm_normal$seurat_clusters)
mm_normal <- RunUMAP(mm_normal,dims = 1:20) ##umap/tsne
DimPlot(mm_normal, reduction = 'umap',label = T,cols = mycol_all,split.by  = 'orig.ident')
DimPlot(mm_normal, reduction = 'umap',label = T,cols = resp.col,group.by  = 'orig.ident')
##
DefaultAssay(mm_normal) <- "RNA"
mm_normal <- NormalizeData(mm_normal,normalization.method = "LogNormalize")
mm_normal <- FindVariableFeatures(mm_normal,selection.method = 'vst',nfeatures = 2000, verbose = FALSE)
mm_normal <- ScaleData(mm_normal,features = VariableFeatures(mm_normal))
p1 <- DimPlot(mm_normal, reduction = 'umap',label = T,cols = mycol_all1,group.by = 'seurat_clusters')+NoLegend()
Idents(mm_normal)<-'seurat_clusters'
p2<- DotPlot(mm_normal,features = big.marker %>% unique(),cluster.idents = T,group.by = 'seurat_clusters')+
  theme(axis.text.x = element_text(size=10,angle = 60,vjust = 1, hjust=1)) + NoLegend() 
p1|p2
FeaturePlot(mm_normal,features = c('LYZ','GZMB'),split.by  = 'orig.ident')
#saveRDS(mm_normal,file = '*/mm_normal1.rds')

#mm_normal <- subset(mm_normal,seurat_clusters %in% setdiff(0:52,c(45,48)))

#mm_normal <- subset(mm_normal,seurat_clusters %in% setdiff(0:51,c(38)))
saveRDS(mm_normal,file = '*/mm_normal_2.rds')


DefaultAssay(mm_normal) <- "integrated"
mm_normal <- ScaleData(mm_normal,features = VariableFeatures(mm_normal))
mm_normal <- RunPCA(mm_normal,features = VariableFeatures(mm_normal))
LabelPoints(VariableFeaturePlot(mm_normal),points = head(VariableFeatures(mm_normal),10),repel = T)
DimPlot(mm_normal,reduction = 'pca',group.by = 'orig.ident')
ElbowPlot(mm_normal,ndims = 50)
mm_normal <- FindNeighbors(mm_normal,dims = 1:20)
mm_normal <- FindClusters(mm_normal,resolution = 3)
table(mm_normal$seurat_clusters)
mm_normal <- RunUMAP(mm_normal,dims = 1:20) ##umap/tsne
DimPlot(mm_normal, reduction = 'umap',label = T,split.by  = 'orig.ident')
DimPlot(mm_normal, reduction = 'umap',label = T,cols = resp.col,group.by  = 'orig.ident',
        raster = F)

DefaultAssay(mm_normal) <- "RNA"
mm_normal <- NormalizeData(mm_normal,normalization.method = "LogNormalize")
mm_normal <- FindVariableFeatures(mm_normal,selection.method = 'vst',nfeatures = 2000, verbose = FALSE)
mm_normal <- ScaleData(mm_normal,features = VariableFeatures(mm_normal))
p1 <- DimPlot(mm_normal, reduction = 'umap',label = T,group.by = 'seurat_clusters')+NoLegend()
Idents(mm_normal)<-'seurat_clusters'
p2<- DotPlot(mm_normal,features = big.marker %>% unique(),cluster.idents = T,group.by = 'seurat_clusters')+
  theme(axis.text.x = element_text(size=10,angle = 60,vjust = 1, hjust=1)) + NoLegend() 
p1|p2
FeaturePlot(mm_normal,features = c('CLEC9A','CPNE3','WDFY4','XCR1',"CLEC10A",'FCER1A','CD1C'))
###############cluster annotation
n=length(unique(mm_normal$seurat_clusters))
celltype=data.frame(ClusterID=0:(n-1),celltype=0:(n-1),stringsAsFactors = F)
celltype[celltype$ClusterID %in% c(53,32,55,37,24),2]='HSPC'
celltype[celltype$ClusterID %in% c(49,29,50,22,62),2]='EPC'
celltype[celltype$ClusterID %in% c(61),2]='Stromal cells'
celltype[celltype$ClusterID %in% c(5,1,3,27,18,59,0,6,2),2]='CD4_T'
celltype[celltype$ClusterID %in% c(16,9,4,8,45,11),2]='CD8_T'
celltype[celltype$ClusterID %in% c(13),2]='NK'
celltype[celltype$ClusterID %in% c(39),2]='Pro B'
celltype[celltype$ClusterID %in% c(41),2]='Pre B'
celltype[celltype$ClusterID %in% c(7),2]='Im B'
celltype[celltype$ClusterID %in% c(23),2]='Mature B'
celltype[celltype$ClusterID %in% c(17,15,19,12,10,28,46,63),2]='cMo'
celltype[celltype$ClusterID %in% c(21),2]='CD16+ Mo'
celltype[celltype$ClusterID %in% c(31),2]='GMP'
celltype[celltype$ClusterID %in% c(30,58),2]='preMono'
celltype[celltype$ClusterID %in% c(34),2]='prolMono'
celltype[celltype$ClusterID %in% c(38,43),2]='DC'
celltype[celltype$ClusterID %in% c(52,40,25,42,56,26,47,48,60,54,14,36,33,44,57,51,20,35),2]='MMPC'
mm_normal@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  mm_normal@meta.data[which(mm_normal@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(mm_normal$celltype)
mm_normal$celltype <- ifelse(mm_normal$orig.ident=='normal BM' & mm_normal$celltype=='MMPC',
                             'nPC',mm_normal$celltype) 

levels_ctype_all=c('HSPC','EPC','Stromal cells', 'CD4_T','CD8_T','NK',
                   'Pro B','Pre B','Im B','Mature B','cMo','CD16+ Mo','GMP','preMono',
                   'prolMono','DC','nPC','MMPC')
mm_normal$celltype <-factor(mm_normal$celltype, levels=levels_ctype_all)

mycol_all <- colorRampPalette(c(brewer.pal(n=12,name='Paired'),
                                brewer.pal(n=5,name='Accent'),
                                brewer.pal(n=5,name = 'Set1')))(36)
mycol_all <- c('HSPC'="#A6CEE3",'EPC'="#559AC6",
               'Stromal cells'="#3C8CAB",'CD4_T'="#94CA92",
               'CD8_T'="#33A02C",'NK'="#AB9C6D",
               'Pro B'="#F68080",'Pre B'="#E73334",'Im B'="#ED5B3D",'Mature B'="#FDBF6F",
               'cMo'="#FE982C",'CD16+ Mo'="#D4A7AB",'GMP'="#A383BE",
               'preMono'="#C3B199",'prolMono'="#EFDD82",'DC'="#C07A3E",#'pDC'="#9C854A",
               'nPC'="#bc8f8f",'MMPC'="#A4B8B1")

p1 <- DimPlot(mm_normal,reduction = 'umap',label = T,cols =mycol_all,
              group.by ='celltype',repel = T,raster = T)
p4 <- DotPlot(mm_normal, features = big.marker %>% unique(),
              group.by = 'celltype') +
  NoLegend()+RotatedAxis()+
  xlab('')+ylab('')+
  scale_color_gradient(low = "lightgrey", high = "#C13533")+
  gg.theme+
  theme(axis.text.x=element_text(angle=90,hjust=1))
cowplot::plot_grid(p1,p4,rel_widths = c(0.7,1))


##change normal samples lables
table(mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396161','A',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396162','B',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient %in% c('GSM3396163','GSM3396164','GSM3396165'),'C',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396166','E',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396167','F',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396168','G',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396169','H',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396170','J',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396171','K',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396172','L',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396173','M',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396174','N',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396175','O',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396176','P',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396177','Q',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396178','R',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient %in% c('GSM3396179','GSM3396180','GSM3396181','GSM3396182'),'S',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396183','T',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396184','U',mm_normal$Patient)
mm_normal$Patient <- ifelse(mm_normal$Patient=='GSM3396185','W',mm_normal$Patient)


##Step--------------------------------------
mm_normal$celltype_big <- as.character(mm_normal$celltype)
mm_normal$celltype_big <- ifelse(mm_normal$celltype %in% c('CD4_T','CD8_T','NK'),
                                 'T & NK',mm_normal$celltype_big)
mm_normal$celltype_big <- ifelse(mm_normal$celltype %in% 
                                   c('cMo','CD16+ Mo','GMP','preMono','prolMono','DC'),
                                 'Myeloid cells',mm_normal$celltype_big)  
mm_normal$celltype_big <- ifelse(mm_normal$celltype %in% c('Pro B','Pre B','Im B','Mature B'),
                                 'B cells',mm_normal$celltype_big)                                                       
table(mm_normal$celltype_big )
celltype_big_levles <- c( 'HSPC','EPC','Stromal cells','T & NK','B cells','Myeloid cells', 'nPC','MMPC')
mm_normal$celltype_big <-factor(mm_normal$celltype_big, levels=celltype_big_levles)
bigcol <- c('HSPC'="#A6CEE3",'EPC'="#559AC6",'Stromal cells'="#33A02C",'T & NK'="#9D9E6C",
            'B cells'="#D199B9",'Myeloid cells'="#F7B56B",'nPC'="#bc8f8f",'MMPC'="#A4B8B1"
)
DimPlot(mm_normal,reduction = 'umap',label = T,cols =bigcol,
        group.by ='celltype_big',repel = T,raster = T)+labs(title = 'celltype')
saveRDS(mm_normal,file = '*/mm_normal_common.rds')


