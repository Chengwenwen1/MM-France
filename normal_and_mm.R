setwd('/data/jinwf/chengww/project/multiple_myeloma/raw_data/GSE120221_RAW_scRNA_20_health_BM/raw_counts/')
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
############查看umap坐标下的batch effect
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
##去除mix细胞群#后再次运行聚类
normal.bm <- subset(normal.bm,seurat_clusters %in% setdiff(0:30,27:29))
saveRDS(normal.bm,file = '/home/chengww/data/project/multiple_myeloma/midel_result_R4/normal_BM.rds')
##已完成处理normal BM的过滤-------------------------------------------------------

##Step---integration for normal & MM cells--------------------------------------
#read normal & MM data
normal.bm <- readRDS('/home/chengww/data/project/multiple_myeloma/midel_result_R4/normal_BM.rds')
MM_obj <- readRDS('/home/chengww/data/project/multiple_myeloma/midel_result_R4/MM_obj_celltype.RDS')
#整合分析
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
##基于原始的RNA矩阵查看基因表达情况
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
#saveRDS(mm_normal,file = '/home/chengww/data/project/multiple_myeloma/midel_result_R4/mm_normal1.rds')
##第一次清洗
#mm_normal <- subset(mm_normal,seurat_clusters %in% setdiff(0:52,c(45,48)))
##第二次清洗#高表达cd4和GZMB 不能定义
#mm_normal <- subset(mm_normal,seurat_clusters %in% setdiff(0:51,c(38)))
saveRDS(mm_normal,file = '/home/chengww/data/project/multiple_myeloma/midel_result_R4/mm_normal_2.rds')


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
##基于原始的RNA矩阵查看基因表达情况
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


##Step--重新生成一列为了，按照大群定义亚群--------------------------------------
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
saveRDS(mm_normal,file = '/home/chengww/data/project/multiple_myeloma/midel_result_R4/normal_mm/mm_normal_common.rds')


##Step--为了分析细胞总数的一致性，再次清理大细胞群按照以下方式------------------
#1.先对各个亚群进行更细亚群的鉴定，基于mm_normal_common.rds
#2.对所有细亚群划分后删除了一些混杂群，将混杂群返回在mm_normal_common_all.rds中删除
#3.保存的final data是剔除干净的数据 ##最终用来大亚群结果的展示，as Fig1
mm_normal <- readRDS('/home/chengww/data/project/multiple_myeloma/midel_result_R4/normal_mm/mm_normal_common.rds')
mm_myeloid <- readRDS('./mm_myeloid_celltype.RDS')
MM_t <- readRDS('./MM_t_celltype_inte.RDS')
mm_pc <- readRDS('./mm_pc_add_normal.RDS')

##提取所有细胞的ID
all.cells <- c(colnames(mm_myeloid),colnames(MM_t),colnames(mm_pc),
               colnames(subset(mm_normal,celltype_big %in% c('HSPC','EPC','Stromal cells','B cells'))))
mm_normal$cellid.m <- rownames(mm_normal[[]])
mm_normal_all <- subset(mm_normal,cellid.m %in% all.cells)

DimPlot(mm_normal_all,reduction = 'umap',label = T,cols =bigcol,
        group.by ='celltype_big',repel = T,raster = T)+
  labs(title = 'celltype')+gg.theme+NoLegend()


##修改metadata，变更临床信息（性别 年龄 临床癌症阶段）----------------------
##change sample tags
mm_normal_all@meta.data$Patient <- ifelse(mm_normal_all$orig.ident=='normal BM',mm_normal_all@meta.data$Patient,
                                          paste('A',gsub('-','',mm_normal_all@meta.data$Patient),sep = ''))
mm_normal_all@meta.data$cellid.m <- rownames(mm_normal_all@meta.data)#给列名
##read new meta data
meta <- fread('/data/jinwf/chengww/project/multiple_myeloma/raw_data/mm_metadata.csv') %>% 
  select(Patient,Response,Age,Gender,SCISS)
mm_normal_all@meta.data <- mm_normal_all@meta.data %>% select(-c(Response,Age,Gender))
mm_normal_all@meta.data <- left_join(mm_normal_all@meta.data,meta,by='Patient')
rownames(mm_normal_all@meta.data)<- mm_normal_all@meta.data$cellid.m#完成
##修改noraml BM 为HD
mm_normal_all$Response <- as.character(mm_normal_all$Response)
mm_normal_all$Response <- ifelse(mm_normal_all$Response=='normal BM','HD',mm_normal_all$Response)
mm_normal_all$Response <- factor(mm_normal_all$Response,levels = c('CR','NR','HD'))
mm_normal_all$Patient <- factor(mm_normal_all$Patient,levels = unique(mm_normal_all$Patient))


##对新筛选过后的细胞矩阵做标准化处理
DefaultAssay(mm_normal_all) <- "integrated"
mm_normal_all <- ScaleData(mm_normal_all,features = VariableFeatures(mm_normal_all))
mm_normal_all <- RunPCA(mm_normal_all,features = VariableFeatures(mm_normal_all))
DimPlot(mm_normal_all,reduction = 'pca',group.by = 'orig.ident')
ElbowPlot(mm_normal_all,ndims = 50)
mm_normal_all <- FindNeighbors(mm_normal_all,dims = 1:16)
mm_normal_all <- FindClusters(mm_normal_all,resolution = 2)
mm_normal_all <- RunUMAP(mm_normal_all,dims = 1:16) ##umap/tsne
DimPlot(mm_normal_all,reduction = 'umap',label = T,repel = T,raster = T,cols =mycol_all,
        group.by = 'celltype')
DimPlot(mm_normal_all,reduction = 'umap',label = T,repel = T,raster = F,cols =bigcol,
        group.by = 'celltype_big')+labs(title = 'Celltype (N = 149365)')+gg.theme

##基于原始的RNA矩阵查看基因表达情况
DefaultAssay(mm_normal_all) <- "RNA"
mm_normal_all <- NormalizeData(mm_normal_all,normalization.method = "LogNormalize")
mm_normal_all <- FindVariableFeatures(mm_normal_all,selection.method = 'vst',nfeatures = 2000, verbose = FALSE)
mm_normal_all <- ScaleData(mm_normal_all,features = VariableFeatures(mm_normal_all))
DimPlot(mm_normal_all,reduction = 'umap',label = T,cols =mycol_all,
        group.by ='celltype',repel = T,raster = T)

DotPlot(mm_normal_all, features = all.marker %>% unique(),
        group.by = 'celltype_big') +
  NoLegend()+RotatedAxis()+
  xlab('')+ylab('')+
  scale_color_gradient(low = "lightgrey", high = "#C13533")+
  gg.theme+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))

##查看抗体的表达情况
DefaultAssay(mm_normal_all) <- "ADT"
mm_all <-subset(mm_normal_all,Response !='HD')
mm_all <- NormalizeData(mm_all, normalization.method = "CLR", margin = 2)
adt.markers <-  c(      
  "CD34:8G12-CD34-AHS0182-pAbO", #HSPC
  "CD45-PTPRC-AHS0040-pAbO",           
  "CD4:SK3-CD4-AHS0032-pAbO","CD8:RPA-T8-CD8A-AHS0027-pAbO",
  "CD161:HP-3G10-KLRB1-AHS0205-pAbO",#T &nk
  "CD19:HIB19-CD19-AHS0161-pAbO",#b
  "CD11b:M1-70-ITGAM-AHS0005-pAbO",#髓系
  "CD38:HIT2-CD38-AHS0022-pAbO",#pc
  'CD183-CXCR3-AHS0031-pAbO'
)
DotPlot(mm_all, features = adt.markers %>% unique(),
        group.by = 'celltype_big') +
  NoLegend()+RotatedAxis()+
  xlab('')+ylab('')+
  scale_color_gradient(low = "lightgrey", high = "#C13533")+
  gg.theme+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))


FeaturePlot(mm_all,features = adt.markers[-8],#order = T,
            raster = T,
            cols =c("lightgrey",'#BE766E','#C13533'),ncol = 4)&
  theme_bw() &
  theme(#legend.position = c(0.1,0.88),
    panel.border =element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    text = element_text(size=12,face="bold"),
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_blank(), 
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.line = element_blank())&NoLegend()
##save final data obj
saveRDS(mm_normal_all,
        file = '/home/chengww/data3/project/multiple_myeloma/midel_result_R4/normal_mm/mm_normal_common_all.rds')



##read final mm_normal_all data
mm_normal_all <- 
  readRDS('/home/chengww/data3/project/multiple_myeloma/midel_result_R4/normal_mm/mm_normal_common_all.rds')
obj <- subset(mm_normal_all,celltype_big %in% c('HSPC','Stromal cells','T & NK','B cells','Myeloid cells'))
obj$celltype_big <-factor(obj$celltype_big, levels=c('HSPC','Stromal cells','T & NK','B cells','Myeloid cells'))
tmp.data <- table(obj$Response,obj$celltype_big) %>% apply(1,prop.table) %>% t()
pheatmap::pheatmap(tmp.data,
                   scale = 'none',
                   display_numbers = TRUE,         #热图格子中显示相应的数值
                   number_color = "black",         #字体颜色为黑色
                   fontsize=13,
                   number_format = "%.4f",
                   cluster_row = F, border=FALSE,
                   cluster_cols = FALSE,angle_col = 0,
                   color = colorRampPalette(colors = c("#f8f2ea",'#d08e36','#bf812d',"#aa7328",
                                                       "#744e1b"))(100),
                   legend = F)



##Heatmap-----------------------------------------------------------------------
Idents(mm_normal_all)<- 'celltype_big'
find.gene <- FindAllMarkers(mm_normal_all,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
data.table::fwrite(find.gene,'/home/chengww/data3/project/multiple_myeloma/midel_result_R4/normal_mm/normal_and_mm_find.genes.csv')
#top2 <- find.gene %>% group_by(cluster) %>% filter(avg_log2FC>0.25 & p_val_adj<0.05& pct.1>0.2 & pct.2 <0.3)
#top2 <- top2[duplicated(top2$gene),]
top2 <- find.gene %>% group_by(cluster) %>% top_n(5,avg_log2FC)
##细胞数有16万，画不出来热图，所以要抽取样本
set.seed(1995)
subobj <- subset(mm_normal_all, downsample = 300)#抽样
subobj <- ScaleData(subobj,features =rownames(subobj))#标准化抽样样本
DoHeatmap(subobj,features = c(top2$gene[1:33],'MZB1','IGHG4','SDC1'),size=4,draw.lines=F,group.by = 'celltype_big',
                group.colors = bigcol,angle = 40) + NoLegend()+
  scale_fill_gradientn(colors= colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(100))
ggsave(p1,filename = '/home/chengww/data/project/multiple_myeloma/midel_result_R4/normal_mm/fig.pdf',width = 7,height = 9)
rm(subobj)
##
DefaultAssay(mm_normal_all) <-'ADT'
Idents(mm_normal_all)<- 'celltype_big'
sub<- subset(mm_normal_all,Response!='HD')
find.gene <- FindAllMarkers(sub,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
#data.table::fwrite(find.gene,'/home/chengww/data3/project/multiple_myeloma/midel_result_R4/normal_mm/normal_and_mm_find.genes.csv')
#top2 <- find.gene %>% group_by(cluster) %>% filter(avg_log2FC>0.25 & p_val_adj<0.05& pct.1>0.2 & pct.2 <0.3)
#top2 <- top2[duplicated(top2$gene),]
top2 <- find.gene %>% group_by(cluster) %>% top_n(5,avg_log2FC)
##细胞数有16万，画不出来热图，所以要抽取样本
set.seed(1995)
subobj <- subset(sub, downsample = 300)#抽样
subobj <- ScaleData(subobj,features =rownames(subobj))#标准化抽样样本
DoHeatmap(subobj,features = c(top2$gene),size=4,draw.lines=F,group.by = 'celltype_big',
          group.colors = bigcol,angle = 40) + NoLegend()+
  scale_fill_gradientn(colors= colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(100))


##绘图看样本组差异
DimPlot(mm_normal_all,reduction = 'umap',label = F,split.by = 'Response',
        ncol = 1,
        cols =bigcol,repel = T,raster = F)+
  #labs(title = 'Group')+
  theme_bw() +
  gg.theme+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank())+
  scale_y_continuous(breaks = NULL)+
  NoLegend()

DimPlot(mm_normal_all,reduction = 'umap',label = F,group.by = 'Gender',
        cols =patie.col[c(20,9)],# c("#F1784B", "#889D5A")
        raster = T)+
  theme_bw() +
  theme(legend.position = c(0.1,0.84),
        panel.border =element_rect(color = 'black',fill = NA,size=1),
        axis.text.x=element_text(),
        axis.text.y=element_text(),
        text = element_text(size=12,face="bold"),
        plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(), 
        strip.text = element_text(),
        strip.background = element_blank(),
        axis.line = element_blank())+
  labs(title = 'Gender')
mm_normal_all$orig.ident <-ifelse(mm_normal_all$orig.ident=='normal BM','HD',mm_normal_all$orig.ident)
DimPlot(mm_normal_all,reduction = 'umap',label = F,group.by = 'orig.ident',
       cols =patie.col[c(10,6)],# c("#F1784B", "#889D5A")
       raster = T)+
  theme_bw() +
  theme(legend.position = c(0.1,0.84),
        panel.border =element_rect(color = 'black',fill = NA,size=1),
        axis.text.x=element_text(),
        axis.text.y=element_text(),
        text = element_text(size=12,face="bold"),
        plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_blank(), 
        strip.text = element_text(),
        strip.background = element_blank(),
        axis.line = element_blank())+
  labs(title = 'Source')

DimPlot(mm_normal_all,reduction = 'umap',label = F,group.by = 'Patient',
        cols = patie.col,
        raster = T)+
  theme_bw() +
  theme(#legend.position = c(0.1,0.88),
    panel.border =element_rect(color = 'black',fill = NA,size=1),
    axis.text.x=element_text(),
    axis.text.y=element_text(),
    text = element_text(size=12,face="bold"),
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_blank(), 
    strip.text = element_text(),
    strip.background = element_blank(),
    axis.line = element_blank())+
  labs(title = 'Patients')




data <- mm_normal_all@assays$RNA@data[c('IL1B','TNF'),]  %>% 
  as.data.frame() %>% t() %>% as.data.frame()
data$celltype <- mm_normal_all$celltype_big
data$Reponse <- mm_normal_all$Response
head(data)
comp <- list(c('T & NK','B cells'),
             c('Myeloid cells','MMPC')
)
p3 <- ggplot(subset(data,TNF!=0),aes(x=celltype,y=TNF,fill=celltype))+
  geom_violin(trim=T,scale='width')+
  geom_boxplot(width=0.1,position = position_dodge(0.9),
               outlier.shape = NA)+
  scale_fill_manual(values = bigcol)+
  ggpubr::stat_compare_means(comparisons = comp,label = 'p.signif')+
  gg.theme+NoLegend()+
  labs(x='',y='Expression levels',title = 'TNF')+
  ggpubr::stat_compare_means()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))
p4 <-ggplot(subset(data,IL1B!=0),aes(x=celltype,y=IL1B,fill=celltype))+
  geom_violin(trim=T,scale='width')+
  geom_boxplot(width=0.1,position = position_dodge(0.9),
               outlier.shape = NA)+
  scale_fill_manual(values = bigcol)+
  ggpubr::stat_compare_means(comparisons = comp,label = 'p.signif')+
  ggpubr::stat_compare_means()+
  labs(x='',y='Expression levels',title = 'IL1B')+
  gg.theme+NoLegend()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))
p4|p3


##part4 专门画图------------------------------------
##########柱状图展示各群样本分布
gg.theme <- theme_bw() +
  theme(#legend.position = "right",
    axis.text.x=element_text(angle=0,hjust=0.5),
    text = element_text(size=12,face="bold"),
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(), 
    strip.text = element_text(size=12,face="bold"),
    strip.background = element_blank(),
    axis.line = element_line(color = 'black'))

table(mm_normal_all$celltype_big,mm_normal_all$Response) %>% as.data.frame() %>%
  ggplot(aes(x=Var2,y=Freq,fill=Var1))+
  geom_bar(stat = 'identity',position = 'fill',width = 0.8)+
  gg.theme+
  labs(x='',y='Number of cells',fill='')+
  scale_fill_manual(values = bigcol)+
  theme(axis.text.x = element_text(size = 10,angle = 60,vjust = 1,hjust = 1))
##统计单个样本的细胞亚群比例
obj <- subset(mm_normal_all,celltype_big != 'nPC' & celltype_big != 'MMPC' )
library(tibble)
library(tidyr)
tmp_data <- table(obj$celltype_big,obj$Patient)
tmp_data <-apply(tmp_data,2,prop.table) %>% as.data.frame() %>% rownames_to_column(.,'celltype_big')
tmp_data <-gather(tmp_data, Patient,Freq,-celltype_big)
colnames(tmp_data) <- c('celltype_big','Patient','Freq')
tmp_data <- left_join(tmp_data,obj[[]] %>% dplyr::select(Patient,Response) %>%unique(),by=c('Patient'))
tmp_data <- tmp_data[order(tmp_data$Response),]
tmp_data$Patient <- paste(tmp_data$Response,tmp_data$Patient,sep = '-')
tmp_data$celltype_big <- factor(tmp_data$celltype_big,levels =celltype_big_levles[-c(7:8)] )
ggplot(tmp_data,aes(x=Patient,y=Freq,fill=celltype_big))+
  geom_bar(stat='identity',width = 0.8,position = 'fill')+
  labs(x='',title='Sample',fill = "")+
  scale_fill_manual(values = bigcol[-c(7:8)])+
  gg.theme+
  #facet_wrap(~Response,ncol = 1)+
  theme(axis.text.x=element_text(angle=90,hjust=1))


##样本组细胞数量比例-柱状图
resp.data <- table(obj$Response,obj$celltype_big)
resp.data <- apply(resp.data, 1, prop.table) %>%as.data.frame() %>%
  rownames_to_column(var = 'celltype_big') 
resp.data <- gather(resp.data, Response,Freq,-celltype_big)
resp.data$celltype_big<-factor(resp.data$celltype_big,levels = celltype_big_levles[-c(7:8)])
resp.data$Response <- factor(resp.data$Response,levels = c('CR','NR','HD'))
ggplot(resp.data %>% na.omit(),
       aes(x=Response,y=Freq,fill=celltype_big,stratum=celltype_big,alluvium=celltype_big))+
  #geom_bar(stat='identity',width = 0.8,position = 'fill')+
  geom_col(width = 0.5,color=NA)+
  ggalluvial::geom_flow(width=0.5,alpha=0.4,knot.pos=0)+
  labs(x='',title='',fill = "")+
  scale_fill_manual(values = bigcol[-c(7:8)])+
  gg.theme


##各亚群的细胞数量统计
tmp <- table(mm_normal_all$celltype_big) %>%as.data.frame()
ggplot(tmp,aes(x=reorder(Var1,Freq),y=Freq,fill=Var1))+
  geom_text(aes(label=Freq, y=Freq+4500), position=position_dodge(1), vjust=0.5) +
  geom_bar(stat='identity',width = 0.6)+
  labs(x='celltype',title='cell number of each cluster',fill='')+
  scale_fill_manual(values = bigcol)+
  gg.theme+NoLegend()+theme(axis.text.x = element_text(size=10,angle = 0,vjust = 0.5, hjust=0.5))+
  coord_flip()

##样本组间细胞比例差异检验
resp.data <- table(obj$Patient,obj$celltype_big)
resp.data <- apply(resp.data, 1, function(x) x/sum(x)) %>% as.data.frame()
resp.data <- resp.data %>% rownames_to_column(var = 'celltype_big') %>% 
  pivot_longer( cols =  c("A006003":"W"),
                names_to = 'Patient',
                values_to = 'freq')
resp.data <- left_join(resp.data,obj@meta.data %>% dplyr::select('Patient','Response','orig.ident') %>%unique(),by = 'Patient')
resp.data$celltype_big <-factor(resp.data$celltype_big, levels=celltype_big_levles[-c(7:8)])
p1<-ggpubr::ggboxplot(resp.data, x = "celltype_big", y = "freq", color = "Response",palette=resp.col)+
  ggpubr::stat_compare_means(aes(group=Response),label = "p.signif")+#p.format
  gg.theme+
  theme(legend.position = "top",axis.text.x=element_text(angle=60,hjust=1))+
  labs(title = '',x='',y='Fraction',color='')
p1
##多组差异显著性比较
group <-list(c('CR','NR'),c('CR','normal BM'),c('NR','normal BM'))
ggpubr::ggboxplot(resp.data, x = "Response", y = "freq", 
                  color = "Response",facet.by = "celltype_big",palette=resp.col,ncol=20)+
  ggpubr::stat_compare_means(comparisons = group,label = "p.signif")+
  labs(x='',title = '',color='')+
  theme(legend.position = "right",#p.signif
        axis.text.x=element_text(angle=60,hjust=1),
        text = element_text(size=12,face="bold"),
        plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.text.y = element_blank(), 
        strip.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        axis.line = element_line(color = 'black'))


##某个基因的组间差异检验
plot<- list()
for (i in as.character(unique(mm_normal_all$celltype_big))){
  obj<-subset(mm_normal_all,celltype_big==i)
  exp_gata3_skin <- data.frame(CD38.exp=obj@assays$RNA@data['SDC1',],
                              # CNV.score=mm_normal_all$cnv_score,
                               celltype_big = obj@meta.data$celltype_big,
                               Response=obj@meta.data$Response)
  exp_gata3_skin <- exp_gata3_skin[which(exp_gata3_skin$CD38.exp!=0),]
  plot[[i]] <- ggplot(exp_gata3_skin,aes(Response,CD38.exp,fill=Response)) + 
    geom_violin(size=0.1,trim=T)+
    xlab('')+
    #ggforce::geom_sina(size=0.2,alpha=0.2)+
    geom_boxplot(width=0.1,fill="white")+ #绘制箱线图
    ggpubr::stat_compare_means(comparisons = group, method = 't.test',label = "p.signif")+
    # geom_signif(comparisons = list(c('CR','NR')),
    #             y_position=c(seq(4,6,0.3)), map_signif_level=TRUE)+
    scale_fill_manual(values=resp.col)+
    labs(title = paste('CD38 exp in ',i,sep = '')) + theme_classic()+NoLegend()+gg.theme
}
plot[[1]]
cowplot::plot_grid(plot[[1]],plot[[2]],plot[[4]],plot[[5]],plot[[6]],plot[[7]],
                   ncol = 3)

#检查个别基因的组间差异——例如CD38
exp_gata3_skin <- data.frame(CD38.exp=mm_normal_all@assays$RNA@data['CD38',],
                             celltype_big=mm_normal_all$celltype_big,
                             Response=mm_normal_all@meta.data$Response)
exp_gata3_skin <- exp_gata3_skin[which(exp_gata3_skin$Response!='HD' &
                                         exp_gata3_skin$celltype_big=='Myeloid cells'  ),]
#exp_gata3_skin <- exp_gata3_skin[which(exp_gata3_skin$CD38.exp!=0),]
ggplot(exp_gata3_skin, aes(x=factor(Response, level = c("CR", "NR")), y=CD38.exp)) + 
  geom_boxplot(outlier.size = 0.4, aes(fill = Response)) +
  facet_wrap(vars(factor(celltype_big, levels = celltype_big_levles)), 
             scales = "free_y", ncol = 8) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4) +
  labs(x='Response group')+
  theme_bw() +
  theme(strip.text.x = element_text(size = 10)) +
  theme(axis.text.x = element_text(size=10))+
  theme(axis.text.y = element_text(size=10))+ 
  ggpubr::stat_compare_means(size=4,label = 'p.signif')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values =resp.col)+
  theme(legend.position = "none",text = element_text(size=12,face="bold"),
        strip.text = element_text(size=12,face="bold"),
        axis.line = element_line(color = 'black'))
DefaultAssay(mm_normal_all)<-'RNA'
p1<-VlnPlot(mm_normal_all,features = c("CD38"),
        split.by = 'celltype_big',pt.size = 0)&
  xlab('')&
  #geom_boxplot(width=0.1,position = position_dodge(0.9))& #绘制箱线图
  scale_fill_manual(values=bigcol)&
  #ggpubr::stat_compare_means( label = "p.signif")&
  labs(title ='CD38') &
  gg.theme&
  theme(axis.text.x=element_text(angle=60,hjust=1))
p2<-VlnPlot(mm_normal_all,features = c("SDC1"),
            split.by = 'celltype_big',pt.size = 0)&
  xlab('')&
  #geom_boxplot(width=0.1,position = position_dodge(0.9))& #绘制箱线图
  scale_fill_manual(values=bigcol)&
  #ggpubr::stat_compare_means( label = "p.signif")&
  labs(title ='SDC1') &
  gg.theme&
  theme(axis.text.x=element_text(angle=60,hjust=1))
p1+NoLegend()|p2


##计算细胞亚群的相关性
mm_normal_all$resp_celltype_big <- ifelse(mm_normal_all$Response=='CR',
                                   paste('CR_',mm_normal_all$celltype_big,sep = ''),
                                   paste('NR_',mm_normal_all$celltype_big,sep = ''))
Idents(mm_normal_all) <-'celltype_big' #'resp_celltype_big'
tmp_data <- AverageExpression(mm_normal_all,group.by='celltype_big',assays = 'RNA')
tmp_data <- tmp_data[[1]]
tmp_gene <- names(tail(sort(apply(tmp_data, 1, sd)),1000))
#View(tmp_data[tmp_gene,])
cor <- cor(tmp_data[tmp_gene,],method = 'spearman')
# annotation <- MM_obj@meta.data[,c('Patient','Response')] %>% unique() 
# rownames(annotation)<- annotation$Patient
# annotation <- annotation %>% select(Response)
pheatmap::pheatmap(cor,border=FALSE)#,#annotation_col = annotation,color = colorRampPalette(c("navy", "white", "firebrick3"))(500)

##infercnv
subobj <- subset(mm_normal_all,celltype_big %in% c('nPC','MMPC'))
table(subobj$celltype_big)
subobj$celltype_big <-as.character(subobj$celltype_big)
saveRDS(subobj,'/home/chengww/data/project/multiple_myeloma/midel_result_R4/normal_mm/infercnv/myeloma_sub__infercnv.RDS') 


###AddModuleScore计算免疫得分------------------1---------------------------
library(Seurat)
library(RColorBrewer)
library(ggplot2)
##炎性得分
inflammatory_genes <- read.csv('/home/chengww/data3/project/multiple_myeloma/code/inflammatory_genes.csv')
head(inflammatory_genes)
inflammatory_genes <- as.list(inflammatory_genes)
obj <- mm_normal_all
DefaultAssay(obj) <-'RNA'
obj <- AddModuleScore(object = obj,features = inflammatory_genes,
                      ctrl = 100, name = 'inflammatory_genes_score')
colnames(obj@meta.data)[ncol(obj@meta.data)] <-'inflammatory_score'


featurePlotCols=c("lightgrey","lightgrey","lightgrey","#ffffcc","#ffeda0",
                  "#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026","#800026")
FeaturePlot(obj,features = 'inflammatory_score',order = T,max.cutoff = 0.21,split.by = 'Response')&
  scale_color_gradientn(colours =featurePlotCols)&
  gg.theme&NoLegend()
FeaturePlot(obj,features = 'inflammatory_score',order = T,max.cutoff = 0.21)&
  scale_color_gradientn(colours =featurePlotCols)&
  gg.theme&NoLegend()

tmp_data <- obj@meta.data[,c('Patient','Response','celltype_big','inflammatory_score')]
group <-list(c('CR','NR'),c('CR','HD'),c('NR','HD'))
ggplot(tmp_data,aes(Response,inflammatory_score)) + 
  geom_violin(aes(fill=Response),size=0.1,trim=T)+
  xlab('')+
  geom_boxplot(width=0.1,outlier.alpha = 0)+ #绘制箱线图
  ggpubr::stat_compare_means(comparisons =group, method = 't.test',label = "p.signif")+
  scale_fill_manual(values=resp.col)+
  theme_classic()+NoLegend()+gg.theme

tmp_data<-subset(tmp_data,Response!='HD')
tmp_data$Patient_resp <- paste0(tmp_data$Response,'_',tmp_data$Patient)

ggplot(tmp_data,aes(Patient_resp,inflammatory_score)) + 
  geom_violin(aes(fill=Patient_resp),size=0.1,trim=T)+
  xlab('')+
  geom_boxplot(width=0.1,outlier.alpha = 0)+ #绘制箱线图
  #ggpubr::stat_compare_means(comparisons =group, method = 't.test',label = "p.signif")+
  scale_fill_manual(values=patie.col)+
  theme_classic()+labs(title = 'Inflammation score in all cells')+
  theme(legend.position = "right",#p.signif
        axis.text.x=element_text(angle=60,hjust=1),
        text = element_text(size=12,face="bold"),
        plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.text.y = element_blank(), 
        strip.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        axis.line = element_line(color = 'black'))+NoLegend()

ggplot(subset(tmp_data,celltype_big!='MMPC'),aes(Patient_resp,inflammatory_score)) + 
  geom_violin(aes(fill=Patient_resp),size=0.1,trim=T)+
  xlab('')+
  geom_boxplot(width=0.1,outlier.alpha = 0)+ #绘制箱线图
  #ggpubr::stat_compare_means(comparisons =group, method = 't.test',label = "p.signif")+
  scale_fill_manual(values=patie.col)+
  theme_classic()+labs(title = 'Inflammation score in immune cells')+
  theme(legend.position = "right",#p.signif
        axis.text.x=element_text(angle=60,hjust=1),
        text = element_text(size=12,face="bold"),
        plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.text.y = element_blank(), 
        strip.text = element_text(size=12,face="bold"),
        strip.background = element_blank(),
        axis.line = element_line(color = 'black'))+NoLegend()


figs<-list()
for (i in unique(tmp_data$celltype_big)){
  tmp_data2<-subset(tmp_data,celltype_big==i)
  figs[[i]]<-ggplot(tmp_data2,aes(Patient_resp,inflammatory_score)) + 
    geom_violin(aes(fill=Patient_resp),size=0.1,trim=T)+
    xlab('')+
    geom_boxplot(width=0.1,outlier.alpha = 0)+ #绘制箱线图
    #ggpubr::stat_compare_means(comparisons =group, method = 't.test',label = "p.signif")+
    scale_fill_manual(values=patie.col)+
    theme_classic()+labs(title = paste('Inflammation score in',i))+
    theme(legend.position = "right",#p.signif
          axis.text.x=element_text(angle=60,hjust=1),
          text = element_text(size=12,face="bold"),
          plot.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          strip.text.y = element_blank(), 
          strip.text = element_text(size=12,face="bold"),
          strip.background = element_blank(),
          axis.line = element_line(color = 'black'))+NoLegend()
}
cowplot::plot_grid(
                   figs[[4]],
                   figs[[6]],
                   figs[[7]],
                   ncol = 1)


VlnPlot(obj,features = 'inflammatory_score',pt.size = 0,
      cols = resp.col,group.by = 'Response')+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA)+
  NoLegend()+
  ggpubr::stat_compare_means(label = "p.signif")+
  theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust = 1))+
  labs(x='')
##细胞因子得分
DefaultAssay(obj)
go_genes1 <- getGO("GO:0001819")#positive regulation of cytokine production
obj <- AddModuleScore(obj, features = go_genes1, name = "positive regulation of cytokine production")
tmp_data <- obj@meta.data[,c('Patient','Response','celltype_big','positive.regulation.of.cytokine.production1')]
group <-list(c('CR','NR'),c('CR','HD'),c('NR','HD'))
ggplot(subset(tmp_data,celltype_big=='Myeloid cells'),aes(Response,positive.regulation.of.cytokine.production1)) + 
  geom_violin(aes(fill=Response),size=0.1,trim=T)+
  xlab('')+
  geom_boxplot(width=0.1,outlier.alpha = 0)+ #绘制箱线图
  ggpubr::stat_compare_means(comparisons =group, method = 't.test',label = "p.signif")+
  scale_fill_manual(values=resp.col)+
  theme_classic()+NoLegend()+gg.theme
VlnPlot(obj,features = 'positive.regulation.of.cytokine.production1',
        pt.size = 0,cols = resp.col,group.by = 'Response')+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA)+
  NoLegend()+
  ggpubr::stat_compare_means(label = "p.signif")+
  theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust = 1))+
  labs(x='',y='Cytokine score',title = 'Cytokine score')
FeaturePlot(obj,features = 'positive.regulation.of.cytokine.production1',order = T,max.cutoff = 0.21)&
  scale_color_gradientn(colours =featurePlotCols)&
  gg.theme&NoLegend()&
  labs(title = 'Cytokine score')

##interferon socre--------------------------------------------------------------
interferon_genes <- read.csv('/home/chengww/data3/project/multiple_myeloma/midel_result_R4/normal_mm/interferon_genes.txt')
head(interferon_genes)
interferon_genes <- as.list(interferon_genes)

DefaultAssay(obj) <-'RNA'
obj <- AddModuleScore(object = obj,features = interferon_genes,
                      ctrl = 100, name = 'interferon_genes_score')
colnames(obj@meta.data)[ncol(obj@meta.data)] <-'interferon_score'
VlnPlot(obj,features = 'interferon_score',pt.size = 0)+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA)
FeaturePlot(obj,features = 'interferon_score',order = T,max.cutoff = 0.21)&
  scale_color_gradientn(colours =featurePlotCols)&
  gg.theme&NoLegend()&
  labs(title = 'Interferon score')

set.seed(1995)
subobj <- subset(mm_normal_all, downsample = 300)#抽样
subobj <- ScaleData(subobj,features =rownames(subobj))#标准化抽样样本
DoHeatmap(subobj,features = c(interferon_genes$ges),size=4,draw.lines=F,group.by = 'celltype_big',
          group.colors = bigcol,angle = 40) + NoLegend()+
  scale_fill_gradientn(colors= colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(100))

avetmp <- AverageExpression(mm_normal_all,group.by = 'celltype_big',features = c(interferon_genes$ges),slot = 'data')
avetmp <- as.data.frame(t(avetmp$RNA)) %>% scale() %>% t()
head(avetmp)
library(ComplexHeatmap)
##列注释
ha = HeatmapAnnotation(
  celltype=unlist(lapply(colnames(avetmp),function(x) strsplit(x,split = '-')[[1]][1])),
  Response=unlist(lapply(colnames(avetmp),function(x) strsplit(x,split = '-')[[1]][2])),
  col = list(celltype = c(bigcol),
             Response = c('CR'="#35978F",'NR'="#BF812D",'HD'="#85658D")
  ),
  na_col = "black"
)
Heatmap(avetmp,
        col =colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(100),
        cluster_rows = F,
        cluster_columns = F,
        top_annotation = ha,
        show_row_names = T,
        row_names_gp = gpar(fontsize = 10),
        show_column_names = F,
        # row_split = c(rep('CD8 TM',20),rep('CD8 tox',20),
        #               rep('CR',12),rep('NR',4)),
        heatmap_legend_param = list(
          title= "Zscore", title_position = "topcenter", 
          legend_height=unit(4,"cm"), legend_direction="vertical"),
        #column_title = "Loading genes in PC2", 
        column_title_gp = gpar(fontsize = 15, fontface = "bold")
)
##plasm score
pc_gens <- data.frame(pcgenes=c('IGHG1','MZB1','SDC1','CCND1'))
pc_gens <- as.list(pc_gens)
obj <- mm_normal_all
DefaultAssay(obj) <-'RNA'
obj <- AddModuleScore(object = obj,features = pc_gens,
                      ctrl = 100, name = 'pc_gens')

featurePlotCols=c("lightgrey","lightgrey",
                  "#CD5D5C","#C13533",'#9A2A29',"#8a2524",'#6b1d1c')
FeaturePlot(subset(obj,Response!='normal BM'),features = 'pc_gens1',order = F,raster = F)&
  scale_color_gradientn(colours =featurePlotCols)&
  #scale_color_gradient(low = "lightgrey", high = "#C13533")&
  labs(title = 'Tumor score')&
  gg.theme2



