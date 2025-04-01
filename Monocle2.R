setwd('*/')
library(monocle)
library(Seurat)
library(dplyr)

mm_pc <- readRDS('*/mm_pc_plasma.RDS')
patie.col_m=c("NR-A052009"="#81571e",#"041-011-NR"="blue","044-004-NR"="blue", 
              "NR-A035004"='#966523',
              "NR-A098004"="#aa7328",
              "NR-A026004"="#bf812d","NR-A173001"="#d08e36", 
              "NR-A059001"="#d49a4b","NR-A068001"="#d9a55f", "NR-A064005"="#eabc81",
              "CR-A166007"="#42bdb3",# "015-005-CR"="red","171-002-CR"="green",
              "CR-A026009"="#3caaa1","CR-A026008"="#35978f",
              "CR-A029001"="#2e847d", "CR-A026007"="#28716b","CR-A099001"="#215e59" 
              #"006-003-CR"="red","015-003-CR"="green", "147-001-CR"="red",
)
gg.theme <-
  theme_bw() +
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
DimPlot(mm_pc,reduction = 'umap',cols = patie.col_m,group.by ='newgroup2')+
  gg.theme

obj <- mm_pc
##data for monocle
RA_matrix <- as(as.matrix(obj@assays$RNA@counts),'sparseMatrix')
p_data <-obj@meta.data
f_data<-data.frame(row.names=rownames(obj),gene_short_name=rownames(obj))

##构建CDS对象
pd<-new("AnnotatedDataFrame", data =p_data)
fd<-new("AnnotatedDataFrame", data =f_data)
cds <- newCellDataSet(RA_matrix,
                      phenoData =pd,
                      featureData =fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily=negbinomial.size())

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds)


expressed_genes <- VariableFeatures(obj) ##gene set1 
##gene set2
# disp_table <- dispersionTable(cds) # 挑有差异的
# unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.5) # 挑表达量不太低的
# expressed_genes <- unsup_clustering_genes$gene_id

cds <- setOrderingFilter(cds, expressed_genes)  
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
#使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点
p1<- plot_cell_trajectory(cds,color_by="Pseudotime", cell_size =0.3,show_backbone=TRUE)+
  gg.theme 
p2 <- plot_cell_trajectory(cds,color_by="celltype", cell_size =0.3,show_backbone=TRUE)+
  scale_color_manual(values = patie.col) + 
  theme(legend.position = "top")
p3 <- plot_cell_trajectory(cds, color_by = "Response",cell_size =0.3,show_backbone=TRUE)+
  scale_color_manual(values = resp.col) + 
  theme(legend.position = "top")+
  facet_wrap(~Response, nrow = 2)+
  gg.theme+NoLegend()
p1/p3
p2
plot_cell_trajectory(cds, color_by = "State",cell_size =0.3,show_backbone=TRUE)+
  facet_wrap(~State, nrow = 1)+gg.theme
saveRDS(cds,'./cds.Rdata')

##统计样本在各个分支的比例
tmp <- table(cds$newgroup2,cds$State) %>% as.matrix()
tmp <- apply(tmp, 1, function(x){x/sum(x)})
tmp_anno <- obj@meta.data %>% dplyr::select('newgroup2','SCISS') %>% unique()
rownames(tmp_anno) <- tmp_anno$newgroup2
tmp_anno<- tmp_anno[colnames(tmp),]##排序
anno.col<- list(newgroup2=c("NR-A052009"="#81571e",#"041-011-NR"="blue","044-004-NR"="blue", 
                            "NR-A035004"='#966523',
                            "NR-A098004"="#aa7328",
                            "NR-A026004"="#bf812d","NR-A173001"="#d08e36", 
                            "NR-A059001"="#d49a4b","NR-A068001"="#d9a55f", "NR-A064005"="#eabc81",
                            "CR-A166007"="#42bdb3",# "015-005-CR"="red","171-002-CR"="green",
                            "CR-A026009"="#3caaa1","CR-A026008"="#35978f",
                            "CR-A029001"="#2e847d", "CR-A026007"="#28716b","CR-A099001"="#215e59" )#,
                #Response=c('CR'="#35978F",'NR'="#BF812D")
                )
labels_row = c("State-1", "State-2", "State-3", "State-4","State-5")
pheatmap::pheatmap(tmp,
                   display_numbers = TRUE,         #热图格子中显示相应的数值
                   number_color = "black",         #字体颜色为黑色
                   fontsize=10,
                   annotation_col = tmp_anno,
                   annotation_colors = anno.col,
                   number_format = "%.2f",
                   gaps_col = c(8),
                   labels_row = labels_row,
                   cluster_row = F, border=FALSE,
                   cluster_cols = FALSE,angle_col = 90,
                   legend = F)

##统计各个分支样本的比例
tmp <- table(cds$Response,cds$State) %>% as.matrix()
tmp <- apply(tmp, 2, function(x){x/sum(x)}) %>% t()
pdf('./pie.pdf')
for (i in 1:nrow(tmp)){
  pie(tmp[i,], border="white",
      labels = c(paste(colnames(tmp),' ',round(100*tmp[i,]/sum(tmp[i,])), "%",sep = '')),
      col = resp.col,family='GB1')
}
dev.off()

##伪时间的亚群密度图
ggplot(pData(cds),aes(Pseudotime,colour=newgroup2,fill=newgroup2))+
  geom_density(bw=0.5,size=1,alpha=0.5)+
  scale_color_manual(values = patie.col_m) +
  scale_fill_manual(values = patie.col_m) +
  gg.theme


####拟时序相关的基因
Time_diff <- differentialGeneTest(cds[expressed_genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- Time_diff[order(Time_diff$qval,decreasing = F),'gene_short_name'][1:150]
p=plot_pseudotime_heatmap(cds[Time_genes,],  show_rownames=T, return_heatmap=T,num_clusters=2,
                        hmcols=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100))

##选择每个cluster的基因，并对基因进行注释
p$tree_row
clusters <- data.frame(cutree(p$tree_row,k=2))
clusters[,1] <- as.character(clusters[,1])
colnames(clusters) <- 'gene_clusters'
table(clusters$gene_clusters)
gene1 <- rownames(clusters)[which(clusters$gene_clusters=='2')]
go.enrich <- clusterProfiler::enrichGO(gene = mapIds(org.Hs.eg.db, 
                                                     keys = gene1,
                                                     keytype = 'SYMBOL',
                                                     column = 'ENTREZID'),
                                       OrgDb = 'org.Hs.eg.db',
                                       ont = 'BP',
                                       pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.05,
                                       readable = T) 

View(go.enrich@result)
#barplot(go.enrich,showCategory = 20)
##cluster1
go.enrich@result$geneID[go.enrich@result$Description %in% 
                          c('humoral immune response',#0.00051
                            'leukocyte migration',#0.0108
                            'positive regulation of cell adhesion'#0.01607
                          )]

c("PRSS2","IGHA2","IGHV3-30","IGHA1","IGHV3-74","PTPRC",
  "JCHAIN","CCL4","CCL3","THY1","BMP5","ITGB7",
  "PRSS2","AGR2","THY1","CD3E","PTPRC") %>%unique()
##cluster2
go.enrich@result$geneID[go.enrich@result$Description %in% 
                          c('positive regulation of lymphocyte activation',#5.150773e-06
                            'I-kappaB kinase/NF-kappaB signaling',#0.01485
                            'phagocytosis, recognition',#6.911105e-07
                            'cytokine-mediated signaling pathway',#0.0488113
                            'immune response-regulating signaling pathway',#2.876486e-06
                            'response to tumor necrosis factor'#0.00788
 )]
c("IGHG2","IGHV2-5","IGHV1-18","IGHM","IGLC3","IGLC2","IGHG4","ID2","IGHV3-23",
  "KLF6","CD74","ZFP36L2","TNFAIP3","PELI1","IGHG1", 
  "IGHG1","TNFAIP3","NFKBIA","RPS3","NOP53","PELI1","IGHG1",
  "CD74","RPS3","PELI1","IGHG1","MAP3K8",               
  "BIRC3","KLF2","ZFP36L2","TNFAIP3","NFKBIA","RPS3","GSTP1",                                                      
  "BIRC3","CD74","TNFAIP3","NFKBIA","TLE1","GSTP1","PELI1"                                                       
   ) %>%unique()


##cluster3
#PRSS2/IGHA2/IGHV3-74/PTPRC/RNASE6/JCHAIN,AGR2/THY1/CD3E/ITGB7/MDK
##cluster4
#IGHG2/IGLC2/IGHG4/PELI1/IGHG1,TXNIP/VIM/NFKBIA/GSTP1/TRIB1
##cluster5
#	CCL4/CCL3/CD9


##选择基因绘制热图
genes<-unique(c('CLEC7A','CCR7','TRIM8','CARD9','SLC44A2','MAP3K3','FLOT1','TLR2',#bulk
                'BCL10','NOD2','ALPK1','RIPK1','MAP3K14','CFLAR','RELA',#pseudo bulk
                'RPS3','EDN1','SPHK1','CD27','PELI1','BIRC3','F2R','TSPAN6','CD74',#sc pc1
                'RPS3','RHOA','CD27','EDN1','LAPTM5','RELA','TERF2IP','PDCD4',
                'BIRC3','EEF1D','CFLAR','CD74','PELI1','F2R','CD40','TMEM9B',
                'BIRC2','VAPA','TIFA'))
time_inflam_gens <- Time_diff[which(Time_diff$qval<0.05),'gene_short_name'][Time_diff[which(Time_diff$qval<0.05),'gene_short_name'] %in% genes]
plot_pseudotime_heatmap(cds[time_inflam_gens,], num_clusters=1, show_rownames=T, return_heatmap=T)
time_inflam_gens
##基因随细胞状态的改变
plot_genes_in_pseudotime(cds[c("EDN1","BIRC3","CD74","LAPTM5","CD27","RPS3",
                               "PELI1","SPHK1","TSPAN6","F2R" ),],color_by = 'newgroup2',ncol = 5)+
  scale_color_manual(values = patie.col_m) +gg.theme

