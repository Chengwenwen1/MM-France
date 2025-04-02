nk<- subset(MM_t,celltype %in% c('NKdim','NKbright'))
nk
DefaultAssay(nk)<- 'RNA'
DimPlot(nk,group.by = 'celltype')

FeaturePlot(nk,features = c('NCAM1','FCGR3A'),order = T)
VlnPlot(nk,features = c('FCGR3A','GZMB','NCAM1'))
##FINDALLMARKER
Idents(nk) <- 'celltype'
find.gene <- FindAllMarkers(nk,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
top2 <- find.gene %>% group_by(cluster) %>% top_n(5,avg_log2FC)
DoHeatmap(nk,features = top2$gene,slot ='scale.data',size=5.5,draw.lines=F,angle = 0,
          group.colors =mycol_t) 
###




DefaultAssay(MM_t)<-'RNA'
Idents(MM_t)<-'celltype'
gene.t <- list()
plot.fire <-list()
gene.data <- data.frame()

##CR vs NR
sub_mm <- subset(MM_t,orig.ident == 'MM')
for (i in c("NKdim","NKbright")){
  gene.t[[i]] <- FindMarkers(sub_mm,ident.1 = 'CR',group.by='Response',min.pct = 0.25,logfc.threshold = 0,subset.ident = i)#group.by = 'orig.ident'
  gene.t[[i]]$cluster <- i
  gene.t[[i]]$threshold = factor(ifelse(gene.t[[i]]$p_val_adj < 0.05 & abs(gene.t[[i]]$avg_log2FC) >= 0.25, 
                                        ifelse(gene.t[[i]]$avg_log2FC >= 0.25 ,'CR high','NR high'),'NS'),
                                 levels=c('CR high','NR high','NS'))
  gene.t[[i]]$geneid <- rownames(gene.t[[i]])
  gene.data<- rbind(gene.data,gene.t[[i]])
  
  plot.fire[[i]] <- ggplot(gene.t[[i]],aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
    geom_point()+
    scale_color_manual(values=c( "#35978F", "#BF812D","#808080"))+
    ggrepel::geom_text_repel(
      data = gene.t[[i]][gene.t[[i]]$p_val_adj<0.05 & abs(gene.t[[i]]$avg_log2FC)>0.25,] %>%
        group_by(threshold) %>% top_n(10,abs(avg_log2FC)),
      aes(label = geneid),
      size = 3,
      segment.color = "black", show.legend = FALSE)+
    ylab('-log10 (p-adj)')+
    xlab('log2 (FoldChange)')+
    labs(color='',title = i)+
    geom_vline(xintercept=c(-0.25,0.25),lty=3,col="black",lwd=0.5) +
    geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+
    gg.theme#+theme(legend.position = c(0.6,0.95))
  #saveRDS(plot.fire,'./plot.fire.rds')
}
plot.fire$NKdim+NoLegend() |plot.fire$NKbright


##go
gene.t <- list()
gsea.data <-list()
go.fig <- list()
go <-list()
z=0
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
for (i in levels_ctype_t[12:13]) {#as.character(levels(MM_t$celltype))[1:8]
  gene.t[[i]] <- FindMarkers(sub_mm,ident.1 = 'NR',group.by='Response',logfc.threshold = 0,subset.ident = i)
  gene.t[[i]]$cluster <-i
  gene.t[[i]]$threshold = factor(ifelse(gene.t[[i]]$p_val_adj < 0.05 & abs(gene.t[[i]]$avg_log2FC) >= 0.25, 
                                        ifelse(gene.t[[i]]$avg_log2FC >= 0.25 ,'NR high','CR high'),'NS'),
                                 levels=c('NR high','CR high','NS'))
  gene.t[[i]]$geneid <- rownames(gene.t[[i]])

  gsea.data[[i]] <- gene.t[[i]][gene.t[[i]]$p_val_adj<0.05 & abs(gene.t[[i]]$avg_log2FC) >0.25,] 
  for (t in unique(gsea.data[[i]]$threshold)){
    gene <- gsea.data[[i]]$geneid[which(gsea.data[[i]]$threshold==t)]
    gene_map <- bitr(gene,fromType="SYMBOL", toType=c("ENTREZID"),OrgDb = 'org.Hs.eg.db')
    gene_map <- unique(gene_map)
    #gsea.data[[j]] <- left_join(gsea.data[[j]],gene_map,c('geneid'='SYMBOL')) %>% na.omit()
    gene <-gene_map$ENTREZID
    ego <- enrichGO(gene = gene,
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 2,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    z=paste(i,t,sep = '_')
    go[[z]] <- clusterProfiler::simplify(
      ego,
      cutoff = 0.7,
      by = "p.adjust",
      select_fun = min,
      measure = "Wang",
      semData = NULL
    )
    go.fig[[z]] <- dotplot(go[[z]] ,showCategory=15,title=paste(i,t,sep = '-'))+gg.theme
  }
}

##----NKdim--------------------------------
View(go$`NKdim_CR high`@result)
View(go$`NKdim_NR high`@result)
go$`NKdim_CR high`@result$group <-'CR'
go$`NKdim_NR high`@result$group <-'NR'
go4 <- rbind(go$`NKdim_CR high`@result,go$`NKdim_NR high`@result)
go4 <- go4[go4$Description %in% c('antibody-dependent cellular cytotoxicity',
                                  'leukocyte mediated cytotoxicity',
                                  'positive regulation of lymphocyte activation',
                                  'leukocyte activation involved in inflammatory response',
                                 # 'phagocytosis',
                                  'defense response to bacterium',##CR
                                  'transcription factor AP-1 complex',
                                  'ER to Golgi transport vesicle membrane',
                                  'transcription repressor complex'
                                  
                                  ),]
go4$GeneRatio1 <- lapply(
  go4$GeneRatio,function(x){
    as.numeric(strsplit(x,'/')[[1]][1])/as.numeric(strsplit(x,'/')[[1]][2])}
) %>% unlist() %>% round(digits = 2)
go4$GeneRatio1 <- ifelse(go4$group=='NR',-go4$GeneRatio1,go4$GeneRatio1)
go4 <- go4[order(go4$GeneRatio1,decreasing=F),]
go4$Description <-factor(go4$Description,levels =go4$Description )
go4$text_x <- ifelse(go4$group=='NR',rep(-0.11,5),rep(-0.11,3))
p1 <- ggplot(data = go4,
       aes(x = GeneRatio1, y = Description)) +
  geom_bar(aes(fill = -log10(p.adjust)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "GeneRatio", y = "Description", title = "Different terms in NKdim") +
  geom_text(aes(x = text_x,
                label = Description),
            hjust = 0)+ #hjust=0，
  gg.theme+ theme(axis.text.y = element_blank())
p1
##----NKbright-----------------------------
View(go$`NKbright_CR high`@result)
View(go$`NKbright_NR high`@result)

go$`NKbright_CR high`@result$group <-'CR'
go$`NKbright_NR high`@result$group <-'NR'
go4 <- rbind(go$`NKbright_CR high`@result,go$`NKbright_NR high`@result)
go4 <- go4[go4$Description %in% c(
  'positive regulation of natural killer cell chemotaxis',
  'cytokine-mediated signaling pathway',
  'positive regulation of lymphocyte migration',
  'positive regulation of leukocyte activation',
  #'phagocytosis',
  'positive regulation of leukocyte mediated cytotoxicity',
  'positive regulation of GTPase activity',#CR
  'negative regulation of leukocyte mediated cytotoxicity',
  'negative regulation of lymphocyte mediated immunity',
  'negative regulation of cell killing'
                                  
),]
go4$GeneRatio1 <- lapply(
  go4$GeneRatio,function(x){
    as.numeric(strsplit(x,'/')[[1]][1])/as.numeric(strsplit(x,'/')[[1]][2])}
) %>% unlist() %>% round(digits = 2)
go4$GeneRatio1 <- ifelse(go4$group=='NR',-go4$GeneRatio1,go4$GeneRatio1)
go4 <- go4[order(go4$GeneRatio1,decreasing=F),]
go4$Description <-factor(go4$Description,levels =go4$Description )
go4$text_x <- ifelse(go4$group=='NR',rep(-0.11,6),rep(-0.11,3))
p2 <- ggplot(data = go4,
       aes(x = GeneRatio1, y = Description)) +
  geom_bar(aes(fill = -log10(p.adjust)), stat = "identity", width = 0.8, alpha = 0.7) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1) +
  labs(x = "GeneRatio", y = "Description", title = "Different terms in NKbright") +
  geom_text(aes(x = text_x, 
                label = Description),
            hjust = 0)+ #hjust=0,
  gg.theme+ theme(axis.text.y = element_blank())
p1|p2

##
VlnPlot(subset(MM_t,Response!='HD'& celltype %in% c('NKdim','NKbright') ),
        features = c('CD226','KLRK1','NCR1'),#CD64,CD32
        pt.size=0,split.by = 'Response',cols = resp.col)&
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA)&
  ggpubr::stat_compare_means(label.x = 2,label = "p.signif")&
  labs(x='')&gg.theme
FeaturePlot(subset(MM_t,Response!='HD'& celltype %in% c('NKdim','NKbright') ),
        features = c('NCAM1','CD226','KLRK1','NCR1'))
