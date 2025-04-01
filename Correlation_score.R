

obj <- nk
Idents(obj) <- 'celltype'
DefaultAssay(obj) <-'RNA'
#GO:0006909 #ADCP
go_genes1 <- getGO("GO:0006909")#ADCP
obj <- AddModuleScore(obj, features = go_genes1, name = "ADCP")
colnames(obj@meta.data)[ncol(obj@meta.data)] <-'ADCP'
FeaturePlot(obj,features = c("ADCP"),order = T,raster = T,
            cols =c("lightgrey",'#BE766E','#C13533'))+
  labs(title ='ADCP score' )+
  gg.theme


#GO:0001788 #"antibody-dependent cellular cytotoxicity"
go_genes1 <- getGO("GO:0001788")#"antibody-dependent cellular cytotoxicity"
obj <- AddModuleScore(obj, features = go_genes1, name = "ADCC")
colnames(obj@meta.data)[ncol(obj@meta.data)] <-'ADCC'
FeaturePlot(obj,features = c("ADCC"),order = T,raster = T,
            cols =c("lightgrey",'#BE766E','#C13533'))+
  labs(title ='ADCC score' )+
  gg.theme



###AddModuleScore计------------------1---------------------------
inflammatory_genes <- read.csv('*/inflammatory_genes.csv')
head(inflammatory_genes)
inflammatory_genes <- as.list(inflammatory_genes)
DefaultAssay(obj) <-'RNA'
obj <- AddModuleScore(object = obj,features = inflammatory_genes,slot ='scale.data',
                      ctrl = 100, name = 'inflammatory_genes_score')
colnames(obj@meta.data)[ncol(obj@meta.data)] <-'inflammatory_score'
VlnPlot(obj,features = 'inflammatory_score',pt.size = 0,cols = mycol_t)+
  geom_boxplot(width=0.1,position = position_dodge(0.9),outlier.shape = NA)+
  ggpubr::stat_compare_means()+
  gg.theme+
  labs(x='',title = 'Inflammation')+
  NoLegend()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))

resp.data <- obj@meta.data[,c('inflammatory_score','ADCC','Response','celltype','Patient')]
resp.data$Response<-factor(resp.data$Response,levels = c('CR','NR','HD'))
# resp.data <-cbind(aggregate(subset(resp.data,celltype %in% c('NKdim','NKbright'))[,'ADCC'],
#                             list(subset(resp.data,celltype %in% c('NKdim','NKbright'))[,'Patient']), FUN=sum),#ADCC Nk patients
#                   aggregate(resp.data$inflammatory_score, list(resp.data$Patient), FUN=sum)
#                   )
# colnames(resp.data) <-c('Patient','ADCC','name2','inflammatory_score')
# resp.data <- left_join(resp.data, obj@meta.data[,c('Response','Patient')] %>% unique(),by='Patient') 

ggplot(subset(resp.data, Response!='HD'),
       aes(x=inflammatory_score,y=ADCC))+
  geom_point(alpha=0.2,size=1)+
  gg.theme+
  scale_color_manual(values = c('CR'="#35978F",'NR'="#BF812D",'HD'="#85658D"))+  
  geom_smooth(method = lm,level=0.9)+
  stat_cor(method = 'pearson')+#label.x = 0.4,label.y = 0.4
  
  labs(title = 'Correlation for ADCC and inflammatory score')

resp.data <- obj@meta.data[,c('inflammatory_score','ADCC','Response','celltype','Patient')]
resp.data$CD16<-obj@assays$RNA@data['FCGR3A',]
resp.data$Response<-factor(resp.data$Response,levels = c('CR','NR','HD'))
ggplot(subset(resp.data,celltype=='NKdim'& Response!='HD'),
       aes(x=CD16,y=inflammatory_score))+
  geom_point(alpha=0.1,size=1)+
  gg.theme+
  scale_color_manual(values = c('CR'="#35978F",'NR'="#BF812D",'HD'="#85658D"))+  
  geom_smooth(method = lm,level=0.9)+
  stat_cor(method = 'pearson')+#label.x = 0.4,label.y = 0.4
  
  labs(title = 'Correlation for inflammatory_score and CD16 in NKdim')

##
obj$celltype <- as.character(obj$celltype)
tmp_data <- table(obj$celltype,obj$Patient)
tmp_data <-apply(tmp_data,2,prop.table) %>% as.data.frame() %>% rownames_to_column(.,'celltype')
tmp_data <-tidyr::gather(tmp_data, Patient,Freq,-celltype)
colnames(tmp_data) <- c('celltype','Patient','Freq')
tmp_data <- left_join(tmp_data,obj[[]] %>% dplyr::select(Patient,Response) %>%unique(),by=c('Patient'))
tmp_data <- tmp_data[order(tmp_data$Response),]
tmp_data$Patient <- paste(tmp_data$Response,tmp_data$Patient,sep = '-')
tmp_data$Patient <- factor(tmp_data$Patient,levels = unique(tmp_data$Patient))
#tmp_data$celltype <-factor(tmp_data$celltype, levels=levels_ctype_t)
tmp_data$Response <-factor(tmp_data$Response)
ggplot(tmp_data, aes(x = celltype, y = Freq, fill = Response)) +
  #facet_wrap(~celltype,ncol = 13)+
  geom_bar(stat = "summary", fun = mean, color = "black", position = position_dodge()) +
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.25,position = position_dodge( .9))+
  stat_compare_means(method = 'anova',method.args = list(var.equal = F),label = "p.signif",label.y = 0.45)+
  scale_fill_manual(values = resp.col)+
  labs(fill='')+
  gg.theme+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))#
##
head(obj@meta.data)
score_data <- obj@meta.data %>% 
  group_by(celltype, Response,Patient) %>% 
  summarize(aa=sum(inflammatory_score)) %>% as.data.frame()
score_data$Patient <- paste(score_data$Response,score_data$Patient,sep = '-') 
type ='NKdim'#NKbright
score_data_dim <- left_join(score_data %>% subset(celltype==type),
                            tmp_data %>% subset(celltype==type) %>% dplyr::select(c(Patient,Freq)))
score_data_dim <-score_data_dim %>% subset(Response!='HD')
ggplot(score_data_dim,aes(x=aa,y=Freq, colour =Response))+
  geom_point(alpha=1,size=1)+
  gg.theme+
  scale_color_manual(values = c('CR'="#35978F",'NR'="#BF812D",'HD'="#85658D"))+  
  geom_smooth(method = lm)+
  stat_cor(method = 'pearson')+#label.x = 0.4,label.y = 0.4
  labs(title = 'Correlation for NKdimFreq and inflammatory score',
       x='inflammatory score',y='NKdimFreq')



resp.data1 <- obj@meta.data[,c('inflammatory_score','ADCP','Response','celltype','Patient')]
resp.data1$Response<-factor(resp.data1$Response,levels = c('CR','NR','HD'))
#resp.data1<- subset(resp.data1,celltype %in% levels_ctype)
ggplot(resp.data1,aes(x=inflammatory_score,y=ADCP, colour =Response))+
  geom_point(alpha=0.1,size=1)+
  gg.theme+
  scale_color_manual(values = c('CR'="#35978F",'NR'="#BF812D",'HD'="#85658D"))+  
  geom_smooth(method = lm)+
  stat_cor(method = 'pearson')+#label.x = 0.4,label.y = 0.4
  
  labs(title = 'Correlation for ADCP and inflammatory score')

##inteferon
library(data.table)
library(Seurat)
interferon_genes <- read.csv('*/interferon_genes.txt')
head(interferon_genes)
interferon_genes <- as.list(interferon_genes)
DefaultAssay(obj) <-'RNA'
obj <- AddModuleScore(object = obj,features = interferon_genes,
                      ctrl = 100, name = 'interferon_genes_score')
colnames(obj@meta.data)[ncol(obj@meta.data)] <-'interferon_score'



resp.data <- obj@meta.data[,c('interferon_score','ADCC','Response','celltype','Patient')]
resp.data$ITGB1<-obj@assays$RNA@data['ITGB1',]
resp.data$Response<-factor(resp.data$Response,levels = c('CR','NR','HD'))
ggplot(subset(resp.data,celltype=='NKdim'),aes(x=ITGB1,y=ADCC, colour =Response))+
  geom_point(alpha=0.1,size=1)+
  gg.theme+
  scale_color_manual(values = c('CR'="#35978F",'NR'="#BF812D",'HD'="#85658D"))+  
  geom_smooth(method = lm,level=0.9)+
  stat_cor(method = 'pearson')+#label.x = 0.4,label.y = 0.4
  
  labs(title = 'Correlation for ADCC and ITGB1')


resp.data1 <- obj@meta.data[,c('interferon_score','ADCP','Response','celltype','Patient')]
resp.data1$KLF6<-obj@assays$RNA@data['KLF6',]#'FCGR1A','FCGR2A','FCER1G'
resp.data1$Response<-factor(resp.data1$Response,levels = c('CR','NR','HD'))
resp.data1<- subset(resp.data1,celltype %in% levels_ctype)
ggplot(resp.data1,aes(x=KLF6,y=ADCP, colour =Response))+
  geom_point(alpha=0.1,size=1)+
  gg.theme+
  scale_color_manual(values = c('CR'="#35978F",'NR'="#BF812D",'HD'="#85658D"))+  
  geom_smooth(method = lm)+
  stat_cor(method = 'pearson')+#label.x = 0.4,label.y = 0.4
  labs(title = 'Correlation for ADCP and FCGR2A')

