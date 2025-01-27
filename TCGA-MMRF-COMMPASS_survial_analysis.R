setwd("/data3/chengww/project/multiple_myeloma/midel_result_R4/normal_mm/COMMPASS")

#install.packages("rjson")
library("rjson")
json <- jsonlite::fromJSON("./metadata.cart.2023-03-08.json")
View(json)
#id <- json$associated_entities[[1]][,1]
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#获取gdc_download文件夹下的所有TSV表达文件的 路径+文件名
count_file <- list.files('./rna_seq_data/',pattern = '*.tsv',recursive = TRUE)
#在count_file中分割出文件名
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[1]})
matrix = data.frame(matrix(nrow=60660,ncol=0))
for (i in 1:length(count_file)){
  path = paste0('./rna_seq_data/',count_file[i])
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  #data <- data[3]   #取出unstranded列（第3列），即count数据，对应其它数据
  data <- data[6] #第六列是TPM
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}
#write.csv(matrix,'COUNT_matrix.csv',row.names = TRUE)##生成count文件
#write.csv(matrix,'COUNT_matrix.csv',row.names = TRUE)##生成tpm文件
##增加部分：设置Gene Symbol为列名的矩阵（前面得到的是Ensembl ID）------------------------------------------
path = paste0('./rna_seq_data/',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
matrix0 <- cbind(gene_name,matrix)
#将gene_name列去除重复的基因，保留每个基因最大表达量结果
matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)    
#将gene_name列设为行名
rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-1]
#write.csv(matrix0,'COUNT_matrix_symbol.csv',row.names = TRUE)##生成count文件
write.csv(matrix0,'TPM_matrix_symbol.csv',row.names = TRUE)##生成tpm文件

#分为normal和tumor矩阵--------------------------
sample <- colnames(matrix0)+
normal <- c()
tumor <- c()
for (i in 1:length(sample)){
  if((substring(colnames(matrix0)[i],14,15)>10)){    #14、15位置大于10的为normal样本
    normal <- append(normal,sample[i])
  } else {
    tumor <- append(tumor,sample[i])
  }
}

tumor_matrix <- matrix0[,tumor]
normal_matrix <- matrix0[,normal]

#临床数据整合
#install.packages("rjson")
library("rjson")
json <- jsonlite::fromJSON("./metadata.cart.2023-03-08.json")
View(json)
entity_submitter_id <- sapply(json$associated_entities,function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
sample_case <- t(rbind(entity_submitter_id,case_id))
clinical <- read.delim('./clinical.tsv',header = T)
clinical <- as.data.frame(clinical[duplicated(clinical$case_id),])
clinical_matrix <- merge(sample_case,clinical,by="case_id",all.x=T)
clinical_matrix <- clinical_matrix[,-1]
clinical_matrix <- clinical_matrix[!duplicated(clinical_matrix$entity_submitter_id),]
fwrite(clinical_matrix,'./clinical_matrix.csv')##生存信息文件


####从这里开始----TCGA---------------------------------------------------------------
matrix0 <- read.csv('./COUNT_matrix_symbol.csv',row.names = 1)
clinical_matrix <- read.csv('./clinical_matrix.csv')
##做一个简单的生存曲线
library(survival)
library(survminer)
clinical_matrix$OS_MONTHS <- round(clinical_matrix$days_to_last_follow_up/30,2)
table(clinical_matrix$vital_status)
table(clinical_matrix$OS_MONTHS >0)
table(clinical_matrix$gender)
all.genes <-unique(c('CLEC7A','CCR7','TRIM8','CARD9','SLC44A2','MAP3K3','FLOT1','TLR2',#bulk
                     'BCL10','NOD2','ALPK1','RIPK1','MAP3K14','CFLAR','RELA',#pseudo bulk
                     'RPS3','EDN1','SPHK1','CD27','PELI1','BIRC3','F2R','TSPAN6','CD74',#sc pc1
                     'RPS3','RHOA','CD27','EDN1','LAPTM5','RELA','TERF2IP','PDCD4',
                     'BIRC3','EEF1D','CFLAR','CD74','PELI1','F2R','CD40','TMEM9B',
                     'BIRC2','VAPA','TIFA'#sc degs
))


myeloma <- left_join(t(matrix0[all.genes,]) %>% as.data.frame() %>% mutate(id=rownames(t(matrix0))) ,
                     clinical_matrix %>% select(OS_MONTHS,vital_status,entity_submitter_id) %>% na.omit(),
                     by=c('id'='entity_submitter_id') ) %>% na.omit()
myeloma$vital_status <- ifelse(myeloma$vital_status=='Dead',1,0)

##所有基因加和的分组结果生存分析
test <- data.frame(all=apply(myeloma[,1:34],1, sum),myeloma[,35:37])
cutoff <- surv_cutpoint(test, #数据集
                        time = "OS_MONTHS", #生存状态
                        event = "vital_status", #生存时间
                        variables =  'all'
)
highname <- test$id[test$all > summary(cutoff)[1,1]]
test$group <- ifelse(test$id %in% highname,'high','low')
surv <- Surv(test$OS_MONTHS,test$vital_status==1)
kmfit1 <- survfit(surv~test$group)
ggsurvplot(kmfit1,data = test,pval = T)

###各个基因单独计算显著性
cutoff <- surv_cutpoint(myeloma, #数据集
                        time = "OS_MONTHS", #生存时间
                        event = "vital_status", #生存状态
                        variables =  all.genes
)
my.surv <- Surv(myeloma$OS_MONTHS,myeloma$vital_status==1)
p<-list()
for (i in all.genes){
  cutoff_value <- summary(cutoff)[i,1]
  highname <- myeloma$id[myeloma[,i]> cutoff_value]
  myeloma[,i] <- ifelse(myeloma$id %in% highname,'high','low')
  myeloma$group <- ifelse(myeloma$id %in% highname,'high','low')
  kmfit1 <- survfit(my.surv~myeloma$group)
  p[[i]] <- ggsurvplot(kmfit1,data = myeloma,pval = T)+labs(title = i)
}
p

# cowplot::plot_grid(p[['RIPK1']]$plot,p[['TLR2']]$plot,p[['FLOT1']]$plot,
#                    p[['BCL10']]$plot,p[['CLEC7A']]$plot,p[['TRIM8']]$plot,
#                    p[['MAP3K3']]$plot,p[['RHOA']]$plot,p[['TMEM9B']]$plot,
#                    p[['BIRC2']]$plot,p[['VAPA']]$plot,
#                    p[['TIFA']]$plot,p[['TERF2IP']]$plot,ncol = 7)

#添加多个参数,定义函数
UniCox<-function(x){
  FML<-as.formula(paste0("my.surv~",x))
  Cox<-coxph(FML,data = myeloma) 
  Sum<-summary(Cox)
  CI.low<-round(Sum$conf.int[,3:4],2)[1]
  CI.high<-round(Sum$conf.int[,3:4],2)[2]
  Pvalue<-round(Sum$coefficients[,5],3)
  HR<-round(Sum$coefficients[,2],2)
  Unicox<-data.frame("Characteristics"=x,
                     "Hazard Ratio"=HR,
                     "95CI.low"=CI.low,
                     "95CI.high"=CI.high,
                     'HR2'= paste(HR," (",CI.low,"-",CI.high,")",sep=""),
                     stringsAsFactors = F,
                     "P value"=Pvalue)
  return(Unicox)
}

#检查所有基因
Univar <-lapply(all.genes, UniCox) 
Univar<-plyr::ldply(Univar,data.frame)
Univar <- Univar[Univar$P.value<0.05 & Univar$Hazard.Ratio<1,]
Univar <- Univar[order(Univar$P.value),]
Univar[nrow(Univar)+1,c(5,6)] <- c('HR2','P.value')
forestplot::forestplot(labeltext=Univar,Univar[,c(1,6,5)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                       mean=Univar[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                       lower=Univar[,3],  #告诉函数表格第3列为5%CI，
                       upper=Univar[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                       zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                       boxsize=0.5,       #设置小黑块的大小
                       ##fpColors函数设置颜色
                       col=forestplot::fpColors(box="#1c61b6", lines="#1c61b6", zero = "lightgray"),
                       #箱线图中基准线的位置
                       cex=0.9, lineheight = "auto",
                       colgap=unit(8,"mm"),
                       #箱子大小，线的宽度
                       lwd.ci=2, 
                       #箱线图两端添加小竖线，高度
                       ci.vertices=TRUE, ci.vertices.height = 0.4,
                       graph.pos=4)       #森林图应插在图形第2列


###选取Lenalidomide用药的病人组
meta <- subset(clinical_matrix,therapeutic_agents %in% c('Lenalidomide','Daratumumab'))#'Daratumumab',
matrix <- matrix0 %>% dplyr::select(meta$entity_submitter_id)
# p<-list()
# for (i in all.genes){
#   my.surv <- Surv(meta$OS_MONTHS,meta$vital_status=='Dead')
#   highname <- colnames(matrix)[matrix[i,]> mean(as.numeric(matrix[i,]))]
#   meta$group <- ifelse(meta$entity_submitter_id %in% highname,
#                        'high','low')
#   kmfit1 <- survfit(my.surv~meta$group)
#   p[[i]] <- ggsurvplot(kmfit1,data = meta,pval = T)+labs(title = i)
# }
# cowplot::plot_grid(p[['RELA']]$plot,p[['TLR2']]$plot,ncol = 2)
myeloma_drug<- left_join(t(matrix[all.genes,]) %>% as.data.frame() %>% mutate(id=rownames(t(matrix))) ,
                         meta %>% select(OS_MONTHS,vital_status,entity_submitter_id) %>% na.omit(),
                         by=c('id'='entity_submitter_id') ) %>% na.omit()

myeloma_drug$vital_status <- ifelse(myeloma_drug$vital_status=='Dead',1,0)
cutoff <- surv_cutpoint(myeloma_drug, #数据集
                        time = "OS_MONTHS", #生存状态
                        event = "vital_status", #生存时间
                        variables =  all.genes
)
my.surv <- Surv(myeloma_drug$OS_MONTHS,myeloma_drug$vital_status==1)
p<-list()
for (i in all.genes){
  cutoff_value <- summary(cutoff)[i,1]
  highname <- myeloma_drug$id[myeloma_drug[,i]> cutoff_value]
  myeloma_drug[,i] <- ifelse(myeloma_drug$id %in% highname,'high','low')
  myeloma_drug$group <- ifelse(myeloma_drug$id %in% highname,'high','low')
  kmfit1 <- survfit(my.surv~myeloma_drug$group)
  p[[i]] <- ggsurvplot(kmfit1,data = myeloma_drug,pval = T)+labs(title = i)
}
cowplot::plot_grid(p[['RIPK1']]$plot,p[['TLR2']]$plot,p[['BCL10']]$plot,
                   p[['FLOT1']]$plot,p[['MAP3K3']]$plot,p[['CARD9']]$plot,
                   p[['CLEC7A']]$plot,
                   p[['RELA']]$plot,p[['SPHK1']]$plot,p[['PELI1']]$plot,
                   p[['RHOA']]$plot,p[['TERF2IP']]$plot,
                   p[['BIRC2']]$plot,p[['VAPA']]$plot,ncol = 7)

#添加多个参数,定义函数
UniCox<-function(x){
  FML<-as.formula(paste0("my.surv~",x))
  Cox<-coxph(FML,data = myeloma_drug) 
  Sum<-summary(Cox)
  CI.low<-round(Sum$conf.int[,3:4],2)[1]
  CI.high<-round(Sum$conf.int[,3:4],2)[2]
  Pvalue<-round(Sum$coefficients[,5],3)
  HR<-round(Sum$coefficients[,2],2)
  Unicox<-data.frame("Characteristics"=x,
                     "Hazard Ratio"=HR,
                     "95CI.low"=CI.low,
                     "95CI.high"=CI.high,
                     'HR2'= paste(HR," (",CI.low,"-",CI.high,")",sep=""),
                     stringsAsFactors = F,
                     "P value"=Pvalue)
  return(Unicox)
}

#如查看
Univar <-lapply(all.genes, UniCox) 
Univar<-plyr::ldply(Univar,data.frame)
Univar <- Univar[Univar$P.value<0.05 & Univar$Hazard.Ratio<1,]
Univar <- Univar[order(Univar$P.value),]
Univar$P.value <- round(as.numeric(Univar$P.value),3)
Univar[nrow(Univar)+1,c(5,6)] <- c('HR2','P.value')
forestplot::forestplot(labeltext=Univar,Univar[,c(1,6,5)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                       mean=Univar[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                       lower=Univar[,3],  #告诉函数表格第3列为5%CI，
                       upper=Univar[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                       zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                       boxsize=0.5,       #设置小黑块的大小
                       ##fpColors函数设置颜色
                       col=forestplot::fpColors(box="#1c61b6", lines="#1c61b6", zero = "lightgray"),
                       #箱线图中基准线的位置
                       cex=0.9, lineheight = "auto",
                       colgap=unit(8,"mm"),
                       #箱子大小，线的宽度
                       lwd.ci=2, 
                       #箱线图两端添加小竖线，高度
                       ci.vertices=TRUE, ci.vertices.height = 0.4,
                       graph.pos=4)       #森林图应插在图形第2列


###自己的数据生存分析-----------------------------------------------------------
library(survival)
library(survminer)
head(bulk.counts)##表达矩阵
clinical.data <- fread('/data/jinwf/chengww/project/multiple_myeloma/raw_data/clinical_21_patients.csv')
head(clinical.data)#临床矩阵
my.surv <- Surv(clinical.data$OSmonths,clinical.data$DEATH==1)
p<-list()
for (i in all.genes){
  highname <- colnames(bulk.counts)[bulk.counts[i,]> mean(as.numeric(bulk.counts[i,]))]
  clinical.data$group <- ifelse(clinical.data$SUBJID %in% substr(highname,4,10),'high','low')
  kmfit1 <- survfit(my.surv~clinical.data$group)
  p[[i]] <- ggsurvplot(kmfit1,data = clinical.data,pval = T)+labs(title = i)
}
cowplot::plot_grid(p[['RIPK1']]$plot,
                   p[['BCL10']]$plot,p[['NFKBIB']]$plot,p[['BCL3']]$plot,ncol = 4)

t <- coxph(my.surv~clinical.data$group, data=clinical.data)
summary(t)

##临床+矩阵数据
meta <- clinical_matrix %>% dplyr::select('entity_submitter_id','OS_MONTHS','vital_status','therapeutic_agents')
matrix1 <- t(matrix0) %>% as.data.frame() %>% rownames_to_column('entity_submitter_id')
survival.data <- left_join(meta,matrix1,by='entity_submitter_id') %>%na.omit()
head(survival.data2[,1:10])
survival.data2 <- apply(survival.data[,5:ncol(survival.data)],2,function(x){ifelse(x > median(x),'high','low')})
survival.data <-cbind(survival.data,survival.data2)
survival.data$vital_status
##按照基因的高低进行排序

##2.批量单因素cox回归分析
covariates <- c('MNDA','HK1','LILRA5','PYHIN1','CD36','MEFV')
#分别对每一个变量，构建生存分析的公式
univ_formulas <- sapply(sc.terms.genes,
                        function(x) as.formula(paste('Surv(survival.data$OS_MONTHS, survival.data$vital_status==\'Dead\')~', x)))
#循环对每一个特征做cox回归分析
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = survival.data)})
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)



##3.多因素cox回归分析
my.surv <- Surv(clinical_matrix$OS_MONTHS,clinical_matrix$vital_status=='Dead')
p<-list()
for (i in all.genes){
  highname <- colnames(matrix0)[ matrix0[i,]> mean(as.numeric(matrix0[i,])) ]
  clinical_matrix[,i] <-ifelse(clinical_matrix$entity_submitter_id %in% highname,'high','low')
  clinical_matrix$group <- ifelse(clinical_matrix$entity_submitter_id %in% highname,'high','low')
  kmfit1 <- survfit(my.surv~clinical_matrix$group)
  p[[i]] <- ggsurvplot(kmfit1,data = clinical_matrix,pval = T)+labs(title = i)
}

my.surv <- Surv(clinical_matrix$OS_MONTHS,clinical_matrix$vital_status=='Dead')
clinical_matrix$status <- ifelse(clinical_matrix$vital_status=='Dead',1,0)
res.cox <- coxph(Surv(OS_MONTHS, status) ~ RIPK1+TLR2+FLOT1+BCL10+CLEC7A+
                   TRIM8+ MAP3K3+ RHOA+ TMEM9B+ TSPAN6+ BIRC2+ VAPA+ TIFA+ TERF2IP,
                 data =  clinical_matrix)
x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(Gene = gsub('low','',names(pvalue)),
                     exp=HR,
                     low.95=low,
                     high.95=high,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F,
                     p.value=pvalue
)
multi_res
forestplot::forestplot(multi_res[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                       mean=multi_res[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                       lower=multi_res[,3],  #告诉函数表格第3列为5%CI，
                       upper=multi_res[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                       zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                       boxsize=0.3,       #设置小黑块的大小
                       graph.pos=2)       #森林图应插在图形第2列



