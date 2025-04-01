setwd("*/")

#install.packages("rjson")
library("rjson")
json <- jsonlite::fromJSON("./metadata.cart.2023-03-08.json")
View(json)
#id <- json$associated_entities[[1]][,1]
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  

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
  #data <- data[3]  
  data <- data[6] 
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}
#write.csv(matrix,'COUNT_matrix.csv',row.names = TRUE)
#write.csv(matrix,'COUNT_matrix.csv',row.names = TRUE)
##------------------------------------------
path = paste0('./rna_seq_data/',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
matrix0 <- cbind(gene_name,matrix)

matrix0 <- aggregate( . ~ gene_name,data=matrix0, max)    

rownames(matrix0) <- matrix0[,1]
matrix0 <- matrix0[,-1]
#write.csv(matrix0,'COUNT_matrix_symbol.csv',row.names = TRUE)
write.csv(matrix0,'TPM_matrix_symbol.csv',row.names = TRUE)


sample <- colnames(matrix0)+
normal <- c()
tumor <- c()
for (i in 1:length(sample)){
  if((substring(colnames(matrix0)[i],14,15)>10)){  
    normal <- append(normal,sample[i])
  } else {
    tumor <- append(tumor,sample[i])
  }
}

tumor_matrix <- matrix0[,tumor]
normal_matrix <- matrix0[,normal]


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
fwrite(clinical_matrix,'./clinical_matrix.csv')



matrix0 <- read.csv('./COUNT_matrix_symbol.csv',row.names = 1)
clinical_matrix <- read.csv('./clinical_matrix.csv')

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

##
test <- data.frame(all=apply(myeloma[,1:34],1, sum),myeloma[,35:37])
cutoff <- surv_cutpoint(test, 
                        time = "OS_MONTHS", 
                        event = "vital_status", 
                        variables =  'all'
)
highname <- test$id[test$all > summary(cutoff)[1,1]]
test$group <- ifelse(test$id %in% highname,'high','low')
surv <- Surv(test$OS_MONTHS,test$vital_status==1)
kmfit1 <- survfit(surv~test$group)
ggsurvplot(kmfit1,data = test,pval = T)

###
cutoff <- surv_cutpoint(myeloma, 
                        time = "OS_MONTHS", #
                        event = "vital_status", #
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

