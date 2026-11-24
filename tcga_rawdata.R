library(TCGAbiolinks)
library(SummarizedExperiment )
library(stringr)
library(RNAAgeCalc)



TCGAbiolinks:::getGDCprojects()$project_id

data.category = "Gene expression"
#下载rna-seq的counts数据
data.type = "Gene Expression Quantification"
# #下载miRNA数据
# data.type = "miRNA Expression Quantification"
# #下载Copy Number Variation数据
# data.type = "Copy Number Segment

project <- "TCGA-LAML"

# 构建查询
query <- GDCquery(
  project = project,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)


GDCdownload(query)

GDCdownload(query)
data <- GDCprepare(query)
expr = assay(data)

ensembleID = rownames(expr) %>% str_split('[.]')
ensembleID = lapply(ensembleID, function(tab){
  return(tab[1])
}) %>% unlist()

anno = readRDS('E:/shengao/R/Genecode_anno.Hs.Rds')
symbles = anno$symbol[match(ensembleID,anno$geneID)]
expr = expr %>% as.data.frame()
expr$gene = symbles

expr = expr %>% na.omit()
symbles = expr$gene
expr = expr[,1:(ncol(expr)-1)]
expr = as.matrix(expr)
rownames(expr) = symbles

query_clinical <- GDCquery(
  project = project,
  data.category = "Clinical", 
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab" # 这是一种格式
)

# 下载数据
GDCdownload(query_clinical, method = "api")

# 准备临床数据
clinical <- GDCprepare(query_clinical)

# clinical 现在是一个列表，包含多个相关的临床表格
names(clinical)

# 最重要的两个表格通常是：
# 1. 患者信息表
clinical_patient <- clinical$clinical_patient_laml




df<-subset(clinical_patient,select =c(bcr_patient_barcode,vital_status,death_days_to,last_contact_days_to))
df = df[3:nrow(df),]

df$OS.time = 0
for(i in 1:nrow(df)){
  
  if(df$vital_status[i] == 'Dead'){
    df$OS.time[i] = df$death_days_to[i]
  }else{
    df$OS.time[i] = df$last_contact_days_to[i]
  }
  
}

df = df[which(df$OS.time != '[Not Available]'),]
df$OS.time = df$OS.time %>% as.numeric()

df[df$vital_status=='Dead',]$vital_status <- 1
df[df$vital_status=='Alive',]$vital_status <- 0

df$OS = ifelse(df$vital_status == 'Dead',1,0)

phe = subset(clinical_patient,select =c(submitter_id,age_at_diagnosis,ann_arbor_clinical_stage,gender))

df = merge(df,phe,'submitter_id')

colnames(expr) = substr(colnames(expr),start = 1,stop = 12)
fpkm = count2FPKM(expr)

saveRDS(expr,file = '01.rawdata/tcga_expr.rds')
saveRDS(df,file = '01.rawdata/tcga_phe.rds')
saveRDS(fpkm,file = '01.rawdata/tcga_fpkm.rds')

