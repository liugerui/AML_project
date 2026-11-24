library(GEOquery)
library(Biobase)
library(tidyverse)
library(dplyr)
library(RTN)
library(RTN)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(snow)


# dat_filepath = '00rawdata/02.GSE13213.csv'
# group_filepath = '00rawdata/02.GSE13213.group.csv'
# phe_filepath = '00rawdata/02.GSE13213.clinical.csv'

#从GEO获取数据

gset<-getGEO(filename = '01.rawdata/GSE6891_series_matrix.txt.gz',getGPL = F)
expr<-as.data.frame(gset@assayData$exprs)
phe = gset@phenoData@data

gpl<-getGEO(filename = '01.rawdata/GPL570.annot.gz',destdir = '.')
gpl<-Table(gpl)

colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene symbol')%>%
  filter('Gene symbol'!='')%>%
  separate('Gene symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)

dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

dat = dat[1:((nrow(dat)*0.9) %>% ceiling()),]#保留表达值前90%的基因

group = data.frame(
  
  sample = phe$geo_accession,
  group = phe$source_name_ch1
  
)

saveRDS(dat,file = '01.rawdata/expr.rds')




