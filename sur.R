library(survival)
library(survminer)
# library(lance)
library(readxl)
library(readr)
library(tidyverse)
library(Seurat)



dat = readRDS('01.rawdata/fpkm.rds')
# dat = readRDS('01.rawdata/tcga_fpkm.rds')
phe = readRDS('01.rawdata/phe.rds')
# phe = readRDS('01.rawdata/tcga_phe.rds')

# risk2 = readRDS('03.cox&lasso/riskScore.RDS')
# 

phe$OS.time = phe$OS.time*365
phe = phe[which(phe$OS.time > 0),]

colnames(dat) = colnames(dat) %>% substr(start = 1,stop = 12)
colnames(dat) = gsub(colnames(dat),replacement = '-',pattern = '[.]')

usedSample = intersect(phe$bcr_patient_barcode,colnames(dat))

dat = dat[,usedSample]
dat = dat %>% na.omit()
phe = phe[which(phe$bcr_patient_barcode %in% usedSample),]

# usesample = intersect(colnames(dat),phe$bcr_patient_barcode)
# phe = phe[which(phe$bcr_patient_barcode %in% usesample),]
# dat = dat[,phe$bcr_patient_barcode]


# seurat ------------------------------------

scobj = readRDS('02.seurat/scobj.rds')
all.marker = FindAllMarkers(scobj,group.by = 'seurat_clusters',logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)

usedGene = all.marker$gene %>% unique()
usedGene =intersect(rownames(dat),usedGene)
surlist = list()
survobj <- with(phe, Surv(OS.time,OS))
for(i in 1:length(usedGene)){
  
  # the_tf = colnames(exprTF)[i]
  # the_df = data.frame(row.names = rownames(exprTF),Level = exprTF[,the_tf])
  
  phe$expr = dat[usedGene[i],phe$bcr_patient_barcode] %>% as.numeric()
  phe$Level = ifelse(phe$expr > median(phe$expr),'High','Low')
  
  if(sum(phe$expr) == 0){
    next
  }
  
  fit0 <- survfit(survobj~ Level,data=phe)
  
  p_value <- surv_pvalue(fit0)$pval
  if (p_value < 0.05) {
    surlist[[usedGene[i]]] = fit0
  }
  
  

  # p = surv_pvalue(fit0, method = "log-rank")$pval %>% print()
  # if(p < 0.05){
  # 
  #   the_df = data.frame(
  #     row.names = rownames(exprTF),
  #     Level = exprTF[,the_tf],
  #     OS = clinical$OS.time
  #   )
  # 
  #   highT = median(the_df$OS[which(the_df$Level == 'High')])
  #   lowT =  median(the_df$OS[which(the_df$Level == 'Low')])
  #   if(highT > lowT){
  #     Time_diff[the_tf,"clincal"] = 'safe'
  #   }else{
  #     Time_diff[the_tf,"clincal"] = 'dangerous'
  #   }
  # 
  # }
  
}

all.marker = all.marker[which(all.marker$gene %in% names(surlist)),]

all.marker %>%
  group_by(cluster) %>%
  mutate(
    specificity_score = (avg_log2FC * 2) 
  ) %>%
  top_n(n = 20, wt = specificity_score) -> alltop_markers




usedGene = "HOXA5,MEIS1,HOXA3,HOXB3,HOXB6,HK3,CD1D,LGALS2,HNMT,CD163,CD109,INPP4B,DLL1,HOPX,FNBP1L,ANKRD9,CLIC2,OPTN,SLC14A1,HBG2,DEFA1B,CEACAM8,S100A12,MMP8,ARG1,PAPSS2,CD1C,SULF2,CLEC5A,S100A9,POU4F1,PALM,TRH,CACNA2D2,TPSAB1,ARPP21,SHD,MYO18B,TRAT1,TRH,FGF13,HGF,TDRD9,LIN7A,KCNE5"
usedGene = usedGene %>% str_split(',') %>% unlist()

pdf(file = '02.seurat/markerGene.pdf',width = 6,height = 4.5)
DoHeatmap(scobj,group.by = 'seurat_clusters',features = usedGene,size = 3,group.colors = ggsci::pal_aaas()(10) )
dev.off()









# 
# usedGene = c('PTGDS',
#              'IGFBP2',
#              'HGF',
#              'LGALS12',
#              'FGF13',
#              'MYO18B',
#              'TSPAN7',
#              'HPGDS',
#              'UMODL1',
#              'ARPP21',
#              'PRAME',
#              'TRH',
#              'CAV1',
#              'POU4F1',
#              'ROBO1',
#              'TGFBI',
#              'MSLN',
#              'CLEC10A',
#              'MYH11',
#              'DEFA1B',
#              'CEACAM8',
#              'S100A12',
#              'MMP8',
#              'ARG1',
#              'HBM',
#              'SELENBP1',
#              'ALAS2',
#              'RHD',
#              'EPB42',
#              'DNTT',
#              'BAALC',
#              'SHANK3',
#              'MN1',
#              'SDK2',
#              'CD14',
#              'LGALS2',
#              'IL1RN',
#              'CD1D',
#              'HK3',
#              'HOXA3',
#              'COL4A5',
#              'NKX2.3',
#              'C3orf80',
#              'HOXA5')

usedGene = alltop_markers$gene
surlist = list()
survobj <- with(phe, Surv(OS.time,OS))
for(i in 1:length(usedGene)){
  
  # the_tf = colnames(exprTF)[i]
  # the_df = data.frame(row.names = rownames(exprTF),Level = exprTF[,the_tf])
  
  phe$expr = dat[usedGene[i],phe$bcr_patient_barcode] %>% as.numeric()
  phe$Level = ifelse(phe$expr > median(phe$expr),'High','Low')
  fit0 <- survfit(survobj~ Level,data=phe)
  
  
  p_value <- surv_pvalue(fit0)$pval
  if (p_value < 0.05) {
    surlist[[usedGene[i]]] = fit0
  }
  
  
  p = ggsurvplot(fit0,
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = F, # 添加风险表
                 risk.table.col = "strata", # 根据分层更改风险表颜色
                 linetype = "strata", # 根据分层更改线型
                 surv.median.line = "hv", # 同时显示垂直和水平参考线
                 ggtheme = theme_bw(), # 更改ggplot2的主题
                 palette = c("#E7B800", "#2E9FDF"))
  p = p$plot + 
    theme(axis.text = element_text(color = 'black'))+
    coord_fixed(max(phe$OS.time)/1)+
    xlab('Time (days)') + ggtitle(usedGene[i])
  p
  ggsave(filename = paste('05.clinicla/surplot/',usedGene[i],'.pdf',sep = ''),plot = p,width = 4,height = 4)
  
  # p = surv_pvalue(fit0, method = "log-rank")$pval %>% print()
  # if(p < 0.05){
  # 
  #   the_df = data.frame(
  #     row.names = rownames(exprTF),
  #     Level = exprTF[,the_tf],
  #     OS = clinical$OS.time
  #   )
  # 
  #   highT = median(the_df$OS[which(the_df$Level == 'High')])
  #   lowT =  median(the_df$OS[which(the_df$Level == 'Low')])
  #   if(highT > lowT){
  #     Time_diff[the_tf,"clincal"] = 'safe'
  #   }else{
  #     Time_diff[the_tf,"clincal"] = 'dangerous'
  #   }
  # 
  # }
  
}




library(survival)
library(dplyr)
library(tibble)

# 确保 phe 中有匹配的样本 ID，且顺序与 dat 列一致
phe <- phe %>% filter(bcr_patient_barcode %in% colnames(dat))
dat <- dat[, phe$bcr_patient_barcode, drop = FALSE]

# 检查
stopifnot(identical(colnames(dat), phe$bcr_patient_barcode))

run_univ_cox <- function(gene_expr, time, event, gene_name = NA) {
  x <- as.numeric(unlist(gene_expr, use.names = FALSE))
  if (length(unique(x[!is.na(x)])) < 2) return(NULL)  # 常数基因跳过
  
  df <- data.frame(
    time  = as.numeric(time),
    event = as.numeric(event),
    expr  = x
  )
  
  fit <- tryCatch(
    coxph(Surv(time, event) ~ expr, data = df),
    error = function(e) { message(gene_name, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(fit)) return(NULL)
  
  s <- summary(fit)
  data.frame(
    gene     = gene_name,
    HR       = s$conf.int["expr", "exp(coef)"],
    HR_lower = s$conf.int["expr", "lower .95"],
    HR_upper = s$conf.int["expr", "upper .95"],
    p_value  = s$coefficients["expr", "Pr(>|z|)"],
    n        = fit$n,
    events   = fit$nevent
  )
}

cox_results <- lapply(rownames(dat), function(g) {
  res <- run_univ_cox(dat[g, ], phe$OS.time, phe$OS)
  if (!is.null(res)) res$gene <- g
  res
}) %>% bind_rows()

# 多重检验校正
cox_results <- cox_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

head(cox_results)



g <- rownames(dat)[1]
gene_expr <- dat[g, ]
str(gene_expr)
length(gene_expr)
length(phe$OS.time)
length(phe$OS)

# 直接拟合看真实报错
df <- data.frame(time = phe$OS.time, event = phe$OS, expr = as.numeric(dat[g, ]))
fit <- coxph(Surv(time, event) ~ expr, data = df)
summary(fit)


res <- run_univ_cox(dat[g, ], phe$OS.time, phe$OS, gene_name = g)
str(res)




cox_results <- lapply(alltop_markers$gene, function(g) {
  run_univ_cox(dat[g, , drop = FALSE], phe$OS.time, phe$OS, gene_name = g)
}) %>% bind_rows()

cox_results <- cox_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

nrow(cox_results)

alltop_markers$HR = cox_results$HR[match(alltop_markers$gene,cox_results$gene)]
alltop_markers$coxPvalue = cox_results$p_value[match(alltop_markers$gene,cox_results$gene)]
alltop_markers = alltop_markers[which(alltop_markers$coxPvalue < 0.05),]

alltop_markers$GeneRisk = ifelse(alltop_markers$HR > 1,'Risk gene','Protective gene')
write.csv(alltop_markers,file = '05.clinicla/alltopmarkers.csv')

result = table(alltop_markers$cluster,alltop_markers$GeneRisk) %>% as.matrix() %>% as.data.frame()

table(alltop_markers$cluster,alltop_markers$GeneRisk)
df_wide <- as.data.frame(table(alltop_markers$cluster, alltop_markers$GeneRisk)) %>%
  setNames(c("cluster", "direction", "n")) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(
    total     = `Protective gene` + `Risk gene`,
    risk_prop = `Risk gene` / total
  )

df_wide = df_wide[,c(1:3)]
write.csv(df_wide,file = '05.clinicla/result.csv')
