library(org.Hs.eg.db)
library(Seurat)
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
library(DESeq2)



dat =  readRDS('01.rawdata/expr.rds') %>% as.matrix()
scobj = readRDS('02.seurat/scobj.rds')

group = scobj@meta.data
group$group = 'High risk'
group$group[which(group$seurat_clusters %in% c('6','7','8'))] = 'Low risk'

table(group$group)
control = 'Low risk'
case = 'High risk'

group$group = factor(group$group, levels = c(control,case))
design <- cbind(control = ifelse(group$group == control, 1, 0), 
                case = ifelse(group$group == case, 1, 0))

contrast.mat = makeContrasts(contrasts='case-control', levels=design)
fit = lmFit(dat, design)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
DEG_test=na.omit(fit)
saveRDS(DEG_test,'06.DEGs/01.limma_out.rds')

############# identify DEGs
Res <- readRDS('06.DEGs/01.limma_out.rds')
logFC_x = 0.5
P_x = 0.05

Res$regulate = 'Normal'
Res$regulate[which(Res$logFC > logFC_x& Res$P.Value<P_x)] = 'Up'
Res$regulate[which(Res$logFC < (-logFC_x) & Res$P.Value<P_x)] = 'Down'

table(Res$regulate)
Res$gene = rownames(Res)

Res = Res[order(Res$regulate,Res$logFC),]

saveRDS(Res,'06.DEGs/01.limma_out.rds')
write.csv(Res,'06.DEGs/01.limma_out.csv')

LR = read.csv(
  '06.DEGs/gene_input.csv'
)
LR = LR$gene_name

sig_diff <- subset(Res,
                   Res$regulate != 'Normal')
# dat_rep<- rbind(head(sig_diff,100),tail(sig_diff,100))

usedGene = c('FTH1','TNFSF13','EREG','LILRB1','LILRB2',
             'CD34','EPCAM','VEGFA','PTPRCAP','CD244' )
dat_rep = sig_diff[which(sig_diff$gene %in% usedGene),]

maxfc = max(abs(Res$logFC))

volcano_plot<- ggplot(data = Res, 
                      aes(x = logFC,
                          y = -log10(P.Value), 
                          color =regulate)) +
  
  
  geom_point(size = 2)+
  geom_point(data = dat_rep,size = 4,aes(x = logFC,
                                         y = -log10(P.Value)), color ='black')+
  
  geom_point(data = dat_rep,size = 3,aes(x = logFC,
                                         y = -log10(P.Value)), color ='white')+
  
  geom_point(data = dat_rep,size = 2,aes(x = logFC,
                                         y = -log10(P.Value), color = regulate))+
  
  scale_color_manual(values = c( "#3B4992FF","darkgray","#EE0000FF")) +
  
  theme_bw(base_size = 12) +
  geom_vline(xintercept = c(-logFC_x,logFC_x),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(P_x),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(
          color="black",
          size=10),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(
          color = "black",
          size = 10),
        axis.text.y = element_text(
          color = "black",
          size = 10),
        axis.title.x = element_text(
          color = "black",
          size = 10),
        axis.title.y = element_text(
          color = "black",
          size = 10)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 3,
    box.padding = unit(0.3, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.3, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "Log2(Fold Change)")+
  scale_x_continuous(limits = c(-maxfc,maxfc))+
  # scale_y_continuous(limits = c(0,15))+
  ylab('-Log10(P-value)')+
  theme(legend.position = 'top')+
  coord_fixed(9/80)
volcano_plot

ggsave('06.DEGs//01.volcano_plot.png', volcano_plot,width = 3.74, height = 4.23)
ggsave('06.DEGs/01.volcano_plot.pdf', volcano_plot,width = 3.74, height = 4.23)

# heatmap -------------------------------

sig_diff <- sig_diff[which(sig_diff$gene %in% LR),]

usedGene = c('GCLM', 'GCH1','FTH1', 'ITGAL' ,'SELPLG' ,'SELL' ,'PECAM1' ,'IL4R', 'IL6R','TREM1' ,'LILRB1', 'LILRB2', 'LILRB3','CD33' ,
             'CD34' ,'KIT' ,'MPL' ,'VEGFA' ,'TIE1' ,'PLAU' ,'NCAM1' ,'NTRK1','SLC18A2'  )





group$title = factor(group$group,levels = c('Low risk','High risk'))
group = group[order(group$title,decreasing = F),]
group$title = as.character(group$title)

rowanno <- data.frame(group = group$title)
rownames(rowanno) <- group$cellname

candidate_Gene = usedGene

plot_dat <- dat[candidate_Gene,group$cellname] %>% as.matrix()

# for (i in 1:nrow(plot_dat)) {
# 
#   tmp = plot_dat[i,] %>% scale()
#   plot_dat[i,] = tmp
# 
# }
library(pheatmap)

# yorder = LR$value

# plot_dat = plot_dat[yorder,group$sample]
p2 <- pheatmap(plot_dat,  
               color = colorRampPalette(c('blue','white','red'))(100), 
               border_color = NA,  
               scale = "row", 
               cluster_rows = F, 
               cluster_cols = FALSE, 
               legend = TRUE, 
               # legend_breaks = c(-0.5, 0, 0.5),
               # legend_labels = c("-0.5","0","0.5"), 
               show_rownames = TRUE, 
               show_colnames = FALSE, 
               fontsize = 12,
               annotation_col = rowanno,   #添加列注释信息
               annotation_legend = TRUE,    #是否显示注释信息图例
               annotation_names_row = TRUE, #是否显示行注释的名称
               annotation_names_col = TRUE #是否显示列注释的名称
)
ggsave(filename = '06.DEGs//04.heatmap.pdf',p2,w=6,h=8)
ggsave(filename = '06.DEGs//04.heatmap.png',p2,w=6,h=8)

group$ISRGs = scobj$ISRGs[group$cellname]
group$FE = scobj$FE[group$cellname]
group$IFE = scobj$IFE[group$cellname]



# 首先进行统计检验
stat.test <- group %>%
  t_test(ISRGs ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_xy_position(x = "dose")

library(ggpubr)
p = ggplot(group,aes(x = group,y = FE,fill = group))+
  geom_boxplot(width = 0.5)+
  stat_compare_means(label.y = 2,label="p.signif", method="wilcox.test",label.x = 1.5)+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(color = 'black'),
        legend.position = '')
p  

ggsave('06.DEGs/FE.pdf',p,width = 4,height = 4)

p = ggplot(group,aes(x = group,y = IFE,fill = group))+
  geom_boxplot(width = 0.5)+
  stat_compare_means(label.y = 2,label="p.signif", method="wilcox.test",label.x = 1.5)+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(color = 'black'),
        legend.position = '')
p  

ggsave('06.DEGs/IFE.pdf',p,width = 4,height = 4)


saveRDS(scobj,file = '02.seurat/scobj.rds')


















