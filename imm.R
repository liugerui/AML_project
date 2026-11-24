# imm -------------------------


library(tidyverse)
library(dplyr)
library(GSVA)
library(stringi)
library(stringr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(CIBERSORT)
library(IOBR)

control = 'control'
case = 'case'

dat = readRDS(file = '01.rawdata/expr.rds')

scobj = readRDS('02.seurat/scobj.rds')

group = scobj@meta.data
group$group = 'High risk'
group$group[which(group$seurat_clusters %in% c('6','7','8'))] = 'Low risk'

table(group$group)
control = 'Low risk'
case = 'High risk'

dat = dat[,group$cellname]

sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
results <- cibersort(sig_matrix, dat %>% as.matrix(),QN = T) %>% as.data.frame() %>% t()
results = results[1:22,]
results = as.data.frame(results)

res.cibersort <- results
colnames_sum <- colnames(results)
colnames(res.cibersort) <- colnames_sum
res.cibersort$cell_type <- rownames(res.cibersort)
rownames(res.cibersort) <- NULL



## 差异-------
dat.cibersort <- res.cibersort %>% 
  tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.cibersort$sample = dat.cibersort$sample %>% gsub(replacement = '-',pattern = '[.]')
group <- data.frame(sample = group$cellname,title = group$group)
dat.cibersort <- merge(group, dat.cibersort, by = "sample")
###删除所有值均为0的列 purrr::discard(df, ~all(.x == 0)) purrr::keep(df, ~any(.x != 0))
dat.cibersort <- dat.cibersort[,apply(dat.cibersort,2,function(dat.cibersort) !all(dat.cibersort==0))]
dat.cibersort2 <- tidyr::gather(dat.cibersort, ImmuneCell, Score, -c("sample", "title"))

colnames(dat.cibersort2)
stat_cibersort <- dat.cibersort2 %>% 
  group_by(ImmuneCell) %>%
  t_test(Score ~ title) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr   ###"holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
  add_significance("p")
write.csv(stat_cibersort,file = '05.imm///stat.cell.csv')
DE.cibersort<-stat_cibersort[which(stat_cibersort$p<0.05),]
write.csv(DE.cibersort,file = '05.imm///DE.cell.csv')
colnames(dat.cibersort2)

# violin.cibersort<-dat.cibersort2[dat.cibersort2$ImmuneCell%in%stat_cibersort$ImmuneCell[which(stat_cibersort$p<0.05)],]
violin.cibersort<-dat.cibersort2
violin.cibersort$group <- factor(violin.cibersort$title,levels = c(case,control))

violin_plot <- ggplot(violin.cibersort, aes(x=ImmuneCell, y=Score,fill=group)) +
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  # geom_violin(aes(fill=group),trim=FALSE,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  # geom_point(aes(color = group), position=position_dodge(0.9), size = 1) +
  scale_fill_manual(values= c("#EE0000FF","#3B4992FF"), name = "Group")+
  labs(x="", y = "",size=20) +
  stat_compare_means(data = violin.cibersort,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',label.y.npc = 0.99,hide.ns = T,
                     paired = F) +
  theme_bw()+
  labs(y = "Fraction")+
  theme(plot.title = element_text(hjust =0.5,colour="black",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10), 
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        legend.text=element_text( hjust = 0.5,colour="black", size=8),
        legend.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'top')
violin_plot

ggsave(filename = '08.imm///01.oscc_control.pdf',violin_plot,w=12,h=4)
ggsave(filename = '08.imm///01.oscc_control.png',violin_plot,w=9,h=4)



###heatmap

group = group[order(group$title),]

rowanno <- data.frame(group = group$title)
rownames(rowanno) <- group$sample
plot_dat <- results[,group$sample] %>% as.matrix()

# for (i in 1:nrow(plot_dat)) {
# 
#   tmp = plot_dat[i,] %>% scale()
#   plot_dat[i,] = tmp
# 
# }
library(pheatmap)
group = group[order(group$title),]
# plot_dat = plot_dat[DE.cibersort$ImmuneCell,group$sample]
p2 <- pheatmap(plot_dat,  
               color = colorRampPalette(c('blue','white','red'))(100), 
               border_color = NA,  
               scale = "row", 
               cluster_rows = T, 
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
               annotation_names_col = TRUE  #是否显示列注释的名称
)
p2
ggsave(filename = '08.imm//04.heatmap.pdf',p2,w=8,h=2)
ggsave(filename = '08.imm//04.heatmap.png',p2,w=8,h=2)

# cor

library(psych)

results_t = t(results)
dat_t = t(dat)
dat_t = dat_t[rownames(results_t),]

hubgene = c( 'SPSB4',
             'STAG3L3',
             'GUSBP9',
             'PMS2P5',
             'BMS1P6' )
dat_t = dat_t[,hubgene]
corMat = cor( results_t , dat_t )
corMat = corr.test( results_t , dat_t )

r = corMat$r
p = corMat$p

plot_data = reshape2::melt(r)
tmp = reshape2::melt(p)
plot_data$p = tmp$value
colnames(plot_data) = c('cell','gene','cor','p')

ggplot(plot_data,aes(x = gene,y = cell,color = cor,size = -log10(p)))+
  geom_point()+
  scale_color_gradient2(low = 'blue',mid = 'white',high = 'red')+
  scale_y_discrete(limit = r[,1] %>% sort() %>% names())+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 90,vjust = 0.5))




d <- corr.test(x,y,use="complete",method = 'pearson')




