library(GSVA)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(stringr)

library(WGCNA)
library(GSVA)
library(readr)
library(dplyr)
library(survival)
library(clusterProfiler)
library(ggpubr)



# GeneSET = read.gmt('03.WGCNA/GOBP_LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS.v2024.1.Hs.gmt')
scobj = readRDS('02.seurat/scobj.rds')

expr = readRDS('01.rawdata/expr.rds')
expr = na.omit(expr)
# group = readRDS('01.rawdata/heart_group.rds')

Res = readRDS('06.DEGs/01.limma_out.rds')
DEGs = Res$gene[which(Res$adj.P.Val<0.05)]

datExpr = expr[DEGs,] %>% t()
usedGene = c()
for(i in 1:ncol(datExpr)){
  
  if(sum(datExpr[,i]) > 12){
    usedGene = c(usedGene,i)
  }
  
}
datExpr = datExpr[,usedGene]

SampleName <- colnames(expr)
# 一步构建共表达网络和划分模块 
##设置网络构建参数选择范围，计算beta值（为了计算power）
powers <- c(seq(1, 10, by=1),seq(12, 20, by=2)) #设置beta值的取值范围
rsquare.cut=0.85  #设置r2阈值为0.85，一般选择0.85或0.9

#开启多线程，加快运行速度
enableWGCNAThreads()
#disbleWGCNAThreads()     #如果运行pickSoftThreshold报call()相关的错误，需要运行禁用多线程
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
sft$powerEstimate  # 6

#软阈值选择绘图
#左图为每个候选beta对应的R2值,
#绘制阈值线0.9（红线）和0.85（蓝线），如没有超过蓝线，则没找到软阈值，sft$powerEstimate为空
#右图为每个候选beta对应的平均连接度

pdf('07.wgcna/03.softThreshold.pdf',w=12,h=8)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1)) 
cex1=1.5
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2,
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"),
     cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

power <- 4

# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，
# 数值越高，网络越符合无标度特征 (non-scale)。
##自动构建WGCNA模型
#表达矩阵转换成邻接矩阵，然后再将邻接矩阵转换成拓扑矩阵,识别模块
cor <- WGCNA::cor
net = blockwiseModules(
  datExpr,
  power = power,
  minModuleSize =100,                    #每个模块最少的基因数
  deepSplit = 4,                         #剪切树参数，一般设置为2，取值0-4（最灵敏）
  mergeCutHeight = 0.1,                 #模块合并参数，越大模块越少，默认0.25  越小，模块越多
  numericLabels = TRUE,                  # T返回数字，F返回颜色
  networkType  = "signed",               #网络类型，一般不需要修改
  maxBlockSize = ncol(datExpr),     
  pamRespectsDendro = FALSE, 
  saveTOMs = TRUE,
  saveTOMFileBase = "FPKM-TOM",
  loadTOMs = TRUE,
  verbose = 3,
)
cor<-stats::cor
table(net$colors)    
saveRDS(net,file = '07.wgcna/net.rds')


##模块重新排序
#获取模块颜色
moduleColors = labels2colors(net$colors)
#获取每个模块的特征值
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); #重新排序，相似颜色的挨在一起
useMEs=subset(MEs, select = -c(MEgrey))  #去掉未分组的基因(灰色模块）
table(moduleColors)

write.csv(table(moduleColors),'07.wgcna//04module_number.csv',quote = F,row.names = F)
write.csv(data.frame(group = moduleColors,gene = names(net$colors) ),'07.wgcna/gene2mod.csv',quote = F,row.names = F)

#useMEs=subset(MEs)  #去掉未分组的基因(灰色模块） 
## 输出每个基因所在的模块，以及与该模块的KME值
if(file.exists('05All_Gene_KME.txt')) file.remove('05All_Gene_KME.txt')
for(module in substring(colnames(useMEs),3)){
  ME=as.data.frame(useMEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=datExpr[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote = F,row.names = T,append = T,file = "07.wgcna//05All_Gene_KME.txt",col.names = F)
}

##模块绘图
mergedColors = labels2colors(net$colors)
mergedColors[net$blockGenes[[1]]]
pdf("07.wgcna//06wgcna.dendroColors.pdf",height = 5,width = 9)
par(mfrow = c(1,2),pin = c(4,4), mar = c(6,6,6,1)) 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,
                    cex.axis=1.6,  ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
                    cex.main=1.8,   ##标题的缩放倍数。
                    font.axis = 2, 
                    font.lab = 2, 
                    font.main = 2, 
                    font.sub =2)

dev.off()


# 模块与表型的关联分析
library(ggpubr)

plot_data = cbind(useMEs,scobj$title[rownames(useMEs)])
plot_data = cor(useMEs,scobj$title[rownames(useMEs)])

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tibble)

# 假设你的数据框名为df
# 首先将数据转换为适合ggplot的格式
df <- data.frame(
  module = rownames(plot_data),
  value = plot_data[,1]
)

# 创建棒棒糖图
p = ggplot(df, aes(x = reorder(module, value), y = value)) +
  geom_segment(aes(xend = module, yend = 0), color = "gray50") +
  geom_point(aes(color = value), size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "Value") +
  coord_flip() +
  labs(
    title = "Module Eigengene Values",
    x = "Module",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 10,color = 'black'),
    legend.position = "right"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30")+
  ylab('Correlation')+xlab('')
p
ggsave('07.wgcna/cor.pdf',p,width = 4.5,height = 4.5)


# 功能富集 ----------------------------------------

gene2mod = read.table('07.wgcna/05All_Gene_KME.txt')
gene2mod = gene2mod[which(gene2mod$V1 != '--'),]
colnames(gene2mod) = c('gene','cor','group')

0
gene = gene2mod$gene
entrezIDs = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= "org.Hs.eg.db", drop = TRUE)
gene2mod$ID = entrezIDs$ENTREZID[match(gene2mod$gene,entrezIDs$SYMBOL)]
gcSample=split(gene2mod$ID, gene2mod$group) 

xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.5,
                     qvalueCutoff = 0.5
)

tab = xx@compareClusterResult
tab = tab[tab$pvalue < 0.05,]

genesTab = table(tab$Description %>% as.character())
# genesTab = genesTab[which(genesTab < 4)]
genesTab = names(genesTab)

tab = tab[which(tab$Description %in% genesTab),]
xx@compareClusterResult = tab

p <- enrichplot::dotplot(xx,showCategory = 2,color = "pvalue")
p1 <- p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
),legend.position = 'top')
p1

ggsave('07.wgcna/enrich.pdf',p1,width = 10,height = 10)



# -----------------------

plist = list()
tablist = list()
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
geneModuleMembership = geneModuleMembership[which(!rownames(geneModuleMembership) %>% str_detect('X...')),]

DEGs = Res$gene[which(Res$regulate != 'Normal')]

# TF= read.table('03.WGCNA/TF-Target-information.txt',sep = '\t',header = T)

# saveRDS(geneModuleMembership,file = '03.WGCNA/geneModuleMembership.rds')
for (i in 1:ncol(geneModuleMembership)) {
  
  genes = colnames(geneModuleMembership)[i]
  genes = gene2mod$gene[gene2mod$group == gsub(genes,pattern = 'ME',replacement = '')]
  
  tab = geneModuleMembership[,i]
  names(tab) = rownames(geneModuleMembership)
  tab = tab %>% sort(decreasing = T)
  tab = data.frame(tab)
  tab$gene = rownames(tab)
  tab = tab[which(toupper(tab$gene) %in% DEGs),]
  tab = tab[1:10,]
  tab = tab[which(tab$gene %in% genes),]
  
  
  tablist[[i]] = tab
  
  p = ggplot(tab, aes(x = gene, y = tab)) +
    geom_segment(aes(xend = gene), yend = min(tab$tab),
                 linewidth = 1, color = "gray") +  # “棍”：线段从y=0到数据点
    geom_point(size = 4, color =colnames(geneModuleMembership)[i] %>% gsub(pattern = 'ME',replacement = '')) +          # “糖”：数据点
    labs(title = colnames(geneModuleMembership)[i],
         x = "",
         y = "Module Membership (MM)") +
    theme_minimal() +
    theme(axis.text = element_text(color = 'black'),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_x_discrete(limit = tab$gene)
  p
  
  plist[[i]] = p
  
}


p = plist[[1]] + plist[[2]] + plist[[3]] +
  plist[[4]] + plist[[5]] + plist[[6]]+
  plist[[7]] + plist[[8]] + plist[[9]]+
  plist[[10]] + 
  patchwork::plot_layout(ncol = 2)

p
