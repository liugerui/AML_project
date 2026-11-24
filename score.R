library(GEOquery)
library(Biobase)
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
library(Seurat)
library(ggplot2)
library(SingleR)
library(harmony)


scobj = readRDS('02.seurat/scobj.rds')
IS = read.csv('01.rawdata/GeneCards-Immune Escape.csv',row.names = 1)
ISRGS = rownames(IS)[1:50]

dat = readRDS('01.rawdata/expr.rds')

IMRES = readxl::read_excel('E:/李老师骨关节炎/02.DEGs/Table_2.XLSX')
IMRES = IMRES$`Table S2: Different immune factors were downloaded from the TISIDB database.` %>% unique()

# imm escape -----------------------------------------
library(GSVA)
## sample gene set sizes
gs <- list(FE = c("ACSL4", "LPCAT3", "ALOX15", "ALOX12", "POR", "TFRC", "IREB2", "CHAC1", "SAT1", "PEBP1", "ALOX5"),
           IFE = c("GPX4", "SLC7A11", "NFE2L2", "FSP1", "ACSL3", "FTH1"))
gsvaPar <- gsvaParam(dat %>% as.matrix(), gs)
gsva.es <- gsva(gsvaPar, verbose=FALSE) 
gsva.es = as.data.frame(gsva.es) %>% as.matrix()

scobj$cellnames = colnames(scobj)
scobj$FE = gsva.es['FE',scobj$cellnames] %>% scale()
scobj$IFE = gsva.es['IFE',scobj$cellnames] %>% scale()
scobj$title = gsva.es['title',scobj$cellnames] %>% scale()

p = FeaturePlot(scobj,features = c('FE'))+
  scale_color_gradient2(low = 'blue',mid = 'white',high = 'red')+
  coord_fixed(ratio = 10/7.5)
p

ip = FeaturePlot(scobj,features = c('IFE'))+
  scale_color_gradient2(low = 'blue',mid = 'white',high = 'red')+
  coord_fixed(ratio = 10/7.5)
ip

ggsave('04.Score/FE.pdf',p,width =5.5,height = 5.5)
ggsave('04.Score/IFE.pdf',ip,width =5.5,height = 5.5)



plotdata = data.frame(cell = scobj$cellnames,ISRGs = scobj$ISRGs,cluster = scobj$seurat_clusters)
p = ggplot(plotdata,aes(x = cluster,y = ISRGs,fill = cluster))+
  geom_boxplot()+
  ggsci::scale_fill_aaas()+
  theme_bw()+
  theme(panel.border = element_blank(),
        axis.text = element_text(colour = 'black'),legend.position = '')
p  
ggsave('04.Score/boxplot.pdf',p,width = 8,height = 4.5)
  
  

all.marker = FindAllMarkers(scobj,group.by = 'seurat_clusters',logfc.threshold = 0.25,min.pct = 0.1,only.pos = T)
all.marker %>%
  group_by(cluster) %>%
  mutate(
    specificity_score = (avg_log2FC * 2) 
  ) %>%
  top_n(n = 5, wt = specificity_score) -> alltop_markers

pdf(file = '02.seurat/markerGene.pdf',width = 6,height = 4.5)
DoHeatmap(scobj,group.by = 'seurat_clusters',features = alltop_markers$gene,size = 3,group.colors = ggsci::pal_aaas()(10) )
dev.off()

p = DimPlot(scobj,pt.size = 3)+
  ggsci::scale_color_aaas()+
  coord_fixed(ratio = 10/7.5)
  
ggsave('02.seurat/Dimplot.pdf',p,width = 5.5,height = 5.2)

intersect( ISRGS , all.marker$gene[which(all.marker$cluster %in% c(0,5))])




# marker enrichment ---------------------------
library(clusterProfiler)

all.marker = all.marker[which(all.marker$p_val < 0.0001),]
all.marker$cluster %>% table

gene = all.marker$gene %>% as.character()
entrezIDs = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= "org.Hs.eg.db", drop = TRUE)
#> Warning in bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb =
#> "org.Hs.eg.db", : 7.1% of input gene IDs are fail to map...
gene<- entrezIDs$ENTREZID

marker_new = data.frame(
  gene = all.marker$gene,
  ID = 0,
  cluster = all.marker$cluster,
  fc = all.marker$avg_log2FC
)
marker_new$ID = entrezIDs$ENTREZID[match(marker_new$gene,entrezIDs$SYMBOL)]
marker_new = marker_new[which(marker_new$fc > 0.5),]

gcSample=split(marker_new$ID, marker_new$cluster) 

xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Hs.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.5,
                     qvalueCutoff = 0.5
)
tab = xx@compareClusterResult
# tab$Cluster = factor(tab$Cluster,levels = c('MOL-2','MOL-1','MFOL','NFOL'))
tab = tab[which(tab$pvalue<0.05),]
xx@compareClusterResult = tab

p <- dotplot(xx,color = "pvalue")
# library(ggthemes)
p1 <- p + theme(axis.text.x = element_text(
  angle = 90,hjust = 1,vjust = 0.5
),legend.position = 'top') + 
  coord_flip()

p1


# stk39 -----------------------------------

tab = data.frame(
  
  ISRGs = scobj$ISRGs,
  cell = colnames(scobj),
  STK39 = scobj@assays$RNA$data['STK39',colnames(scobj)] %>% as.numeric()
  
)

p6<-ggscatter(tab,x="STK39",y="ISRGs",
              add = "reg.line",
              conf.int = T)+
  labs(x = "STK39", 
       y = "ISRGs")+
  stat_cor(label.x.npc = 0.2, label.y.npc = 0.1)
p6
ggsave('oldresult/cor.pdf',p6,width = 4.5,height = 4.5)














