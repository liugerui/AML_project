
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

# anno = readRDS('E:/shengao/R/Genecode_anno.Hs.Rds')

expr <- readRDS('01.rawdata/expr.rds')
expr = expr[intersect(rownames(expr),anno$symbol),]
dim(expr)

scobj <- CreateSeuratObject(counts = expr, 
                            project = 'NPC', 
                            min.cells = 0, 
                            min.features = 200)

scobj$cellname = colnames(scobj)
scobj <- SetAssayData(scobj, 
                          assay = "RNA", 
                          layer = "data", 
                          new.data = log2(1 + scobj@assays$RNA$counts)) 
  
  
scobj <- FindVariableFeatures(scobj,
                              selection.method = "vst",
                              nfeatures = 4000)

all.genes <- rownames(scobj)
scobj <- ScaleData(scobj, features = all.genes)
scobj <- RunPCA(scobj,features = VariableFeatures(object = scobj))
scobj <- FindNeighbors(scobj, reduction = 'pca',dims = 1:15) # 前10个维度
scobj <- FindClusters(scobj, reduction = 'pca',resolution = 1.3)
scobj <- RunUMAP(scobj,dims = 1:15, reduction = 'pca')
scobj <- RunTSNE(scobj,dims = 1:15, reduction = 'pca')

DimPlot(scobj,label = T,reduction = 'tsne')
DimPlot(scobj,label = T,label.color = 'white',label.size = 5,reduction = 'umap',
        pt.size = 4,
        group.by = 'seurat_clusters')+
  ggsci::scale_color_aaas()+
  coord_fixed()+
  xlab('UMAP1')+ylab('UMAP2')

# all.marker = FindAllMarkers(scobj,group.by = 'seurat_clusters',logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
# 
# write.table(all.marker[,c('avg_log2FC','p_val','gene','cluster')],file = '02.seurat/allmarker.txt',col.names = T,row.names = F,sep = '\t',quote = F)
# 
# FeaturePlot(scobj,features = all.marker$gene[which(all.marker$cluster == '1')])

saveRDS(scobj,file = '02.seurat/scobj.rds')


# TCGA ---------------------------------------------------------












