library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(SeuratData)
library(patchwork)
library(Rserve)
set.seed(2)
lung <- readRDS(file = "/07_co_lung/01_ref/04_celltype/Combination_celltype.rds")
lung <- subset(lung, celltype %in% c("AT1","AT2","Ciliated","Clu_Club","EC","Stroma"))

int.list <- SplitObject(lung, split.by = "sample")

# normalize and identify variable features for each dataset independently
int.list <- lapply(X = int.list, FUN = function(x) {
    DefaultAssay(x) <- "RNA"
	x <- SCTransform(x, verbose = FALSE)
})

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 20000 * 1024^2)

int.features <- SelectIntegrationFeatures(object.list = int.list, nfeatures = 3000)
int.list <- PrepSCTIntegration(object.list = int.list, anchor.features = int.features, verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT",
                                      anchor.features = int.features, verbose = FALSE)

plan("multiprocess", workers = 1)
options(future.globals.maxSize = 20000 * 1024^2)

lung.new <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)
DefaultAssay(lung.new) <- "integrated"

lung.new <- RunPCA(lung.new, verbose = FALSE)

pdf('/07_co_lung/01_ref/04_celltype/Stroma/pca_heatmap.pdf', height=20, width=10)
DimHeatmap(lung.new, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

pdf('/07_co_lung/01_ref/04_celltype/Stroma/pca_elbow.pdf', height=5, width=8)
p=ElbowPlot(lung.new, ndims=50)
print(p)
dev.off()
saveRDS(lung.new, file = '/07_co_lung/01_ref/04_celltype/Stroma/beforePC.rds')

pc.num=c(1:11,13)
lung <- readRDS(file = '/07_co_lung/01_ref/04_celltype/Stroma/beforePC.rds')

lung <- RunUMAP(lung, dims = pc.num)# ,min.dist=0.2,spread=3)
lung <- FindNeighbors(lung, reduction = "pca", dims = pc.num)
lung <- FindClusters(lung, resolution = 2.0)#2.2

saveRDS(lung, file = '/07_co_lung/01_ref/04_celltype/Stroma/umap.rds')

plot <- DimPlot(lung, reduction = "umap", label = TRUE,pt.size = 0.3,raster = FALSE)
ggsave("/07_co_lung/01_ref/04_celltype/Stroma/umap.pdf", plot = plot, width = 8, height = 6)

pdf('/07_co_lung/01_ref/04_celltype/Stroma/umap_sample_single.pdf',width=24,height=18)
DimPlot(lung, reduction = "umap", split.by = "sample",pt.size = 0.3,label.size = 6,ncol=3,raster = FALSE)
dev.off()

pdf('/07_co_lung/01_ref/04_celltype/Stroma/UAMP_group.pdf',width=8,height=6)
DimPlot(lung, reduction = "umap", pt.size = 0.3, group.by = "group", cols= c('#a0c3d5','#c9caca', '#4b2415'),raster = FALSE)
DimPlot(lung, reduction = "umap", pt.size = 0.3, group.by = "group", cols= c('#a0c3d5',"NA", "NA"),raster = FALSE)
DimPlot(lung, reduction = "umap", pt.size = 0.3, group.by = "group", cols= c("NA",'#c9caca', "NA"),raster = FALSE)
DimPlot(lung, reduction = "umap", pt.size = 0.3, group.by = "group", cols= c("NA","NA", '#4b2415'),raster = FALSE)
dev.off()

DefaultAssay(lung) <- "RNA"
lung <- NormalizeData(object = lung, normalization.method = "LogNormalize")

celltypes <- c('AT1', 'AT2','Cil', 'Clu', 'Gob', 'EC', 'Fib','SMC')
C1 <- c(10,15,24,28)
C2 <- c(1,9,12,29)
C3 <- c(3,4,6,14,19,20,22)
C4 <- c(8,18,23,26)
C5 <- c(5,25)
C6 <- c(0,2,11,13,16,30)
C7 <- c(7,17,21)
C8 <- c(27)

for (i in 1:8) {
  m <- paste0('C',i)
  tmp <- get(m)
  assign(celltypes[i],tmp)
}

new.cluster.ids <- c(rep(NA,length(unique(Idents(lung)) )))
cluster <- function(celltypes,new.cluster.ids){
  for (cell in celltypes) {
    tmp <- get(cell)
    for (i in tmp) {
      new.cluster.ids[i+1] <- cell
    }
  }
  return(new.cluster.ids)

}

new.cluster.ids <- cluster(celltypes,new.cluster.ids)

lung <- RenameIdents(lung,c("0"=new.cluster.ids[1], "1"=new.cluster.ids[2], "2"=new.cluster.ids[3], "3"=new.cluster.ids[4],"4"= new.cluster.ids[5], "5"=new.cluster.ids[6],"6"=new.cluster.ids[7], "7"=new.cluster.ids[8],"8"= new.cluster.ids[9],"9"= new.cluster.ids[10],
                               "10"= new.cluster.ids[11], "11"=new.cluster.ids[12], "12"=new.cluster.ids[13], "13"=new.cluster.ids[14], "14"=new.cluster.ids[15], "15"=new.cluster.ids[16],"16"= new.cluster.ids[17],"17"= new.cluster.ids[18],"18"= new.cluster.ids[19],"19"=new.cluster.ids[20],
                               "20"=new.cluster.ids[21], "21"=new.cluster.ids[22], "22"=new.cluster.ids[23], "23"=new.cluster.ids[24],"24"= new.cluster.ids[25], "25"=new.cluster.ids[26],"26"=new.cluster.ids[27],"27"=new.cluster.ids[28],"28"= new.cluster.ids[29],"29"= new.cluster.ids[30],"30"= new.cluster.ids[31]))

Idents(lung) <- factor(Idents(lung), levels=celltypes)
lung$celltype <- Idents(lung)
CT.col <- c('#c72e29','#d5695d','#E5D0B0','#fb832d','#DB9078','#65a479',"#F0E68C","#C09736")

UMAP <- DimPlot(lung, reduction = "umap", cols=CT.col, label = TRUE, pt.size = 0.1,label.size = 3.0,raster = FALSE)
ggsave("/01_result/07_lung/00_scRNA_new/04_celltype/Stroma/celltype_UMAP.pdf", plot = UMAP, width = 6, height = 5)
UMAP <- DimPlot(lung, reduction = "umap", cols=CT.col, label = FALSE, pt.size = 0.1,raster = FALSE)
ggsave("/01_result/07_lung/00_scRNA_new/04_celltype/Stroma/celltype_UMAP2.pdf", plot = UMAP, width = 6, height = 5)

#markers <- c('CLIC5', 'AGER', 'CAV1','SFTPB', 'SFTPC', 'ABCA3', 'SCGB3A2','MUC4','CLIC6',"MUC5B", 'MUC5AC', "SPDEF","TPPP3",'FOXJ1','TP73', 'CCDC78',"CPE","GJA5","BMX", "DKK2","CA4",'EDNRB','PROX1','PDPN','DCN','LUM','PDGFRA','COL1A1','ABCA10',"ACTA2", "TAGLN")
markers <- c('CLIC5', 'AGER', 'CAV1','SFTPB', 'SFTPC', 'ABCA3', 'SCGB3A2','MUC4','CLIC6',"MUC5B", 'MUC5AC', "SPDEF","TPPP3",'FOXJ1','TP73', 'CCDC78',"GNG11", 'CLDN5','DCN','LUM','COL1A1','PDGFRA',"ACTA2",'MYH11')
Idents(lung) <- factor(lung$celltype, levels=rev(celltypes))
DotPlot<-DotPlot(lung,features= markers, cols=c('grey90', '#C63C3C'))+RotatedAxis()
ggsave("/01_result/07_lung/00_scRNA_new/04_celltype/Stroma/DotPlot.pdf", plot = DotPlot, width = 8, height = 4)
saveRDS(lung, file = "/01_result/07_lung/00_scRNA_new/04_celltype/Stroma/Combination_celltype.rds")
