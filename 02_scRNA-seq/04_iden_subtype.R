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
lung <- subset(lung, celltype %in% c("AT1","AT2","Cil","Clu_Club","EC","Stroma"))

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

pc.num=c(1:12)
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

