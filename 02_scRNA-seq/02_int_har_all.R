library(Seurat)
library(cowplot)
library(harmony)
library(ggplot2)
library(future)
set.seed(2)
samples <- c('Ctrl-1', 'Ctrl-3', 'Ctrl-2', 'Ctrl-4', 'Ctrl-5', 'Ctrl-6', 'Ctrl-7', 'AL_TB-1', 'L_TB-1', 'AL_TB-2', 'L_TB-2', 'AL_TB-3', 'L_TB-3', 'Ctrl-8', 'Ctrl-9', 'Ctrl-10', 'Ctrl-11', 'Ctrl-12', 'Ctrl-13', 'AL_TB-4', 'AL_TB-5', 'AL_TB-6', 'AL_TB-7', 'AL_TB-8', 'AL_TB-9', 'AL_TB-10', 'AL_TB-11', 'AL_TB-12', 'L_TB-4', 'L_TB-5', 'L_TB-6', 'L_TB-7', 'L_TB-8', 'L_TB-9', 'L_TB-10', 'L_TB-11', 'L_TB-12')
for (i in c(1:37)) {
 tmp <- readRDS(file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",samples[i],"/",samples[i],"_final.rds"))
 assign(samples[i],tmp)
}
int.list <- list(get(samples[1]),get(samples[2]),get(samples[3]),get(samples[4]),get(samples[5]),get(samples[6]),get(samples[7]),get(samples[8]),get(samples[9]),get(samples[10]),get(samples[11]),get(samples[12]),get(samples[13]),
                get(samples[14]),get(samples[15]),get(samples[16]),get(samples[17]),get(samples[18]),get(samples[19]),get(samples[20]),get(samples[21]),get(samples[22]),get(samples[23]),get(samples[24]),get(samples[25]),get(samples[26]),
                get(samples[27]),get(samples[28]),get(samples[29]),get(samples[30]),get(samples[31]),get(samples[32]),get(samples[33]),get(samples[34]),get(samples[35]),get(samples[36]),get(samples[37]))

int.features <- SelectIntegrationFeatures(object.list = int.list, nfeatures = 3000)

#在运行Harmony之前，创建一个Seurat对象并按照标准PCA进行分析。
lung <- merge(int.list[[1]],y = int.list[2:length(int.list)], add.cell.ids = samples, project = "scRNA") %>%
    #SCTransform(verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE,features =int.features) #R语言中%>%的含义是什么呢，管道函数啦，就是把左件的值发送给右件的表达式，并作为右件表达式函数的第一个参数。 

p1 <- DimPlot(object = lung, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = lung, features = "PC_1", group.by = "sample", pt.size = .1)
p <- p1+p2
ggsave("/05_Result/07_co_lung/01_ref/02_int/before.pdf", plot = p, width = 12, height = 5)

#pdf('/05_Result/07_co_lung/01_ref/02_int/RunHarmony.pdf', height=2.5, width=6)
lung <- lung %>%
RunHarmony("sample", assay.use = "SCT") #Harmony converged after 8 iterations ,plot_convergence = TRUE
#dev.off()

lung@reductions$harmony

harmony_embeddings <- Embeddings(lung, 'harmony')
#harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = lung, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = lung, features = "harmony_1", group.by = "sample", pt.size = .1)
p <- p1+p2
ggsave("/05_Result/07_co_lung/01_ref/02_int/after.pdf", plot = p, width = 12, height = 5)

#lung <- RunPCA(lung, verbose = FALSE)
pdf('/05_Result/07_co_lung/01_ref/02_int/pca_heatmap.pdf', height=20, width=10)
DimHeatmap(lung, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()
pdf('/05_Result/07_co_lung/01_ref/02_int/pca_elbow.pdf', height=5, width=8)
p=ElbowPlot(lung, ndims = 50)
print(p)
dev.off()
saveRDS(lung, file = '/05_Result/07_co_lung/01_ref/02_int/before_pca.rds')

lung <- lung %>%
    RunUMAP(reduction = "harmony", dims = 1:50) %>%
    FindNeighbors(reduction = "harmony", dims = 1:50) %>%
    FindClusters(resolution = 2.0) %>%
    identity()

pdf('/05_Result/07_co_lung/01_ref/02_int/sample_umap.pdf', height=12, width=45)
DimPlot(lung, reduction = "umap", group.by = "sample", pt.size = .1, split.by = 'sample',raster=FALSE,ncol=7)
dev.off()

pc.num=c(1:25)
lung <- readRDS(file = '/05_Result/07_co_lung/01_ref/02_int/before_pca.rds')

lung <- RunUMAP(lung,reduction = "harmony", dims = pc.num)#,spread = 3)
lung <- FindNeighbors(lung,reduction = "harmony", dims = pc.num)
lung <- FindClusters(lung, resolution = 2.0)

plot <- DimPlot(lung, reduction = "umap", label = TRUE,pt.size = 0.1,raster=FALSE)
ggsave("/05_Result/07_co_lung/01_ref/02_int/umap.pdf", plot = plot, width = 8, height = 6)

pdf('/05_Result/07_co_lung/01_ref/02_int/umap_sample_single.pdf',width=24,height=18)
DimPlot(lung, reduction = "umap", split.by = "sample",pt.size = 0.1,label.size = 6,ncol=4,raster=FALSE)
dev.off()

lung$group <- substr(lung$sample, 1,4)

saveRDS(lung, file = '/05_Result/07_co_lung/01_ref/02_int/umap.rds')

pdf('/05_Result/07_co_lung/01_ref/02_int/UAMP_group.pdf',width=8,height=6)
DimPlot(lung, reduction = "umap", pt.size = 0.1, group.by = "group",raster=FALSE, cols= c('#a0c3d5','#c9caca', '#4b2415'))
DimPlot(lung, reduction = "umap", pt.size = 0.1, group.by = "group",raster=FALSE, cols= c('#a0c3d5',"NA","NA"))
DimPlot(lung, reduction = "umap", pt.size = 0.1, group.by = "group",raster=FALSE, cols= c('NA',"#c9caca","NA"))
DimPlot(lung, reduction = "umap", pt.size = 0.1, group.by = "group", raster=FALSE,cols= c("NA","NA",'#4b2415'))
dev.off()

pdf('/05_Result/07_co_lung/01_ref/02_int/UAMP_sample.pdf',width=8,height=6)
DimPlot(lung, reduction = "umap", pt.size = 0.1, group.by = "sample",raster=FALSE)
dev.off()

