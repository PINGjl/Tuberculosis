library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(DoubletFinder)
library(future)
library(RColorBrewer)
library(SeuratData)
library(patchwork)
set.seed(2)
# Read in data ####
exercise_dir <- '/01_rawdata/01_lung/01_ref/'
read_in_sample <- function(i){
  tmp <- Read10X(paste0(exercise_dir, i, '/outs/filtered_feature_bc_matrix/'))
  tmp <- CreateSeuratObject(counts = tmp, min.features = 200, min.cells = 3, project = i)
  tmp$sample <- i
  tmp$group<-substr(i,1,4)
  tmp[["percent.MT"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  ggsave(paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_vln.pdf"), plot = p, width = 10, height = 6)
  plot1 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "percent.MT")
  plot2 <- FeatureScatter(tmp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- plot1 + plot2
  ggsave(paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_Scatter.pdf"),plot = plot3, width = 8, height = 4)
  saveRDS(tmp, file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_before_qc.rds"))
  
  tmp <- subset(tmp, subset = percent.MT < 10 & nFeature_RNA > 200 )#& nFeature_RNA < 7000
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  ggsave(paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_qc_vln.pdf"), plot = p, width = 10, height = 6)
  saveRDS(tmp, file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_after_qc.rds"))
  tmp <- RenameCells(tmp, add.cell.id = i)
}

sample_normalize <- function(tmp, i){
  print(date())
  print(paste0(i, ': SCTransform started'))
  tmp <- SCTransform(tmp, verbose = FALSE)
  print(date())
  print(paste0(i, ': SCTransform finished'))

  saveRDS(tmp, file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_SCT.rds"))

  tmp <- RunPCA(tmp, verbose=F)

  pdf(paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_pca_heatmap.pdf"), width=10,height=30)
  DimHeatmap(tmp, dims=1:50, cells=500, balanced=T)
  dev.off()

  p<- ElbowPlot(tmp, ndims = 50)
  pdf(paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_ElbowPlot.pdf"), height = 6, width = 8)
  print(p)
  dev.off()
  
  saveRDS(tmp, file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_bfPCR.rds"))

  return(tmp)
  print(paste0(i ,' completed'))
}

for (sample in c('Ctrl-1', 'Ctrl-2', 'Ctrl-3', 'Ctrl-4', 'Ctrl-5', 'Ctrl-6', 'Ctrl-7', 'AL_TB-1', 'AL_TB-2', 'AL_TB-3', 'L_TB-1', 'L_TB-2', 'L_TB-3')){
  tmp <- read_in_sample(sample)
  assign(sample, tmp)
}

for (sample in c('Ctrl-1', 'Ctrl-2', 'Ctrl-3', 'Ctrl-4', 'Ctrl-5', 'Ctrl-6', 'Ctrl-7', 'AL_TB-1', 'AL_TB-2', 'AL_TB-3', 'L_TB-1', 'L_TB-2', 'L_TB-3')){
  tmp <- sample_normalize(get(sample), sample)
  assign(sample, tmp)
}

#SCT
dims <- list(c(1:30), c(1:30), c(1:30),c(1:30),c(1:30), c(1:30), c(1:30),c(1:30),c(1:30), c(1:30), c(1:30),c(1:30),c(1:30))
samples <- c('Ctrl-1', 'Ctrl-2', 'Ctrl-3', 'Ctrl-4', 'Ctrl-5', 'Ctrl-6', 'Ctrl-7', 'AL_TB-1', 'AL_TB-2', 'AL_TB-3', 'L_TB-1', 'L_TB-2', 'L_TB-3')
res <- c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8)

for (i in seq(1:13)){
  tmp <- readRDS(file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",samples[i],"/",samples[i],"_bfPCR.rds"))
  tmp <- RunUMAP(tmp, dims = dims[[i]], verbose=F)
  tmp <- FindNeighbors(tmp, reduction = "pca", dims = dims[[i]])
  tmp <- FindClusters(tmp, res=res[i])
  tmp[["cluster"]] <- Idents(tmp)
  UMAP <- DimPlot(object = tmp, reduction = "umap", label = TRUE)
  ggsave(paste0("/05_Result/07_co_lung/01_ref/01_qc/",samples[i],"/",samples[i],"_umap.pdf"), plot = UMAP, width = 8, height = 6)
  saveRDS(tmp, file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",samples[i],"/",samples[i],"_PCR.rds"))
  assign(samples[i], tmp)
  print(paste0(samples[1], ' completed'))
}

# Identify doublets ####
pK.df <- data.frame(matrix(nrow=0, ncol=3))
colnames(pK.df) <- c("Sample", "Optimal_pK","After_qc")

for (i in c(1:13)){
  sweep.res.list <- paramSweep_v3(get(samples[i]), PCs = dims[[i]], sct = T, num.cores=6)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- arrange(bcmvn, desc(BCmetric))$pK[1]
  tmp <- data.frame(Sample=samples[i], Optimal_pK=pK,After_qc=length(get(samples[i])@meta.data$orig.ident))
  pK.df <- rbind(pK.df, tmp)
  print(bcmvn)
  print(paste0("--------------", samples[i], " completed (", i, "/13)--------------"))
}
write.table(pK.df,file ="/05_Result/07_co_lung/01_ref/01_qc/Optimal_pK.txt", sep = "\t" )

doublet.prop <- data.frame(matrix(nrow=0, ncol=3))
colnames(doublet.prop) <- c("Sample", "Number","Doublet_prop")

dims <- list(c(1:30), c(1:30), c(1:30),c(1:30),c(1:30), c(1:30), c(1:30),c(1:30),c(1:30), c(1:30), c(1:30),c(1:30),c(1:30))
ratio <- c(0.092,0.092,0.092,0.108,0.076,0.092,0.084,0.124,0.084,0.069,0.108,0.084,0.108)

for (i in seq(1:13)) {
  tmp <- readRDS(file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",samples[i],"/",samples[i],"_PCR.rds"))
  pK.use <- as.numeric(as.character(pK.df$Optimal_pK[i]))
  homotypic.prop <- modelHomotypic(tmp@meta.data$cluster)
  nExp_poi <- round(ratio[i]*length(tmp@meta.data$orig.ident))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = dims[[i]], pN = 0.25, pK = pK.use, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  tmp[["doublet"]] <- tmp[[paste("DF.classifications_0.25", pK.use, nExp_poi.adj, sep="_")]]
  prop <- nExp_poi.adj/length(tmp@meta.data$cluster)
  prop.tmp <- data.frame(Sample=samples[i], Number=nExp_poi.adj, Doublet_prop=prop)
  doublet.prop <- rbind(doublet.prop, prop.tmp)
  saveRDS(tmp, file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",samples[i],"/",samples[i],"_doublets.rds"))
  assign(samples[i], tmp)
  print(paste0("--------------", samples[i], " completed (", i, "/13)--------------"))
}

write.table(doublet.prop,file ="/05_Result/07_co_lung/01_ref/01_qc/doublet.prop.txt", sep = "\t" )

for (i in samples){
  doub <- DimPlot(get(i), group.by='doublet', cols=c('firebrick', 'grey90'))
  pdf(paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_UMAP_doublet.pdf"), height = 6, width = 8)
  print(doub)
  dev.off()
}

for (i in samples){
  tmp <- get(i)
  tmp <- subset(tmp, doublet=='Singlet')
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  ggsave(paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_qc_vln1.pdf"), plot = p, width = 10, height = 6)
  tmp <- subset(tmp, subset = nFeature_RNA < 7000 )
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
  ggsave(paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_qc_vln2.pdf"), plot = p, width = 10, height = 6)
  tmp <- SCTransform(tmp, verbose = FALSE)
  assign(i, tmp)
  saveRDS(tmp, file = paste0("/05_Result/07_co_lung/01_ref/01_qc/",i,"/",i,"_final.rds"))
   print(paste0(i, ' completed'))
}
