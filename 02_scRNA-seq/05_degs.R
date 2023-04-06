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
lung <- readRDS(file = "/05_Result/07_co_lung/01_ref/04_celltype/Combination_celltype.rds")
celltypes <- c("AT1","AT2","Ciliated","Clu_Gob", "EC","Stroma",'Mast',"Classical_MC","Nonclassical_MC","Mac1","Mac2","Neu","mDC","pDC",'CD4+TC',"CD8+TC",'NK',"NKT","BC","Pla",'MKI67+_pro')

lung@meta.data$celltype<-lung@active.ident

A_C <- data.frame()
L_C <- data.frame()
L_A <- data.frame()

Idents(lung) <- paste(lung$celltype, lung$group, sep='_')
for (cell in celltypes){
  tmp <- FindMarkers(lung, ident.1=paste0(cell, '_AL_TB'), ident.2=paste0(cell, '_Ctrl'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'AL_TB-Ctrl'
  A_C <- rbind(A_C, tmp)
  print(paste0(cell, ' is finished'))
}
write.csv(A_C,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","AL_TB-Ctrl.csv",sep=""))
A_C.deg <- subset(A_C, p_val_adj<0.05 & abs(avg_log2FC)>0.5)
write.csv(A_C.deg,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","AL_TB-Ctrl_deg.csv",sep=""))

for (cell in celltypes){
  tmp <- FindMarkers(lung, ident.1=paste0(cell, '_L_TB'), ident.2=paste0(cell, '_Ctrl'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'L_TB-Ctrl'
  L_C <- rbind(L_C, tmp)
  print(paste0(cell, ' is finished'))
}
write.csv(L_C,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","L_TB-Ctrl.csv",sep=""))
L_C.deg <- subset(L_C, p_val_adj<0.05 & abs(avg_log2FC)>0.5)
write.csv(L_C.deg,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","L_TB-Ctrl_deg.csv",sep=""))

for (cell in celltypes){
  tmp <- FindMarkers(lung, ident.1=paste0(cell, '_L_TB'), ident.2=paste0(cell, '_AL_TB'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'L_TB-AL_TB'
  L_A <- rbind(L_A, tmp)
  print(paste0(cell, ' is finished'))
}
write.csv(L_A,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","L_TB-AL_TB.csv",sep=""))
L_A.deg <- subset(L_A, p_val_adj<0.05 & abs(avg_log2FC)>0.5)
write.csv(L_A.deg,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","L_TB-AL_TB_deg.csv",sep=""))

