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
celltypes <- c("AT1","AT2","Cil","Clu_Gob", "EC","Stroma",'Mast',"Classical_MC","Nonclassical_MC","Mac1","Mac2","Neu","mDC","pDC",'CD4+TC',"CD8+TC",'NK',"NKT","BC","Pla",'MKI67+_pro')

lung@meta.data$celltype<-lung@active.ident

A_C <- data.frame()
F_C <- data.frame()
F_A <- data.frame()

Idents(lung) <- paste(lung$celltype, lung$group, sep='_')
for (cell in celltypes){
  tmp <- FindMarkers(lung, ident.1=paste0(cell, '_A_TB'), ident.2=paste0(cell, '_Ctrl'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'A_TB-Ctrl'
  A_C <- rbind(A_C, tmp)
  print(paste0(cell, ' is finished'))
}
write.csv(A_C,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","A_TB-Ctrl.csv",sep=""))
A_C.deg <- subset(A_C, p_val_adj<0.05 & abs(avg_log2FC)>0.5)
write.csv(A_C.deg,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","A_TB-Ctrl_deg.csv",sep=""))

for (cell in celltypes){
  tmp <- FindMarkers(lung, ident.1=paste0(cell, '_F_TB'), ident.2=paste0(cell, '_Ctrl'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'F_TB-Ctrl'
  F_C <- rbind(F_C, tmp)
  print(paste0(cell, ' is finished'))
}
write.csv(F_C,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","F_TB-Ctrl.csv",sep=""))
F_C.deg <- subset(F_C, p_val_adj<0.05 & abs(avg_log2FC)>0.5)
write.csv(F_C.deg,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","F_TB-Ctrl_deg.csv",sep=""))

for (cell in celltypes){
  tmp <- FindMarkers(lung, ident.1=paste0(cell, '_F_TB'), ident.2=paste0(cell, '_A_TB'))
  tmp$gene <- rownames(tmp)
  tmp$celltype <- cell
  tmp$compare <- 'F_TB-A_TB'
  F_A <- rbind(F_A, tmp)
  print(paste0(cell, ' is finished'))
}
write.csv(F_A,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","F_TB-A_TB.csv",sep=""))
F_A.deg <- subset(F_A, p_val_adj<0.05 & abs(avg_log2FC)>0.5)
write.csv(F_A.deg,file=paste("/05_Result/07_co_lung/01_ref/05_DEGs/","F_TB-A_TB_deg.csv",sep=""))

