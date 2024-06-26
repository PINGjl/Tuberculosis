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
lung <- readRDS(file = '/05_Result/07_co_lung/01_ref/02_int/umap.rds')
DefaultAssay(lung) <- "RNA"
lung <- NormalizeData(object = lung, normalization.method = "LogNormalize")

##AT1 AT2 Cil Clu Gob EC Stroma Mast Classical_MC Nonclassical_MC Mac Neu mDC pDC CD4+TC CD8+TC NK NKT BC Pla MKI67+_pro
lung <- subset(lung, idents = c(56,61,62,65), invert = TRUE) ## >70% in one sample

celltypes <- c("AT1","AT2","Cil","Clu_Gob", "EC","Stroma",'Mast',"cMC","nMC","Mac1","Mac2","Neu","DC",'CD4+TC',"CD8+TC",'NK',"NKT","BC","Pla",'MKI67+_pro')
C1 <- c(28)
C2 <- c(24,44,57)
C3 <- c(20,31,50)
C4 <- c(22)
C5 <- c(15,42,60)
C6 <- c(25)

C7 <- c(23,59,63)
C8 <- c(9,36,46)
C9 <- c(26)
C10 <- c(17)
C11 <- c(1,8,29,37,39,41,43,49,64)
C12 <- c(18,47,51)

C13 <- c(45)

C14 <- c(0,3,6,10,11,14,21,35,53)
C15 <- c(2,7,30,48)
C16 <- c(4,5,13,16,27,32)
C17 <- c(12,19,33,58)

C18 <- c(6,40,54,55)
C19 <- c(38)
C20 <- c(34,52)

for (i in 1:20) {
  m <- paste0('C',i)
  tmp <- get(m)
  assign(celltypes[i],tmp)
}

new.cluster.ids <- c(rep(NA,length(unique(Idents(lung))  )))
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

lung <- RenameIdents(lung,c("0"=new.cluster.ids[1], "1"=new.cluster.ids[2], "2"=new.cluster.ids[3], "3"=new.cluster.ids[4],"4"= new.cluster.ids[5], "5"=new.cluster.ids[6],"6"=new.cluster.ids[7], "7"=new.cluster.ids[8],"8"= new.cluster.ids[9],"9"= new.cluster.ids[10],"10"= new.cluster.ids[11],"11"=new.cluster.ids[12], "12"=new.cluster.ids[13], "13"=new.cluster.ids[14],"14"= new.cluster.ids[15],"15"=new.cluster.ids[16],"16"=new.cluster.ids[17], "17"=new.cluster.ids[18],"18"= new.cluster.ids[19],"19"= new.cluster.ids[20],"20"= new.cluster.ids[21],"21"=new.cluster.ids[22], "22"=new.cluster.ids[23], "23"=new.cluster.ids[24],"24"= new.cluster.ids[25],"25"=new.cluster.ids[26],"26"=new.cluster.ids[27], "27"=new.cluster.ids[28],"28"= new.cluster.ids[29],"29"= new.cluster.ids[30],"30"= new.cluster.ids[31],"31"=new.cluster.ids[32], "32"=new.cluster.ids[33], "33"=new.cluster.ids[34],"34"= new.cluster.ids[35],"35"=new.cluster.ids[36],"36"=new.cluster.ids[37], "37"=new.cluster.ids[38],"38"= new.cluster.ids[39],"39"= new.cluster.ids[40],"40"= new.cluster.ids[41],"41"=new.cluster.ids[42], "42"=new.cluster.ids[43], "43"=new.cluster.ids[44],"44"= new.cluster.ids[45], "45"=new.cluster.ids[46],"46"=new.cluster.ids[47], "47"=new.cluster.ids[48],"48"= new.cluster.ids[49],"49"= new.cluster.ids[50],"50"= new.cluster.ids[51],"51"=new.cluster.ids[52],"52"=new.cluster.ids[53], "53"=new.cluster.ids[54],"54"= new.cluster.ids[55], "55"=new.cluster.ids[56],"56"=new.cluster.ids[57], "57"=new.cluster.ids[58],"58"= new.cluster.ids[59],"59"= new.cluster.ids[60],"60"= new.cluster.ids[61],"61"=new.cluster.ids[62],"62"=new.cluster.ids[63], "63"=new.cluster.ids[64],"64"= new.cluster.ids[65], "65"=new.cluster.ids[66]))

Idents(lung) <- factor(Idents(lung), levels=celltypes)
lung$celltype <- Idents(lung)

CT.col <- c('#c72e29','#d5695d','#F39B7FB2',"#E64B35B2",
"#F0E68C",'#ddd954',
'#804588',"#e2a461",'#fb832d','#cc759d',"#ebc3d5","#c0d8b9","#80bbb2",
'#a0c167','#4b9748',"#00A087B2",'#59a77f','#779eb9','#346e9a','#955025')

UMAP <- DimPlot(lung, reduction = "umap", cols=CT.col, label = TRUE, pt.size = 0.1,label.size = 3.0,raster = FALSE)
ggsave("/01_result/07_lung/00_scRNA_new/04_celltype/celltype_UMAP.pdf", plot = UMAP, width = 6, height = 5)
UMAP <- DimPlot(lung, reduction = "umap", cols=CT.col, label = FALSE, pt.size = 0.1,raster = FALSE)
ggsave("/01_result/07_lung/00_scRNA_new/04_celltype/celltype_UMAP2.pdf", plot = UMAP, width = 6, height = 5)

markers <- c('AGER','PDPN', 'CLIC5','SFTPC','ABCA3','SFTPD',"TPPP3",'FOXJ1','CCDC78','SCGB3A2',"MUC5B", 'MUC5AC','CLDN5','GNG11','DCN',"COL1A1", 'PDGFRA','ACTA2','PTPRC','SPI1','MS4A2', 'CPA3', 'TPSAB1','CD14','FCGR3A','CD68','CD86','CD80','MRC1','CD163','MARCO', 'MSR1','CMTM2','IFITM2', 'PTGS2','FCGR3B','CD1C','CLEC9A','LILRB4', "IRF8",'CD3D','CD3E','CD4','CD8A','KLRD1','NKG7', 'CD79A','MS4A1', 'CD19','MZB1','CD27','SLAMF7','MKI67')
Idents(lung) <- factor(lung$celltype, levels=rev(celltypes))
DotPlot<-DotPlot(lung,features= rev(markers), cols=c('grey90', '#C63C3C'))+RotatedAxis()
ggsave("/05_Result/07_co_lung/01_ref/04_celltype/DotPlot.pdf", plot = DotPlot, width = 16, height = 6)
saveRDS(lung, file = "/05_Result/07_co_lung/01_ref/04_celltype/Combination_celltype.rds")

##Analysis############################################################################
lung@meta.data$celltype<-lung@active.ident

write.csv(lung@meta.data,"/05_Result/07_co_lung/01_ref/04_celltype/Combination_al.csv", row.names =FALSE)

# cell proportion ####
samples <- c('Ctrl-1', 'Ctrl-3', 'Ctrl-2', 'Ctrl-4', 'Ctrl-5', 'Ctrl-6', 'Ctrl-7', 'AL_TB-1', 'L_TB-1', 'AL_TB-2', 'L_TB-2', 'AL_TB-3', 'L_TB-3', 'Ctrl-8', 'Ctrl-9', 'Ctrl-10', 'Ctrl-11', 'Ctrl-12', 'Ctrl-13', 'AL_TB-4', 'AL_TB-5', 'AL_TB-6', 'AL_TB-7', 'AL_TB-8', 'AL_TB-9', 'AL_TB-10', 'AL_TB-11', 'AL_TB-12', 'L_TB-4', 'L_TB-5', 'L_TB-6', 'L_TB-7', 'L_TB-8', 'L_TB-9', 'L_TB-10', 'L_TB-11', 'L_TB-12')
lung$sample <- factor(lung$sample, levels=samples)
cell_number=table(Idents(lung),lung$sample)
write.csv(cell_number,"/05_Result/07_co_lung/01_ref/04_celltype/cell_number_all_al.csv")
cell.prop <- data.frame(table(lung$celltype,lung$sample))
colnames(cell.prop) <- c('celltype', 'sample', 'number')
cell.prop <- cell.prop %>% group_by(sample) %>% mutate(total=sum(number))
cell.prop$prop <- cell.prop$number/cell.prop$total

cell.prop$sample <- factor(cell.prop$sample,
                           levels=samples, ordered=TRUE)
write.csv(cell.prop,"/05_Result/07_co_lung/01_ref/04_celltype/cell_prop_all_al.csv")
p=ggplot(cell.prop, aes(celltype, prop, fill=cell.prop$sample)) +
  geom_bar(stat='identity', position = 'dodge', width=0.9) +
  #scale_y_continuous(expand = c(0,0.005)) +
  #scale_fill_manual(values = c( 'dodgerblue1', 'orange2')) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'),
        axis.text = element_text(color='black'),axis.text.x=element_text(angle=45, hjust=1, vjust=1))
ggsave("/05_Result/07_co_lung/01_ref/04_celltype/cell_prop_al.pdf", plot = p, width = 24, height = 6)

