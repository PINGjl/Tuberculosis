library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggridges)
library(showtext)

Combination <- readRDS(file = "/05_Result/07_co_lung/01_ref/04_celltype/Combination_celltype.rds")
celltypes <- c("AT1","AT2","Ciliated","Clu_Gob", "EC","Stroma",'Mast',"Classical_MC","Nonclassical_MC","Mac1","Mac2","Neu","mDC","pDC",'CD4+TC',"CD8+TC",'NK',"NKT","BC","Pla",'MKI67+_pro')
gene_list=read.csv("/05_Result/07_co_lung/01_ref/09_gene_score/gene_list.csv")

for (i in c(1:10) ){
  genes <- list(as.character(get(object[i])))
  Combination <- AddModuleScore(Combination, features=genes, seed=15, name="score")

  Combination$celltype <- factor(Combination$celltype,
                             levels=c("AT1","AT2","Ciliated","Clu_Gob", "EC","Stroma",'Mast',"Classical_MC","Nonclassical_MC","Mac1","Mac2","Neu","mDC","pDC",'CD4+TC',"CD8+TC",'NK',"NKT","BC","Pla",'MKI67+_pro'), ordered=TRUE)

  Combination$group <- factor(Combination$group,
                             levels=c("Ctrl","AL_TB", "L_TB"), ordered=TRUE)

  CT.col <- c('#47B1B6', '#CEB7B3', '#E6949A')
  p <- ggplot(Combination[[]],aes(group,score1,colour=group))+
       geom_violin(aes(fill=group),cex=0.5,trim = FALSE)+
	   scale_fill_manual(values=CT.col) +
	   geom_boxplot(width=0.2,cex=0.5,outlier.shape = NA)+
       stat_compare_means(label = "p.format", label.x = 1.5,
                     comparisons=list(c("Ctrl","L_TB"),c("Ctrl","AL_TB"))) +
       facet_wrap(~celltype, scales = "free",nrow=3)+
       scale_color_manual(values=c('#000000','#000000','#000000')) +
       theme(axis.text.x = element_text(angle = 15, hjust = 1))+
       theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), legend.position = 'none',
        axis.text.y = element_text(color='black'), axis.text.x = element_text(color='black'))
  ggsave(paste0("/05_Result/07_co_lung/01_ref/09_gene_score/all/",object[i],".pdf"), plot = p, width = 14, height = 7.5)
}
