library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(monocle)
library(ggridges)
library(pheatmap)
library(viridis)
library(ggpubr)
library(plyr)
library(scales)
library(showtext)
showtext_auto()

set.seed(2)
options(future.globals.maxSize = 26000 * 1024^2)
setwd("/05_Result/07_co_lung/01_ref/12_pseudotime/Mono_Mac/")
cell_order<-c("Classical_MC","Nonclassical_MC","Mac2")
cell_color<-c("Classical_MC"="#e2a461","Nonclassical_MC"='#fb832d',"Mac2"="#ebc3d5") 
age_order<-c("Ctrl","AL_TB","L_TB")
age_color<-c("Ctrl"='#47B1B6', "AL_TB"='#CEB7B3', "L_TB"='#E6949A')
my_theme <- theme(#legend.position = "NULL",axis.ticks.x = element_blank(),
  axis.title = element_blank(),
  panel.border = element_blank(),panel.grid=element_blank(),panel.background=element_blank(),
  axis.text.x = element_text(color = "black",angle = 90,vjust = 0.5,hjust = 1),axis.text.y=element_text(color="black"),
  axis.line = element_line(size = 0.5))

lung <- readRDS(file = "/05_Result/07_co_lung/01_ref/04_celltype/Combination_celltype.rds")
germ_part <- subset(lung, celltype %in% cell_order)
rm(lung)
# 对表达矩阵和meta矩阵进行downsample ####
exp_mat<-as.matrix(GetAssayData(germ_part,slot="counts"))
sample_mat<-exp_mat
sample_mat<-sample_mat[rowSums(sample_mat)!=0,]
sample_germ_part<-germ_part[[]]
# 创建CDS对象 ####
pd<-data.frame(sample_germ_part)
fd<-data.frame(gene_short_name=rownames(sample_mat),row.names=rownames(sample_mat))
pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)
germ_part_monocle <- newCellDataSet(as(as.matrix(sample_mat), "sparseMatrix"),
                                    phenoData = pd, 
                                    featureData = fd,
                                    lowerDetectionLimit=0.5,
                                    expressionFamily=negbinomial.size())

germ_part_monocle<-estimateSizeFactors(germ_part_monocle)
germ_part_monocle<-estimateDispersions(germ_part_monocle)
head(fData(germ_part_monocle))
# 过滤基因 ####
germ_part_monocle<-detectGenes(germ_part_monocle,min_expr=1) 
length(rownames(germ_part_monocle)) # 此时有17273个基因
expressed_genes<-rownames(subset(fData(germ_part_monocle),num_cells_expressed>=50))
length(expressed_genes)
germ_part_monocle<-germ_part_monocle[expressed_genes,]
#过滤掉在小于50个细胞中表达的基因，还剩10357个基因。
saveRDS(germ_part_monocle,"Mono_Mac_celltype_monocle_filter.rds")

# 找高变异基因并排序 ####
germ_part_monocle <- readRDS("Mono_Mac_celltype_monocle_filter.rds")
celltype_marker <- read.csv("celltype_marker.csv")
celltype_marker <- celltype_marker %>% group_by(cluster) %>% top_n(n = 150,wt = avg_log2FC)
virgene <- unique(celltype_marker$gene)
germ_part_monocle<-setOrderingFilter(germ_part_monocle,virgene)

pdf("ordering_gene.pdf")
plot_ordering_genes(germ_part_monocle)
dev.off()
germ_part_monocle<-reduceDimension(germ_part_monocle, max_components=2, method='DDRTree')
germ_part_monocle<-orderCells(germ_part_monocle)
saveRDS(germ_part_monocle,"Mono_Mac_celltype_monocle_marker150.rds")

#重排序并作图 ####
germ_part_monocle<-readRDS("Mono_Mac_celltype_monocle_marker150.rds")
#germ_part_monocle<-orderCells(germ_part_monocle,root_state = 5)
germ_part_monocle$group<-factor(germ_part_monocle$group,levels=c("Ctrl","A_TB","F_TB"))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
pdf("plot_pseudotime_by_pseudotime.pdf")#,width=6,height=5
p=plot_cell_trajectory(germ_part_monocle, color_by = "Pseudotime",cell_size=1.5,show_branch_points = F)
print(p)
dev.off()
pdf("plot_pseudotime_by_celltype.pdf")#,width=6,height=5
p=plot_cell_trajectory(germ_part_monocle, color_by = "celltype",cell_size=1.5,show_branch_points = F)+scale_color_manual(values=cell_color)
print(p)
dev.off()
pdf("plot_pseudotime_by_group.pdf",width=6,height=5)
p=plot_cell_trajectory(germ_part_monocle, color_by = "group",cell_size=1.5,show_branch_points = F)+scale_color_manual(values= age_color)
print(p)
dev.off()

# 查看TFs在伪时间轨迹上的相对表达 ####
setwd("/05_Result/07_co_lung/01_ref/12_pseudotime/Mono_Mac/final/")
tfs_list <- c("PPARG","BHLHE41","TCF7L2")
pic <- plot_genes_in_pseudotime(germ_part_monocle[tfs_list,], color_by = "celltype")
pdf("expression_level_tf_Celltype.pdf")
print(pic)
dev.off()
# 单细胞轨迹的“分支”分析 ####
BEAM_res <- BEAM(germ_part_monocle[virgene,], branch_point = 1) 
#BEAM_res <- BEAM(germ_part_monocle, branch_point = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name","pval", "qval")]
head(BEAM_res)
write.csv(BEAM_res,"BEAM_res.csv",row.names = F)
my_color <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "RdYlGn")))(n=62)
p <- plot_genes_branched_heatmap(germ_part_monocle[row.names(BEAM_res),],
		branch_point = 1, #绘制的是哪个分支
		num_clusters = 5, #分成几个cluster，根据需要调整
		cores = 1,
		use_gene_short_name = T,hmcols = my_color,
		show_rownames = T,return_heatmap = TRUE, branch_colors = c('#DA561D','#0E5B49','#7C5674'))
ggsave("branched_time_heatmap.pdf",p$ph_res, width = 6.5, height = 16)
# 取出cluster及其基因
clusters <- cutree(p$ph_res$tree_row, k = 5)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
write.csv(clustering,"branched_time_heatmap.csv",quote=F)
