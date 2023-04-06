library(Seurat)
library(monocle)
library(htmlwidgets)
library(ggpubr)
library(pheatmap)
library(colorRamps)
library(dplyr)
library(RColorBrewer)
library(magrittr)

options(stringsAsFactors=F)

#### object is an S4 object processed by seurat software
object <- readRDS("/05_Result/07_co_lung/01_ref/04_celltype/Combination_celltype.rds")
Idents(object) <- "Time"

#### input DEGs list
deg_list <- as.data.frame(matrix(nrow=0,ncol=9))
for (i in c("AL_TB-Ctrl","L_TB-AL_TB","L_TB-Ctrl")){
	DEGs_scRNA.list <- read.csv(paste0('/05_Result/07_co_lung/01_ref/05_DEGs/',i,'_deg.csv'),header=T)
	deg_list=rbind(deg_list,DEGs_scRNA.list)
}
time_select_gene <- unique(as.vector(deg_list$gene))

object_matrix <- GetAssayData(object = object)[time_select_gene, ]

#### Make all time points have the same number of cells by downsampling
time <- c("Ctrl", "AL_TB", "L_TB") ### time point
CellID <- c()
Sample_ID <- c()
Sample_Size <- 500 ### Set the sample size
for (m in time) {
  sub_sample <- subset(object@meta.data,group%in%m) %>% rownames(.) %>% as.vector()
  set.seed(176)
  sample_id <- sample(x=sub_sample,size=Sample_Size,replace=T) ### Set the sample size for each sample 
  Sample_ID <- c(Sample_ID,sample_id)
  cellid <- paste0(m,"_",c(1:Sample_Size))
  CellID <- c(CellID,cellid)
}

anno_cell_sample <- data.frame(Sample_ID,CellID)

all_sample_matrix <- c()
all_cell_num <- length(time)*Sample_Size
for (i in 1:all_cell_num) {
  if (i==1) {
    all_sample_matrix <- object_matrix[,which(colnames(object_matrix)%in%Sample_ID[i])] %>% as.matrix()
  }else{
    sample_matrix <- object_matrix[,which(colnames(object_matrix)%in%Sample_ID[i])] %>% as.matrix()
    all_sample_matrix <- cbind(all_sample_matrix,sample_matrix)
  }
}

colnames(all_sample_matrix) <- CellID

meta <- data.frame(Sample_ID=rownames(object@meta.data), Time=object@meta.data$group)
for (i in 1:all_cell_num) {
  if (i==1) {
    all_sample_meta <- meta[which(meta$Sample_ID%in%Sample_ID[i]),] %>% as.data.frame()
  }else{
    sample_meta <- meta[which(meta$Sample_ID%in%Sample_ID[i]),] %>% as.data.frame()
    all_sample_meta <- rbind(all_sample_meta,sample_meta)
  }
}
all_sample_meta <- merge(all_sample_meta, anno_cell_sample, by="Sample_ID") %>% .[!duplicated(.),]
rownames(all_sample_meta) <- as.vector(all_sample_meta$CellID)
all_sample_meta$CellID <- factor(all_sample_meta$CellID, levels = CellID)
all_sample_meta <- all_sample_meta[order(all_sample_meta$CellID,decreasing = F),]

count_matrix <- all_sample_matrix

### Smoothing single-cell sequencing data

##create pd
pd <- all_sample_meta
##create fd
fd <- as.data.frame(rownames(all_sample_matrix))
colnames(fd) <- "gene_short_name"
rownames(fd) <- fd$gene_short_name
my_cds <- newCellDataSet(count_matrix,
                         phenoData = new("AnnotatedDataFrame", data = pd),
                         featureData = new("AnnotatedDataFrame", data = fd))
DelayedArray:::set_verbose_block_processing(TRUE)

#Accessors and generic functions used in the context of count datasets
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)

time_select_gene <- as.data.frame(time_select_gene)
colnames(time_select_gene) <- "gene"

length.out <- nrow(pData(my_cds))
times.df <- seq(0.1,100,length.out = length.out)
pData(my_cds)$age_time <- times.df


hclust_method = "ward.D2"
num_clusters = 5 #### Set the number of clusters
scale_max = 3
scale_min = -3
cores = 20

num_clusters <- min(num_clusters, nrow(my_cds))
pseudocount <- 1
newdata <- data.frame(age_time = seq(min(pData(my_cds)$age_time),
                                     max(pData(my_cds)$age_time), length.out = all_cell_num))

### Fit smooth spline curves and return the response matrix
smooth_mat <- genSmoothCurves(my_cds, cores = cores, trend_formula = "~sm.ns(age_time, df=3)",
                     relative_expr = T, new_data = newdata)
smooth_mat = smooth_mat[!apply(smooth_mat, 1, sum) == 0, ]
smooth_mat = vstExprs(my_cds, expr_matrix = smooth_mat)
smooth_mat = smooth_mat[!apply(smooth_mat, 1, sd) == 0, ]
smooth_mat = Matrix::t(scale(Matrix::t(smooth_mat), center = TRUE))
smooth_mat = smooth_mat[is.na(row.names(smooth_mat)) == FALSE, ]
smooth_mat[is.nan(smooth_mat)] = 0
smooth_mat[smooth_mat > scale_max] = scale_max
smooth_mat[smooth_mat < scale_min] = scale_min

heatmap_matrix <- smooth_mat
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

bks <- seq(-3.1, 3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)
ph <- pheatmap(heatmap_matrix, useRaster = TRUE, cluster_cols = FALSE,
               cluster_rows = TRUE, show_rownames = FALSE, show_colnames = FALSE,
               clustering_distance_rows = row_dist, clustering_method = hclust_method,
               cutree_rows = num_clusters, silent = TRUE, filename = NA,
               breaks = bks, border_color = NA, color = hmcols)
ggsave("/05_Result/07_co_lung/01_ref/05_DEGs/time/ph_plot_6.pdf", plot = ph, width = 5, height = 6)
annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, num_clusters)))

feature_label <- as.character(fData(my_cds)[row.names(heatmap_matrix),"gene_short_name"])
feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
row.names(heatmap_matrix) <- feature_label
colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))

annotation_col <- data.frame(Time = rep(time,each=Sample_Size))
annotation_col$Time <- factor(annotation_col$Time,levels = time)
row.names(annotation_col) <- 1:all_cell_num

cluster_cols <- colorRampPalette(brewer.pal(8,"Dark2"))(num_clusters)
names(cluster_cols) <- 1:num_clusters

ann_colors = list(Time = c("Ctrl"='#c9caca',
                           "AL_TB"='#a0c3d5',
                           "L_TB"='#4b2415'),
                  Cluster = cluster_cols)

ph_res <- pheatmap(heatmap_matrix[, ], useRaster = TRUE, cluster_cols = FALSE,
                   cluster_rows = TRUE, show_rownames = FALSE,
                   show_colnames = FALSE, clustering_distance_rows = row_dist,
                   clustering_method = hclust_method, cutree_rows = num_clusters,
                   annotation_row = annotation_row, annotation_col = annotation_col,
                   annotation_colors = ann_colors,annotation_names_row = FALSE,annotation_names_col = FALSE,
                   treeheight_row = 0, breaks = bks, fontsize = 8, color = hmcols,
                   border_color = NA, silent = TRUE, filename = NA)
pdf("/05_Result/07_co_lung/01_ref/05_DEGs/time/Time_dependent_DEGs_heatmap_6.pdf",height = 6,width = 5)
ph_res
dev.off()

### Get clustering for each gene
clusters <- cutree(ph_res$tree_row, k = num_clusters)
clustering <- data.frame(clusters)
clustering[,1] <- as.numeric(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$Gene <- rownames(clustering)
clustering <- clustering[order(clustering$Gene_Clusters,decreasing = F),]

write.csv(clustering,"/05_Result/07_co_lung/01_ref/05_DEGs/time/Time_dependent_DEGs_clustering_6.csv",row.names = F)

#####曲线图######
library(reshape2)
for (i in c(1:5)){
	tmp <- subset(clustering,Gene_Clusters==i)
	gene <- tmp$Gene
	newdata <- heatmap_matrix[gene,]
	df <- melt(newdata)
	colnames(df) <- c("gene","order",'expression')
	
	p <- ggplot(df,aes(x=order, y=expression))+
		geom_smooth(formula = y~poly(x,15),method = lm,size = 1,alpha = 1,color="#EE3B3B",se=FALSE)+ 
		stat_smooth(aes(group=gene),formula = y~poly(x,15),method= lm, size = 0.02,linetype=2, alpha = 0.05,color="#6495ED",se=FALSE)+
		scale_x_continuous(limits=c(1,1500), breaks=seq(250, 1250, 500),label = c("Ctrl","AL_TB","L_TB"))+
		geom_vline(xintercept=c(1,500,1000,1500),lty=4,col="#000000",lwd=0.6)+
		theme(panel.grid=element_blank())+
		theme_bw()
	pdf(paste0("/05_Result/07_co_lung/01_ref/05_DEGs/time/cluster_",i,".pdf"),height = 3,width = 3.5)
	p
	dev.off()
}
