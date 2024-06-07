library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(Seurat)
library(pheatmap)
setwd("/01_result/07_lung/00_scRNA_new/05_DEGs/up_up-down_down")

#### input DEGs list
deg_list <- as.data.frame(matrix(nrow=0,ncol=9))
for (i in c("A_TB-Ctrl","F_TB-A_TB","F_TB-Ctrl")){
        DEGs_scRNA.list <- read.csv(paste0('/01_result/07_lung/00_scRNA_new/05_DEGs/',i,'_deg_025.csv'),header=T)
        deg_list=rbind(deg_list,DEGs_scRNA.list)
}
deg <- unique(deg_list[,c('gene','celltype')])

lung <- readRDS(file = "/01_result/07_lung/00_scRNA_new/04_celltype/Combination_celltype.rds")
lung$group <- factor(lung$group,
                             levels=c("Ctrl", "A_TB", "F_TB"), ordered=TRUE)

source("loess_smooth_heatmap.R")

celltypes <- c("AT1","AT2","Cil","Clu_Gob", "EC","Stroma",'Mast',"cMC","nMC","Mac1","Mac2","Neu","DC",'CD4+TC',"CD8+TC",'NK',"NKT","BC","Pla",'MKI67+_pro')

#deg$celltype <- str_replace(deg$celltype, " ", "_")

lung$celltype <- factor(lung$celltype,
                             levels=celltypes, ordered=TRUE)

ct_cluster_list <- celltypes %>%
    map(~{
        message(.x)
        ct_data <- lung %>%
            subset(celltype == .x)
        ct_gene_list <- deg %>%
            subset(celltype == .x)

        ct_mat <- get_log2cpm_mat(ct_data,
            group, celltype,
            gene_list = ct_gene_list$gene,
            sample_number = 200)
        ct_mat <- smooth_data(ct_mat)

        setNames(list(ct_mat), .x)
    }) %>%
    unlist(recursive = F)

diff_gene_mat <- ct_cluster_list %>%
    imap_dfr(~{
        mat <- .x
        ct <- .y

        gene_list <- mat$gene_anno %>%
            filter(type != "other")
        mat_up <- mat$mat.smooth[gene_list$id, ]

        colnames(mat_up) <- str_c(substr(colnames(mat_up),1,4), "-", seq(dim(mat_up)[2])) 
        # rownames(mat_up) <- str_c( .y, "_Up_", rownames(mat_up))

        mat_up %>%
            as_tibble(rownames = "id") %>%
            left_join(mat$gene_anno) %>%
            mutate(
                gene = id,
                celltype = ct,
                id = str_c(ct, "_", id))
    })

diff_gene_mat %>%
    write_csv("diff_gene_mat.csv")


diff_gene_mat %>%
    group_by(celltype, type) %>%
    summarise(n = n()) %>%
    write_csv("diff_gene_summary.csv")

source("sc_utils.R")

diff_gene_mat %>%
    select(celltype, type, gene) %>%
    filter(type == "Up") %>%
    tibble_to_lists(celltype, col_name = "gene") %>%
    write_csv("diff_gene_list_up.csv")

diff_gene_mat %>%
    select(celltype, type, gene) %>%
    filter(type == "Down") %>%
    tibble_to_lists(celltype, col_name = "gene") %>%
    write_csv("diff_gene_list_down.csv")

diff_gene_mat %>%
    select(celltype, type, gene) %>%
    write_csv("age_dependent_gene.csv")

read_csv("age_dependent_gene.csv") %>%
    filter(type != "Other") %>%
    write_csv("age_dependent_gene_up_down.csv")

##############################################

celltype_full <- c("AT1","AT2","Cil","Clu_Gob", "EC","Stroma",'Mast',"cMC","nMC","Mac1","Mac2","Neu","DC",'CD4+TC',"CD8+TC",'NK',"NKT","BC","Pla",'MKI67+_pro')
CT.col <- c('#c72e29','#d5695d','#F39B7FB2',"#E64B35B2",
"#F0E68C",'#ddd954',
'#804588',"#e2a461",'#fb832d','#cc759d',"#ebc3d5","#c0d8b9","#80bbb2",
'#a0c167','#4b9748',"#00A087B2",'#59a77f','#779eb9','#346e9a','#955025')
celltype_cols <- CT.col %>%
    setNames(celltype_full)
#celltype_cols <- celltype_cols[-c(10, 14)]

group_cols <- c("Ctrl"='#47B1B6', "A_TB"='#CEB7B3', "F_TB"='#E6949A')
SpectralColors <- colorRampPalette(colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu") %>% rev())

up_mat <- diff_gene_mat %>%
    filter(type == "Up") %>%
    select(-c(type, gene, slope)) %>%
    mutate(celltype = factor(celltype, levels = celltypes)) %>%
    group_nest(celltype) %>%
    arrange(celltype) %>%
    mutate(
        data = map(data, ~{
            d <- .x %>%
                select(-id) %>%
                dist() %>%
                hclust()
            .x[d$order,]
        })
    ) %>%
    unnest(data) %>%
    column_to_rownames("id") %>%
    select(-celltype)


up_row_anno <- diff_gene_mat %>%
    filter(type == "Up") %>%
    column_to_rownames("id") %>%
    select(celltype)

up_col_anno <- data.frame(
    Group = str_extract(colnames(up_mat), "^[YMO]"),
    row.names = colnames(up_mat)
)

up_row_anno$celltype <- factor(up_row_anno$celltype, levels = celltypes) 
gaps_row <- cumsum(table(up_row_anno$celltype))

p <- pheatmap::pheatmap(
    up_mat,
    cluster_rows = F,
    cluster_cols = F,
    show_rownames = F,
    show_colnames = F,
    annotation_names_row = F,
    annotation_row = up_row_anno,
    # annotation_col = up_col_anno,
    annotation_legend = F,
    legend = F,
    scale = "row",
    color = colorRampPalette(colors = brewer.pal(n = 11, name = "RdBu") %>% rev())(100),
    gaps_row = gaps_row[-12],
    annotation_colors = list(
        celltype = celltype_cols,
        Group = group_cols)
)

ggsave("up_gene.png", plot = p, width = 5, height = 12)
ggsave("up_gene.pdf", plot = p, width = 5, height = 12)


down_mat <- diff_gene_mat %>%
    filter(type == "Down") %>%
    select(-c(type, gene, slope)) %>%
    mutate(celltype = factor(celltype, levels = celltypes)) %>%
    group_nest(celltype) %>%
    arrange(celltype) %>%
    mutate(
        data = map(data, ~{
            d <- .x %>%
                select(-id) %>%
                dist() %>%
                hclust()
            .x[d$order,]
        })
    ) %>%
    unnest(data) %>%
    column_to_rownames("id") %>%
    select(-celltype)

down_row_anno <- diff_gene_mat %>%
    filter(type == "Down") %>%
    column_to_rownames("id") %>%
    select(celltype)

down_col_anno <- data.frame(
    Group = str_extract(colnames(down_mat), "^[YMO]"),
    row.names = colnames(down_mat)
)

gaps_row <- table(down_row_anno$celltype)
gaps_row <- cumsum(gaps_row[celltypes])

p <- pheatmap::pheatmap(
    down_mat,
    cluster_rows = F,
    cluster_cols = F,
    show_rownames = F,
    show_colnames = F,
    annotation_names_row = F,
    annotation_row = down_row_anno,
    # annotation_col = down_col_anno,
    annotation_legend = F,
    legend = F,
    scale = "row",
    color = colorRampPalette(colors = brewer.pal(n = 11, name = "RdBu") %>% rev())(100),
    gaps_row = gaps_row[-12],
    annotation_colors = list(
        celltype = celltype_cols,
        Group = group_cols)
)

ggsave("down_gene.png", plot = p, width = 5, height = 12)
ggsave("down_gene.pdf", plot = p, width = 5, height = 12)
