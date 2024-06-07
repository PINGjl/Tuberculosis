#' 
#' @description 获取数据
#' @param .data seurat 对象
#' @param time_series 时间序列变量列
#' @param sample_number 每个时间点抽样数量
#' 
get_log2cpm_mat <- function(.data,
    time_series,
    ...,
    gene_list = NULL,
    sample_number = 200,
    assay = "RNA",
    replace = FALSE) {
    group_vars = enquos(..., .named = TRUE)

    time_series <- enquo(time_series)

    mat <- GetAssayData(.data, slot = "counts", assay = assay)
    if (!is.null(gene_list)) {
        mat <- mat[gene_list, ]
    }
    mat <- mat[rowSums(mat) != 0, ]
    mat_col_sums <- colSums(mat)

    cell_list <- .data@meta.data %>%
        as_tibble(rownames = "id") %>%
        inner_join(
            tibble(
                id = names(mat_col_sums),
                mat_col_sums = unname(mat_col_sums)
            )
        ) %>%
        filter(mat_col_sums != 0) %>%
        group_by(!!time_series) %>%
        sample_n(sample_number, replace = replace) %>%
        ungroup() %>%
        arrange(!!time_series)

    cell_anno <- cell_list %>%
        mutate(time_order := !!time_series) %>%
        select(id, time_order, !!time_series, !!!group_vars)

    if(typeof(cell_anno$time_order) != "integer") {
        cell_anno$time_order <- factor(cell_anno$time_order)
    }
        
    print("time order:")
    print(levels(cell_anno$time_order))

    mat <- mat[,cell_list$id]
    mat <- mat[rowSums(mat) != 0, ]

    mat.cpm <- apply(mat, 2, function(x) { 1e6 * x / sum(x) })
    mat.log2cpm <- apply(mat.cpm, 2, function(x) { log2(x + 1) })
    mat.scaled <- t(scale(t(mat.log2cpm)))

    print("log2cpm matrix:")
    print(dim(mat.log2cpm))

    list(
        "mat.log2cpm" = mat.log2cpm,
        "mat.scaled" = mat.scaled,
        "cell_anno" = cell_anno
    )

}

#' 
#' @description  平滑数据 loess
#' 
#' @param mat.obj 数据矩阵
#' @param span 平滑参数
#' @param scale_factor 阈值缩放因子
#' 
smooth_data <- function(mat.obj, span = 1, adj.r.squared.cutoff = 0) {

    mat.smooth <- mat.obj$mat.scaled %>%
        apply(1, function(x) {
            predict(loess(x ~ seq(1, length(x)), span = span))
        }) %>%
        t()

    colnames(mat.smooth) <- colnames(mat.obj$mat.scaled)
    mat.obj$mat.smooth <- mat.smooth

    res_t <- mat.obj$mat.scaled %>%
        apply(1, function(exp) {
            y <- unname(exp)
            x <- seq_along(exp)

            line.model <- lm(y ~ x)

            if( summary(line.model)$coefficients[2,2] < 0.05 & 
                summary(line.model)$adj.r.squared > adj.r.squared.cutoff
            ) {
                res <- line.model$coefficients[2]
            } else {
                res <- NaN
            }

            res
        })


    # reorder
    mat.obj$gene_anno <- tibble(
        id = names(res_t),
        slope = unname(res_t)
    ) %>%
    mutate(
        type = case_when(
            slope > 0 ~ "Up",
            slope < 0 ~ "Down",
            TRUE ~ "Other"
        )
    ) %>%
    arrange(-slope)

    print(table(mat.obj$gene_anno$type))

    mat.obj
}

# 聚类热图

