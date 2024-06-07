
tibble_to_lists <- function(diff_gene, ..., col_name = "GeneName", gene_filter = T, fill = "") {
    group_vars <- enquos(..., .named = TRUE)
    col <- sym(col_name)

    if (gene_filter) {
        diff_gene <- diff_gene %>%
            filter(!grepl("^ENSMFAG|^ENSMUSG|^ENSG", .[[{{col_name}}]]))
    }
    diff_gene %>%
        select(col, ...) %>%
        group_nest(...) %>%
        mutate(
            group_name = str_c(!!!group_vars, sep = "_"),
            max_length = map_dbl(data, nrow) %>% max()
        ) %>%
        select(group_name, data, max_length) %>%
        pmap_dfc(~ {
            t <- ..2[[{{col_name}}]]
            length(t) <- ..3
            t <- replace_na(t, fill)
            setNames(as_tibble(t), ..1)
        })
}

