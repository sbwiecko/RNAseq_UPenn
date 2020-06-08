profile_function <- function(data, samples) {
    # the original code chunck comes from the step2 script
    # one dedicated file for each function

    data_count <- cpm(data, log = TRUE)
    data_count_df <- as_tibble(data_count, rownames = "geneID")
    colnames(data_count_df) <- c("geneID", samples)
    data_count_df_pivot <- pivot_longer(data_count_df,
                                        cols = -1, # all except col 1
                                        names_to = "samples",
                                        values_to = "expression")

    ggplot(data_count_df_pivot) +
        aes(x = samples, y = expression, fill = samples) +
        geom_violin(trim = FALSE, show.legend = FALSE) +
        stat_summary(fun = "median",
                     geom = "point",
                     shape = 95,
                     size = 10,
                     color = "black",
                     show.legend = FALSE) +
        labs(y = "log2 expression", x = "sample",
             title = "Log2 Counts per Million (CPM)",
             caption = paste0("produced on ", Sys.time())) +
        theme_bw()
}