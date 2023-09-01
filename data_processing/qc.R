knee_df = function(mtx,expt_name){
    df = as.data.frame(colSums(mtx))
    colnames(df) = c("nUMI")
    df <- tibble(total = df$nUMI,
               rank = row_number(dplyr::desc(total))) %>%
    distinct() %>%
    arrange(rank)
    df$experiment = expt_name
    out = df
}