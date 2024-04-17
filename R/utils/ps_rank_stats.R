# Creates a dataframe of the total number of ASVs found at different taxonomical level for all samples


ps_rank_stats <- function(ps) {
  ranks <- base::setdiff(phyloseq::rank_names(ps), c("Sequence", "ASV", "unique"))
  taxtab <- as.data.frame(phyloseq::tax_table(ps))
  
  # define internal helper function
  compute_stats_for_a_rank <- function(taxtab, rank) {
    rank_vec <- taxtab[[rank]]
    is_na <- is.na(rank_vec)
    is_stump <- base::grepl("^[A-Za-z]_{1,4}$", rank_vec)
    is_zero_char <- !nzchar(rank_vec)
    is_classified <- !is_na & !is_stump & !is_zero_char
    classified_vec <- rank_vec[is_classified]
    n_ASVs <- length(rank_vec)
    n_ASVs_classified_at_rank <- length(classified_vec)
    prop_ASVs_classified_at_rank <- n_ASVs_classified_at_rank / n_ASVs
    n_unique_at_rank <- length(unique(classified_vec))
    df <- data.frame(
      rank = rank, n_ASVs = n_ASVs,
      n_ASVs_classified_at_rank = n_ASVs_classified_at_rank,
      prop_ASVs_classified_at_rank = prop_ASVs_classified_at_rank,
      n_unique_at_rank = n_unique_at_rank
    )
    return(df)
  }
  # apply function to each rank
  rank_stats_df <- ranks %>%
    purrr::map(function(x) compute_stats_for_a_rank(taxtab, rank = x)) %>%
    purrr::list_rbind() %>%
    dplyr::mutate(rank = fct_rev(fct_inorder(rank)))
  
  return(as_tibble(rank_stats_df))
}
