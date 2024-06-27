ps_rankstat <- function(ps) {
  # Collects rank names & taxonomy table from phyloseq object
  ranks <- phyloseq::rank_names(ps)
  taxtab <- as.data.frame(phyloseq::tax_table(ps))
  # tmp is required to fill the first row
  df_ranks <- data.frame(tmp = 0)
  for (rank in ranks) {
    df_ranks <- df_ranks %>% 
      mutate(
        !!rank := taxtab %>% select( all_of ({{ rank }}) ) %>% drop_na() %>% nrow()
      )
  }
  
  # Pivots table longer, and factorizes order of ranks
  df_ranks.long <- df_ranks %>% 
    pivot_longer(cols = ranks,
                 names_to = "variable",
                 values_to = "counts")
  df_ranks.long$variable <- factor(df_ranks.long$variable, levels = base::rev(ranks))
  
  # Creates simple barplot
  plt <- df_ranks.long %>% 
    ggplot(mapping = aes(x = variable,
                         y = counts)) +
    geom_col(fill = "grey", 
             colour = "grey15", 
             linewidth = 0.25) +
    coord_flip() +
    geom_text(mapping = aes(label = counts),
              hjust = -0.1, 
              fontface = "bold") +
    ylim(0, max(df_ranks)*1.10) +
    theme_bw() +
    labs(x = "Rank",
         y = "Number of ASVs classified")
  
  return(plt)  
}
