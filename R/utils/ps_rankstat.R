# Plots the total number of ASVs found at different taxonomical levels from all samples

ps_rankstat <- function(ps, title) {
  nASVs <- phyloseq::ntaxa(ps)
  rank_stats <- ps_rank_stats(ps)
  
  plot <- rank_stats %>%
    ggplot(aes(x = n_ASVs_classified_at_rank, y = rank)) +
    geom_vline(xintercept = nASVs, colour = "grey15", linewidth = 0.25) +
    geom_col(fill = "grey", colour = "grey15", linewidth = 0.25) +
    geom_text(
      mapping = aes(label = n_ASVs_classified_at_rank),
      hjust = 1.5, nudge_x = nASVs / 100, fontface = "bold"
    ) +
    scale_x_continuous(
      name = "Number of ASVs classified",
      sec.axis = sec_axis(
        trans = function(x) x / nASVs, labels = scales::label_percent()
      )
    ) +
    theme_bw()
  return(plot)
}