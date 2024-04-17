relative_vs_absolute_plot <- function(ps, tax_target, sample_format=FALSE) {
  if (sample_format == TRUE) {
    # Renaming sample names
    sample_names(ps) <- sample_names(ps) %>% 
      stringr::str_extract("(^[^_]*_[^_]*_?\\d)|(C\\d?)")
  }
  # Creates barplot of total sum of counts per sample
  df_labels <- ps %>% sample_data %>% unclass() %>% as.data.frame()
  nASVs_plot <- df_labels %>% 
    ggplot(mapping = aes(x = row.names(df_labels),
                         y = nASVs)) +
    geom_bar(stat = "identity") +
    geom_col(just = 0.5) +
    geom_text(mapping = aes(label = nASVs),
              hjust = -0.1, 
              fontface = "bold") +
    ylim(0, max(df_labels$nASVs)*1.10) +
    coord_flip() +
    theme(text=element_text(size=20),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    labs(title = "Total counts per sample",
         y = "total counts") +
    xlab(NULL)
  
  # Creates stacked barplot of relative abundances per sample
  per_group_composition <- ps %>%
    tax_fix(anon_unique = FALSE, 
            verbose = FALSE,
            unknowns = c("uncultured")) %>% 
    comp_barplot(tax_level = tax_target, 
                 n_taxa = 15, 
                 merge_other = FALSE, 
                 sample_order = "asis") +
    coord_flip() +
    guides(fill = guide_legend(ncol = 2, title.position = "top")) +
    theme(
      text=element_text(size=20),
      legend.position = "bottom", legend.key.size = unit(5, "mm")) +
    labs(title = "Relative abundance per sample",
         x = "Sample IDs",
         y = "Relative Abundance")
  
  # combines plot
  comb_plot <- (per_group_composition + nASVs_plot) +
    plot_layout(axis_titles = "collect")
  
  return(comb_plot)
}
