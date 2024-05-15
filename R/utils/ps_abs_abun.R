ps_abs_abun <- function(ps, tax_level) { 
  # Extract sample sums from phyloseq object 
  nASVs <- get_otu(ps, tax_target = tax_level) %>%
    t() %>% colSums() %>% as.data.frame()
  
  # Creates barplot of total sum of counts per sample
  nASVs_plot <- nASVs %>% 
    ggplot(mapping = aes(x = rownames(nASVs),
                         y = nASVs$.)) +
    geom_bar(stat = "identity") +
    geom_col(just = 0.5) +
    geom_text(mapping = aes(label = nASVs$.),
              hjust = -0.1, 
              fontface = "bold") +
    ylim(0, max(nASVs$.)*1.10) +
    coord_flip() +
    theme_bw() +
    theme(text=element_text(size=12),
          legend.position = "none",
          axis.title.y = element_text(size = 12),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          plot.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    labs(y = "Abs. Abun.",
         x = NULL)
  
  return(nASVs_plot)
}