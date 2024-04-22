ps_composition <- function(ps, tax_level, col_name, taxa_n=15, col_order="unique") {
  # General function for relative abundance plot
  rel_Abun_plot <- function(ps, tax_level, col_name, group_order, taxa_n) {
    plt <- target_ps %>%
      tax_fix(anon_unique = FALSE, 
              verbose = FALSE,
              unknowns = c("uncultured"))
    
    if(group_order == "unique") {
      plt <- plt %>% merge_samples(group = col_name)
    }
    plt <- plt %>% 
      comp_barplot(tax_level = tax_level, 
                   n_taxa = taxa_n, 
                   merge_other = FALSE, 
                   sample_order = group_order) +
      coord_flip() +
      guides(fill = guide_legend(ncol = 2, title.position = "top")) +
      theme(
        text=element_text(size=20),
        legend.position = "bottom", legend.key.size = unit(5, "mm")) +
      labs(x = "Sample IDs",
           y = "Rel. Abun.")
    return(plt)
  }
  
  # Compute number of ASVs
  target_ps <- ps %>% 
    ps_mutate(nASVs = sample_sums(ps))
  
  # Creates barplot of total sum of counts per sample
  meta_tab <- get_meta(target_ps)
  nASVs_plot <- meta_tab %>% 
    ggplot(mapping = aes(x = row.names(meta_tab),
                         y = nASVs)) +
    geom_bar(stat = "identity") +
    geom_col(just = 0.5) +
    geom_text(mapping = aes(label = nASVs),
              hjust = -0.1, 
              fontface = "bold") +
    ylim(0, max(meta_tab$nASVs)*1.10) +
    coord_flip() +
    theme(text=element_text(size=20),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    labs(y = "Abs. Abun.") +
    xlab(NULL)
  
  # Switch case to compute different dissimilarity matrices
  if (col_name %in% colnames(meta_tab)) {
    group_order <- switch(
      col_order,
      "all" = "asis",
      "unique" = unique(sample_data(target_ps)[[ {{ col_name }} ]])
    )
  } else warning("col_name is not specified, please input a valid column name from metadata table, all or unique")
  
  # Creates stacked barplot of relative abundances per sample
  composition_plot <- rel_Abun_plot(ps = target_ps, 
                                    tax_level = tax_level, 
                                    col_name = col_name,
                                    group_order = group_order,
                                    taxa_n = taxa_n)
  
  # combines plot
  comb_plot <- (composition_plot + nASVs_plot) + 
    plot_layout(axis_titles = "collect",
                widths = c(5,2))
  
  return(comb_plot)
}
