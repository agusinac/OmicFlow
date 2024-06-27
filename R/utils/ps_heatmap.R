# Creates 2-fold heatmap from phyloseq object
# 
# Requires: paired patient data
# TO DO: add non-paired option
# 
# 
ps_heatmap <- function(ps, taxa_n = 20, taxa_rank = "Genus", col_id = "PATIENT.ID", col_group = "RANKSTAT_treatment", condition_A, condition_B) {
  # Fetch otus and filter out less representative Genera
  otu_tab <- get_otu(ps, tax_target = taxa_rank, top_n = taxa_n) %>% t()
  meta_tab <- get_meta(ps)
  
  # log2 scale
  otu_tab.trans <- otu_tab %>%
    as.matrix() %>% 
    as("sparseMatrix") %>% 
    log2() %>% 
    as.matrix()
  
  # replace infinite by zero
  otu_tab.trans[!is.finite(otu_tab.trans)] <- 0
  
  # Creates dataframe & orders taxa abundances
  df <- as.data.frame(otu_tab.trans)
  df$Taxa <- rownames(df) 
  df$Total <- rowSums(df[, -ncol(df)])
  df <- df[order(-df$Total, decreasing=TRUE), ]
  
  # Pivot longer table by specified id.vars
  df_melt <- reshape2::melt(df, id.vars = c("Taxa", "Total"))
 
  # Removes uncultured genera, cleans sample names
  df_final <- df_melt %>% 
    rowwise() %>% 
    mutate(
      # Collects sample group names from metadata 
      group = meta_tab[[ {{ col_group }} ]][stringr::str_detect(as.character(meta_tab$SAMPLE.ID), as.character(variable))],
      # Metadata should have a column of identical sample names for paired samples!
      sample.id = meta_tab[[ {{ col_id }} ]][stringr::str_detect(as.character(meta_tab$SAMPLE.ID), as.character(variable))]
    ) %>% 
    group_by(sample.id) %>% 
    filter(any(group %in% {{condition_A}}) & any(group %in% {{condition_B}}))
  
  # Computes log2(A) - log2(B)
  # Supports multiple inputs for A and B.
  # For example A = T1, T2 and B = H1, H2
  for (i in 1:length(condition_A)) {
    name <- paste0("diff_",i)
    df_diff <- df_final %>% 
      group_by(sample.id, Taxa) %>% 
      summarise(!!name := value[group == condition_A[i]] - value[group == condition_B[i]])
    
    df_final <- df_final %>% 
      left_join(df_diff, by = c("Taxa", "sample.id"))
  }
  # Collects diff_ columns
  diff_columns <- sum(grepl("^diff_", names(df_final)))
  
  # Generate heatmap plot with df_diff data
  heatmap_plot <- df_final %>% 
    ggplot(mapping = aes(x = as.factor(sample.id),
                         y = Taxa))
  
  # If there is only one column uses default settings
  if (diff_columns == 1) {
    heatmap_plot <- heatmap_plot +
      geom_tile(aes(fill = diff_1))
  } else {
    # Adds geom_tile for number of diff_columns
    for (i in 1:diff_columns) {
      if (i == 1) {
        heatmap_plot <- heatmap_plot +
          geom_tile(aes(fill = !!sym(paste0("diff_", i))), width = 0.45)
      } else {
        heatmap_plot <- heatmap_plot +
          geom_tile(aes(fill = !!sym(paste0("diff_", i))), width = 0.45,
                    position = position_nudge(x = 0.5))
      }
      
    }
  }
  # Finishes heatmap plot
  heatmap_plot <- heatmap_plot +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=12),
          axis.text = element_text(size=12),
          text = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
    scale_fill_gradient2(name = paste0("log2( A / B )"),
                         low = "blue",
                         mid = "white",
                         high = "red",
                         na.value = "grey80") +
    scale_y_discrete(limits = rev(levels(as.factor(df_final$Taxa)))) +
    labs(x = "Samples", 
         y = "Taxa")
  
  # Fetch otu table for boxplot and reshapes into long table
  stats_tab <- as.data.frame(otu_tab)
  stats_tab$Taxa = rownames(stats_tab)
  
  # Pivot longer
  stats_melt <- reshape2::melt(stats_tab, id.vars = c("Taxa"))
  
  # Validate numeric zero's instead of NAs
  stats_melt$value <- as.numeric(stats_melt$value)
  stats_melt$value[is.na(stats_melt$value)] <- 0
  
  stats_final <- stats_melt %>% 
    rowwise() %>% 
    mutate(
      # Collects sample group names from metadata 
      sample.id = meta_tab[[ {{ col_id }} ]][stringr::str_detect(as.character(meta_tab$SAMPLE.ID), as.character(variable))],
      group = meta_tab[[ {{ col_group }} ]][stringr::str_detect(as.character(meta_tab$SAMPLE.ID), as.character(variable))]
    ) %>% 
  group_by(sample.id) %>% 
  filter(any(group %in% {{condition_A}}) & any(group %in% {{condition_B}}))
  
  # Creates boxplot from relative abundances
  rel_abun_plot <- stats_final %>% 
    ggplot(mapping = aes(x = value,
                         y = Taxa)) +
    geom_boxplot() +
    facet_wrap(~group, ncol = length(condition_A) + length(condition_B)) +
    theme_bw() +
    theme(text=element_text(size=12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.x = unit(1, "lines")) +
    # scale_x_continuous(trans = scales::log_trans()) +
    scale_y_discrete(limits = rev(levels(as.factor(df_final$Taxa)))) +
    labs(x = "Log10( Rel. Abun. )")
  
  # Combines plots
  comb_plot <- (heatmap_plot + rel_abun_plot) +
    plot_layout(guides = "collect",
                axes = "collect",
                widths = c(6,4))
  
  return(comb_plot)
}
