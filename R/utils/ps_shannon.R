ps_shannon <- function(ps, df_shannon, sample_order=NULL, col_name = "RANKSTAT_treatment", Brewer.palID="Set2") {
  # Fetch Shannon file and metadata
  meta_tab <- get_meta(ps)
  
  # Pivots into long table
  df_shannon_long <- df_shannon %>% 
    pivot_longer(cols = starts_with("depth-"),
                 names_to = "iters",
                 values_to = "alpha_div")
  
  # Adds relevant metadata
  df_shannon_final <- add_metadata(df_long = df_shannon_long,
                                   meta_tab = meta_tab,
                                   meta_col.id = "SAMPLE.ID",
                                   meta_col.add = c(col_name)
  )
  
  # Creating a color palette
  unique_groups <- unique(meta_tab[[col_name]])
  chosen_palette <- RColorBrewer::brewer.pal(length(unique_groups), Brewer.palID)
  colors <- stats::setNames(chosen_palette, unique_groups)
  
  if (!is.null(sample_order)) {
    df_shannon_final[[col_name]] <- factor(df_shannon_final[[col_name]], levels=sample_order)
  }
  # Creates shannon plot
  shannon_plot <- df_shannon_final %>%
    ggplot(mapping = aes(x = base::get(col_name, df_shannon_final),
                         y = alpha_div)) +
    geom_violin(width = 1.4, aes(fill = base::get(col_name, df_shannon_final))) +
    geom_boxplot(width = 0.1) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 12, color = "black")) + 
    scale_fill_manual(name = "", 
                      values = colors) +
    labs(title = NULL,
         subtitle = paste0("selected column: ", col_name),
         x = "sample groups",
         y = "Shannon Index")
  
  return(shannon_plot)
}
