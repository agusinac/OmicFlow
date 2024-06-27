ps_shannon <- function(ps, shannon_file, sample_order=NULL, metadata.columns = c("RANKSTAT_treatment"), Brewer.palID="Set2") {
  # Fetch Shannon file and metadata
  df_shannon <- read_csv(file = shannon_file)
  meta_tab <- get_meta(ps_abs_bac)
  
  # Pivots into long table
  df_shannon_long <- df_shannon %>% 
    pivot_longer(cols = starts_with("depth-"),
                 names_to = "iters",
                 values_to = "alpha_div")
  
  # Adds relevant metadata
  df_shannon_final <- add_metadata(df_long = df_shannon_long,
                                   meta_tab = meta_tab,
                                   meta_col.id = "SAMPLE.ID",
                                   meta_col.add = metadata.columns
  )
  
  # Creating a color palette
  unique_groups <- unique(meta_tab$RANKSTAT_treatment)
  chosen_palette <- RColorBrewer::brewer.pal(length(unique_groups), Brewer.palID)
  colors <- stats::setNames(chosen_palette, unique_groups)
  
  if (!is.null(sample_order)) {
    df_shannon_final$RANKSTAT_treatment <- factor(df_shannon_final$RANKSTAT_treatment, levels=sample_order)
  }
  # Creates shannon plot
  shannon_plot <- df_shannon_final %>%
    ggplot(mapping = aes(x = RANKSTAT_treatment,
                         y = alpha_div)) +
    geom_violin(width = 1.4, aes(fill = RANKSTAT_treatment)) +
    geom_boxplot(width = 0.1) +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 12, color = "black")) + 
    scale_fill_manual(name = "", 
                      values = colors) +
    labs(title = NULL,
         x = "sample groups",
         y = "Shannon Index")
  
  return(shannon_plot)
}