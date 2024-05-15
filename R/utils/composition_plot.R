composition_plot <- function(df, palette, tax_level = "Genus", title_name = "") {
  # Generates a stacked barplot as base with custome palette
  plt <- df %>% 
    ggplot(mapping = aes(y = value,
                         x = variable,
                         fill = Taxa)) +
    geom_bar(position = "fill",
             stat = "identity") +
    coord_flip() +
    theme_bw() +
    theme(
      legend.position = "right",
      text=element_text(size=12),
      axis.text.x = element_text(angle = 90, size = 12,
                                 vjust = 0.5, hjust=1,
                                 colour = "black"),
      axis.title.y = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(colour = "black", size = 12)
    ) +
    scale_fill_manual(values = palette, name = tax_level) +
    labs(y = "Rel. Abun.",
         x = NULL,
         title = "")
  
  return(plt)
}