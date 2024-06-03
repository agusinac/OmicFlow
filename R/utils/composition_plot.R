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
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 18),
      axis.text.x = element_text(angle = 90, size = 18,
                                 vjust = 0.5, hjust=1,
                                 colour = "black"),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15, vjust=0.5),
      legend.title = element_text(size = 15, face = "bold"),
      legend.text = element_text(size = 15, colour = "black"),
      axis.text.y = element_text(colour = "black", size = 15)
    ) +
    scale_x_discrete(limits = rev(levels(as.factor(df$variable)))) +
    scale_fill_manual(values = palette, name = tax_level) +
    labs(y = "Rel. Abun.",
         x = NULL) +
    ggtitle(title_name)
  
  return(plt)
}
