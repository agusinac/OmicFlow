composition_plot <- function(data, palette, feature_rank = "Genus", title_name = "", group_by = FALSE) {
  # Generates a stacked barplot as base with custome palette
  if (group_by != FALSE) {
    plt <- data %>% 
      ggplot(mapping = aes(y = value,
                           x = base::get(group_by, data),
                           fill = Taxa))
  } else {
    plt <- data %>% 
      ggplot(mapping = aes(y = value,
                           x = variable,
                           fill = Taxa))
  }
  # Required for stacked barplot
  plt <- plt +
    geom_bar(position = "fill",
             stat = "identity")
  
  if (group_by == FALSE) {
    plt <- plt +
      coord_flip()
  }
  plt <- plt +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 13),
      axis.text.x = element_text(angle = 90, size = 12,
                                 vjust = 0.5, hjust=1,
                                 colour = "black"),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12, vjust=0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(colour = "black", size = 12)
    )
  
  if (group_by == FALSE) {
    plt <- plt +
      scale_x_discrete(limits = rev(levels(as.factor(data$variable))))
  }
  plt <- plt +
    scale_fill_manual(values = palette, name = feature_rank) +
    labs(y = "Rel. Abun.",
         x = NULL) +
    ggtitle(title_name)
  
  return(plt)
}
