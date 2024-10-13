diversity_plot <- function(dt, values, col_name, palette, method) {
  return(dt %>%
    ggplot(mapping = aes(x = .data[[ col_name ]],
                         y = .data[[ values ]])) +
    gghalves::geom_half_boxplot() +
    gghalves::geom_half_point_panel(aes(color = .data[[ col_name ]])) +
    theme_bw() +
    theme(legend.position = "none",
          text=element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12)) + 
    scale_colour_manual(name = "groups", 
                        values = palette) +
    labs(title = NULL,
         subtitle = paste0("selected column: ", col_name),
         x = "sample groups",
         y = paste0("Alpha diversity: ", method))
  )
}
