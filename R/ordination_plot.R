ordination_plot <- function(df_pcs, pcs, pair, dist.metric) {
  # Creates labels
  plot_title = paste0("Distance metric used: ", dist.metric)
  if (!is.null(pcs$eig_norm)) {
    x_label = paste0(pair[1], " (", round(as.numeric(pcs$eig_norm[1]), 2), "%)")
    y_label = paste0(pair[2], " (", round(as.numeric(pcs$eig_norm[2]), 2), "%)")
  } else {
    x_label = paste0(pair[1])
    y_label = paste0(pair[2])
  }
  
  return(
    df_pcs %>% 
      ggplot(mapping = aes(x = .data[[ pair[1] ]],
                           y = .data[[ pair[2] ]],
                           color = groups,
                           linetype = groups)) +
      geom_point(alpha = 5) +
      stat_ellipse(type = "norm") +
      theme_bw() +
      theme(text=element_text(size=14),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            axis.text = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.text.x = element_text(size=12)
      ) +
      labs(title = NULL,
           subtitle = NULL,
           x = x_label,
           y = y_label
      )
  )
}