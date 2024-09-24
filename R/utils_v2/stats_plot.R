stats_plot <- function(data, X, Y, Label, Y_title, plot.title) {
  return(
    data %>%
      ggplot(mapping=aes(x = base::get(X, data),
                         y = base::get(Y, data),
                         label = base::get(Label, data))) +
      geom_bar(stat = "identity",
               fill = "blue") +
      geom_label(nudge_y = 0) +
      labs(title = plot.title,
           subtitle = "P adjusted significant scores are shown above each bar",
           x = "groups",
           y = Y_title) +
      theme_bw()
  )
}