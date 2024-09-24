fold_plot <- function(df, X, title, method, taxa_labels = FALSE, pvalues = NULL, pvalues.col) {
  if (method == "barplot") {
    plt <- df %>% 
      aggregate(formula(paste(X, " ~ Taxa")), sum) %>% 
      ggplot(mapping = aes(x = base::get(X),
                           y = Taxa,
                           fill = base::get(X))) +
      geom_bar(stat = "identity")
    if (!is.null(pvalues)) {
      plt <- plt +
        geom_text(aes(label = ifelse(!is.na(pvalues[, pvalues.col]) & pvalues[, pvalues.col] < 0.05, "*", "")),
                  fontface = "bold",
                  position = position_dodge(width = 1),
                  size = 6)
    }
  } else if (method == "boxplot") {
    plt <- df %>% 
      ggplot(mapping = aes(x = base::get(X),
                           y = Taxa)) +
      geom_boxplot()
  }
  
  plt <- plt +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=12),
          axis.text = element_text(size=12),
          text = element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          legend.position = "none",
          axis.title.y = element_blank(),
          strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))
  if (!taxa_labels) {
    plt <- plt + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  plt <- plt +
    scale_fill_gradient2(name = paste0("log2( A / B )"),
                         low = "blue",
                         mid = "white",
                         high = "red",
                         na.value = "grey80") +
    scale_y_discrete(limits = rev(levels(as.factor(df$Taxa)))) +
    labs(x = NULL,
         y = "Taxa") +
    ggtitle(title)
  
  return(plt)
}