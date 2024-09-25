fold_plot <- function(dt, X, Y, method, title = NULL, taxa_labels = FALSE, pvalues = NULL, pvalue.col = "pvalue", pvalue.threshold = 0.01) {
  if (method == "barplot") {
    new_dt <- dt %>% 
      stats::aggregate(formula(paste(X, " ~ ", Y)), sum)
    if (!is.null(pvalues)) {
      new_dt <- base::merge(new_dt, pvalues, by = Y, all.x = TRUE)
    }
    
    plt <- new_dt %>% 
      ggplot(mapping = aes(x = .data[[ X ]],
                           y = .data[[ Y ]],
                           fill = .data[[ X ]])) +
      geom_bar(stat = "identity")
    
    if (!is.null(pvalues)) {
      plt <- plt +
        geom_text(aes(label = ifelse(!is.na(new_dt[[pvalue.col]]) & new_dt[[pvalue.col]] < pvalue.threshold, "*  ", "")),
                  fontface = "bold",
                  position = position_dodge(width = 0.2),
                  size = 6)
    }
  } else if (method == "boxplot") {
    plt <- dt %>% 
      ggplot(mapping = aes(x = .data[[ X ]],
                           y = .data[[ Y ]])) +
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
    scale_y_discrete(limits = rev(levels(as.factor(dt[[ Y ]])))) +
    labs(x = paste(X),
         y = paste(Y)) +
    ggtitle(title)
  
  return(plt)
}
