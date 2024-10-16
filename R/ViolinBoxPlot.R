ViolinBoxPlot <- function(dt, X, Y, diff, feature_rank, logfold.threshold = 0.06, pvalue.threshold = 0.05) {
  # Subset significant groups
  significant_subset <- dt$volcano[dt$volcano[, base::get(Y) < pvalue.threshold & (base::get(X) > logfold.threshold  | base::get(X) < -logfold.threshold)], ]
  significant_final <- dt$data[dt$data[[feature_rank]] %in% significant_subset[[feature_rank]], ]
  
  # Use significant values to create half box & violin plot
  return(
    significant_final %>% 
      ggplot(aes(y = .data[[ diff ]],
                 x = .data[[ feature_rank ]])) +
      gghalves::geom_half_boxplot(side = "l") +
      gghalves::geom_half_violin(side = "r", nudge = 0.1) +
      coord_flip() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
            axis.text.y = element_text(size=12),
            axis.text = element_text(size=12),
            text = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            legend.position = "none",
            axis.title.y = element_blank(),
            strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
      scale_x_discrete(limits = rev(levels(as.factor(significant_final[[ feature_rank ]])))) +
      labs(x = NULL,
           y = NULL)
  )
}
