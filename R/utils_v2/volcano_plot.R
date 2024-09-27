# Creates volcano plot
volcano_plot <- function(dt, feature_rank, pvalue.threshold = 0.05, logfold.threshold = 0.6) {
  # copies data.table
  tmp_dt <- data.table::copy(dt)
  
  # Creates labels for significant and non-significant differential expression
  tmp_dt[, pvalue := -log10(dt$pvalue)]
  tmp_dt[, diffexpressed := ifelse(foldchange > logfold.threshold & pvalue > -log10(pvalue.threshold), "Upregulated",
                                   ifelse(foldchange < -logfold.threshold & pvalue > -log10(pvalue.threshold), "Downregulated", "non-significant"))]
  tmp_dt[, diffexpressed_labels := ifelse(diffexpressed != "non-significant", base::get(feature_rank), "")]
  
  max_pvalue <- max(tmp_dt$pvalue)  
  
  return(
    tmp_dt %>% 
      ggplot(mapping = aes(x = foldchange,
                           y = pvalue,
                           label = diffexpressed_labels,
                           color = foldchange,
                           size = ifelse(diffexpressed != "non-significant", Mean, ""))) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
            axis.text.y = element_text(size=12),
            axis.text = element_text(size=12),
            text = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
      geom_rect(aes(xmin = -logfold.threshold, xmax = -Inf, 
                    ymin = -log10(pvalue.threshold), ymax = Inf), 
                fill = "#C8E7F1", alpha = 0.1, color = NA) +
      geom_rect(aes(xmin = logfold.threshold, xmax = Inf, 
                    ymin = -log10(pvalue.threshold), ymax = Inf),
                fill = "#FFF1F3", alpha = 0.1, color = NA) +
      annotate("text", x = -logfold.threshold*10 , y = max_pvalue+1.5, 
               label = "Significant\ndecrease",
               vjust = 2, size = 5, color = "black") +
      annotate("text", x = logfold.threshold*10 , y = max_pvalue+1.5, 
               label = "Significant\nincrease", 
               vjust = 2, size = 5, color = "black") +
      geom_vline(xintercept = c(-logfold.threshold, logfold.threshold), 
                 col = "black", linetype = 'dashed') +
      geom_hline(yintercept = -log10(pvalue.threshold), 
                 col = "black", linetype = 'dashed') + 
      scale_color_gradient2(name = "",
                            low = "blue",
                            mid = "black",
                            high = "red",
                            na.value = "grey80") +
      guides(color = "none") +
      ggrepel::geom_label_repel(show.legend = FALSE, 
                                color = "black") +
      geom_point(shape = 16, alpha = 0.5) +
      scale_size_discrete(name = "Percentage Mean Rel. Abun.") +
      labs(x = "Fold Change log2(response / non-responders)",
           y = "-log10(Adjusted P-value)") +
      ylim(0, max_pvalue + 1.5)
  )
  # removes tmp 
  base::gc(tmp_dt)
}
