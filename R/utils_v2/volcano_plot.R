# Creates volcano plot
# TO be optimizes
volcano_plot <- function(dt, feature_rank, pvalue.threshold = 0.05, logfold.threshold = 0.6) {
  dt$pvalue <- -log10(dt$pvalue)
  dt$diffexpressed <- "non-significant"
  dt$diffexpressed[dt$Mean > logfold.threshold & dt$pvalue > -log10(pvalue.threshold)] <- "Upregulated"
  dt$diffexpressed[dt$Mean < -logfold.threshold & dt$pvalue > -log10(pvalue.threshold)] <- "Downregulated"
  dt$diffexpressed_labels <- ifelse(dt$diffexpressed != "non-significant", dt[[ feature_rank ]], "")
  
  return(
    dt %>% 
      ggplot(mapping = aes(x = Mean,
                           y = pvalue,
                           col = diffexpressed,
                           label = diffexpressed_labels)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
            axis.text.y = element_text(size=12),
            axis.text = element_text(size=12),
            text = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
      ggrepel::geom_label_repel() +
      geom_point(shape = 16, alpha = 0.5) +
      geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
      scale_colour_manual(values = c("Downregulated" = "blue", 
                                     "non-significant" = "black", 
                                     "Upregulated" = "red")) +
      labs(x = "Fold change (log2)",
           y = "P-value")
  )
}
