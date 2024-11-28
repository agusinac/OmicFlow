#' Volcano plot with ggplot2
#'
#' @description Creates a Volcano plot. This function is built into the \code{differential_feature_expression} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param dt A \code{data.frame} or \code{data.table}.
#' @param X A column name of a continuous variable.
#' @param Y A column name of a continuous variable.
#' @param feature_rank A character variable of the feature column.
#' @param logfold.threshold A Log2(A/B) Fold Change threshold, default is 0.6.
#' @param pvalue.threshold A P-value threshold, default is 0.05.
#' @return A \link[ggplot2]{ggplot2} object to be further modified.
#'
#' @export

volcano_plot <- function(dt, X, Y, feature_rank, pvalue.threshold = 0.05, logfold.threshold = 0.6, label_A = "A", label_B = "B") {
  # copies data.table
  tmpdt <- data.table::copy(dt)

  # Creates labels for significant and non-significant differential expression
  tmpdt[, (Y) := -log10(base::get(Y))]
  tmpdt[, diffexpressed := ifelse(base::get(X) > logfold.threshold & base::get(Y) > -log10(pvalue.threshold), "Upregulated",
                                   ifelse(base::get(X) < -logfold.threshold & base::get(Y) > -log10(pvalue.threshold), "Downregulated", "non-significant"))]
  tmpdt[, diffexpressed_labels := ifelse(diffexpressed != "non-significant", base::get(feature_rank), "")]

  max_pvalue <- max(tmpdt[[ Y ]])

  return(
    tmpdt %>%
      ggplot(mapping = aes(x = .data [[ X ]],
                           y = .data [[ Y ]],
                           label = diffexpressed_labels,
                           color = .data [[ X ]])) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
            axis.text.y = element_text(size=12),
            axis.text = element_text(size=12),
            text = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
      # geom_rect(aes(xmin = -logfold.threshold, xmax = -Inf,
      #               ymin = -log10(pvalue.threshold), ymax = Inf),
      #           fill = "#C8E7F1", alpha = 0.1, color = NA) +
      # geom_rect(aes(xmin = logfold.threshold, xmax = Inf,
      #               ymin = -log10(pvalue.threshold), ymax = Inf),
      #           fill = "#FFF1F3", alpha = 0.1, color = NA) +
      # annotate("text", x = -logfold.threshold*4 , y = max_pvalue+1.5,
      #          label = "Significant\ndecrease",
      #          vjust = 2, size = 5, color = "black") +
      # annotate("text", x = logfold.threshold*4 , y = max_pvalue+1.5,
      #          label = "Significant\nincrease",
      #          vjust = 2, size = 5, color = "black") +
      geom_vline(xintercept = c(-logfold.threshold, logfold.threshold),
                 col = "black", linetype = 'dashed') +
      geom_hline(yintercept = -log10(pvalue.threshold),
                 col = "black", linetype = 'dashed') +
      scale_color_gradient2(name = "foldchange",
                            low = "blue",
                            mid = "black",
                            high = "red",
                            na.value = "grey80") +
      ggrepel::geom_label_repel(show.legend = FALSE,
                                max.overlaps = getOption("ggrepel.max.overlaps", default = Inf),
                                color = "black") +
      geom_point(aes(size = ifelse(diffexpressed != "non-significant", rel_abun, 0)),
                 shape = 16, alpha = 0.5) +
      scale_size_continuous(name = "Mean Rel. Abun.") +
      labs(x = paste0("Fold Change log2( ",label_A," / ",label_B," )"),
           y = "-log10( P-value )")
    # +
    #   ylim(0, max_pvalue + 1.5)
  )
}
