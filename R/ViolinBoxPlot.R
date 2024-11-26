#' Volin & Box plot with ggplot2
#'
#' @description Creates a \link[gghalves]{geom_half_violin} and \link[gghalves]{geom_half_boxplot} plot. This function is built into the \code{differential_feature_expression} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param dt A \code{data.frame} or \code{data.table}.
#' @param X A column name of a categorical variable.
#' @param Y A column name of a continuous variable.
#' @param feature_rank A character variable of the feature column.
#' @param logfold.threshold A Log2(A/B) Fold Change threshold, default is 0.0.
#' @param pvalue.threshold A P-value threshold, default is 0.05.
#' @return A \link[ggplot2]{ggplot} object to be further modified
#'
#' @export

ViolinBoxPlot <- function(dt, X, Y, diff, feature_rank, logfold.threshold = 0.6, pvalue.threshold = 0.05) {
  # Subset significant groups
  significant_subset <- dt[dt[, base::get(Y) < pvalue.threshold & (base::get(X) > logfold.threshold  | base::get(X) < -logfold.threshold)], ]
  significant_final <- dt[dt[[feature_rank]] %in% significant_subset[[feature_rank]], ]

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
