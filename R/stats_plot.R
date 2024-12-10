#' Stats plot with ggplot2
#'
#' @description Creates a Stats plot of \link[OmicFlow]{pairwise_adonis} or \link[OmicFlow]{pairwise_anosim} results. This function is built into the \code{ordination} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param data A \code{data.frame} or \code{data.table}.
#' @param Y A column name of a continuous variable.
#' @param X A column name of a categorical variable.
#' @param Label A column name of a categorical variable to label the bars.
#' @param plot.title A character variable to name the plot title, default is \code{NULL}.
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#'
#' @export

stats_plot <- function(data, X, Y, Label, Y_title=NULL, plot.title=NULL) {
  return(
    data %>%
      ggplot(mapping=aes(x = base::get(X, data),
                         y = base::get(Y, data),
                         label = base::get(Label, data))) +
      geom_bar(stat = "identity",
               fill = "blue") +
      geom_label(nudge_y = 0) +
      labs(title = plot.title,
           subtitle = "Above each bar: P-adjusted Values",
           x = "groups",
           y = Y_title) +
      theme_bw() +
      theme(axis.text = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1))
  )
}
