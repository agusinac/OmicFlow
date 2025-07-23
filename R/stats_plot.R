#' Stats plot with ggplot2
#'
#' @description Creates a Stats plot of \link[OmicFlow]{pairwise_adonis} or \link[OmicFlow]{pairwise_anosim} results. This function is built into the \code{ordination} method from the abstract class \link[OmicFlow]{omics} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metagenomics} and \link[OmicsFLow]{proteomics}.
#'
#' @param data A \code{data.frame} or \code{data.table}.
#' @param Y A column name of a continuous variable.
#' @param X A column name of a categorical variable.
#' @param Label A column name of a categorical variable to label the bars.
#' @param Y_title A character variable to name the Y - axis title (default: NULL).
#' @param plot.title A character variable to name the plot title (default: NULL).
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#'
#' @export

stats_plot <- function(data,
                       X,
                       Y,
                       Label,
                       Y_title=NULL,
                       plot.title=NULL) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") || !inherits(data, "data.table"))
    cli::cli_abort("Data must be a data.frame or data.table.")

  if (!is.character(X) && length(X) != 1) {
    cli::cli_abort("{X} needs to contain characters with length of 1.")
  } else if (!column_exists(X, data)) {
    cli::cli_abort("The {X} column does not exist in the provided data.")
  }

  if (!is.character(Y) && length(Y) != 1) {
    cli::cli_abort("{Y} needs to contain characters with length of 1.")
  } else if (!column_exists(Y, data)) {
    cli::cli_abort("The {Y} column does not exist in the provided data.")
  }

  if (!is.character(Label))
    cli::cli_abort("{Label} needs to contain characters.")

  if (!is.null(Y_title) && !is.character(Y_title))
    cli::cli_abort("{Y_title} needs to contain characters.")

  if (!is.null(plot.title) && !is.character(plot.title))
    cli::cli_abort("{plot.title} needs to contain characters.")

  ## MAIN
  #--------------------------------------------------------------------#

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
