#' Create pairwise stats plot
#'
#' @description Creates a pairwise stats plot from \link{pairwise_adonis} or \link{pairwise_anosim} results. 
#' This function is built into the class \link{omics} with method \code{ordination()} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#'
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table}.
#' @param stats_col A column name of a continuous variable.
#' @param group_col A column name of a categorical variable.
#' @param label_col A column name of a categorical variable to label the bars.
#' @param y_axis_title A character variable to name the Y - axis title (default: NULL).
#' @param plot_title A character variable to name the plot title (default: NULL).
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#' @importFrom ggplot2 ggplot aes geom_bar geom_label theme_bw theme element_text labs
#' @examples 
#' library("ggplot2")
#' 
#' # Create random data
#' set.seed(42)
#' mock_data <- matrix(rnorm(15 * 10), nrow = 15, ncol = 10)
#' 
#' # Create euclidean dissimilarity matrix
#' mock_dist <- dist(mock_data, method = "euclidean")
#' 
#' # Define group labels, should be equal to number of columns and rows to dist
#' mock_groups <- rep(c("A", "B", "C"), each = 5)
#' 
#' # Compute pairwise adonis
#' adonis_res <- pairwise_adonis(x = mock_dist, 
#'                               groups = mock_groups, 
#'                               p.adjust.method = "bonferroni", 
#'                               perm = 99)
#' # Compute pairwise anosim
#' anosim_res <- pairwise_anosim(x = mock_dist, 
#'                               groups = mock_groups, 
#'                               p.adjust.method = "bonferroni", 
#'                               perm = 99)
#' 
#' # Visualize PERMANOVA pairwise stats
#' plot_pairwise_stats(data = adonis_res,
#'                     group_col = "pairs",
#'                     stats_col = "F.Model",
#'                     label_col = "p.adj",
#'                     y_axis_title = "Pseudo F test statistic",
#'                     plot_title = "PERMANOVA")
#' 
#' # Visualize ANOSIM pairwise stats
#' plot_pairwise_stats(data = anosim_res,
#'                     group_col = "pairs",
#'                     stats_col = "anosimR",
#'                     label_col = "p.adj",
#'                     y_axis_title = "ANOSIM R statistic",
#'                     plot_title = "ANOSIM")
#' @export

plot_pairwise_stats <- function(data,
                                stats_col,
                                group_col,
                                label_col,
                                y_axis_title=NULL,
                                plot_title=NULL) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("Data must be a data.frame or data.table.")

  if (!is.character(stats_col) && length(stats_col) != 1) {
    cli::cli_abort("{stats_col} needs to contain characters with length of 1.")
  } else if (!column_exists(stats_col, data)) {
    cli::cli_abort("The {stats_col} column does not exist in the provided data.")
  }

  if (!is.character(group_col) && length(group_col) != 1) {
    cli::cli_abort("{group_col} needs to contain characters with length of 1.")
  } else if (!column_exists(group_col, data)) {
    cli::cli_abort("The {group_col} column does not exist in the provided data.")
  }

  if (!is.character(label_col))
    cli::cli_abort("{label_col} needs to contain characters.")

  if (!is.null(y_axis_title) && !is.character(y_axis_title))
    cli::cli_abort("{y_axis_title} needs to contain characters.")

  if (!is.null(plot_title) && !is.character(plot_title))
    cli::cli_abort("{plot_title} needs to contain characters.")

  ## MAIN
  #--------------------------------------------------------------------#

  return(
    data %>%
      ggplot(mapping=aes(x = base::get(group_col, data),
                         y = base::get(stats_col, data),
                         label = base::get(label_col, data))) +
      geom_bar(stat = "identity",
               fill = "blue") +
      geom_label(nudge_y = 0) +
      labs(title = plot_title,
           subtitle = paste0("Above each bar: ", label_col, " values"),
           x = "groups",
           y = y_axis_title) +
      theme_bw() +
      theme(axis.text = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1))
  )
}
