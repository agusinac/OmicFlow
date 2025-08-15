#' Ordination plot
#'
#' @description Creates an ordination plot pre-computed principal components from \link[vegan]{wcmdscale}. 
#' This function is built into the class \link{omics} with method \code{ordination()} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#' 
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table} of Principal Components as columns and rows as loading scores.
#' @param col_name A categorical variable to color the contrasts (e.g. "groups").
#' @param pair A vector of character variables indicating what dimension names (e.g. PC1, NMDS2).
#' @param dist_explained A vector of numeric values of the percentage dissimilarity explained for the dimension pairs, default is NULL.
#' @param dist_metric A character variable indicating what metric is used (e.g. unifrac, bray-curtis), default is NULL.
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#' 
#' @importFrom ggplot2 ggplot aes .data theme_bw theme element_text stat_ellipse labs geom_point
#' 
#' @examples 
#' library(ggplot2)
#' 
#' # Mock principal component scores
#' set.seed(123)
#' mock_data <- data.frame(
#'   SampleID = paste0("Sample", 1:10),
#'   PC1 = rnorm(10, mean = 0, sd = 1),
#'   PC2 = rnorm(10, mean = 0, sd = 1),
#'   groups = rep(c("Group1", "Group2"), each = 5)
#' )
#' 
#' # Basic usage
#' ordination_plot(
#'   data = mock_data,
#'   col_name = "groups",
#'   pair = c("PC1", "PC2")
#' )
#' 
#' # Adding variance/dissimilarity explained.
#' ordination_plot(
#'   data = mock_data,
#'   col_name = "groups",
#'   pair = c("PC1", "PC2"),
#'   dist_explained = c(45, 22),
#'   dist_metric = "bray-curtis"
#' )
#' @export

ordination_plot <- function(data, 
                            col_name,
                            pair, 
                            dist_explained = NULL, 
                            dist_metric = NULL) {
  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("data must be a data.frame or data.table.")

  if (!is.character(pair) && length(pair) != 2)
    cli::cli_abort("{pair} needs to contain characters with length of 2.")

  if (!is.character(col_name) && length(col_name) != 1) {
    cli::cli_abort("{col_name} needs to contain characters with length of 1.")
  } else if (!column_exists(col_name, data)) {
    cli::cli_abort("The {col_name} column does not exist in the provided data.")
  }

  if (!is.null(dist_explained) && !is.numeric(dist_explained))
    cli::cli_abort("{dist_explained} needs to be a numeric vector.")

  if (!is.null(dist_metric) && !is.character(dist_metric) && length(dist_metric) != 1)
    cli::cli_abort("{dist_metric} needs to contain characters with length of 1.")

  ## MAIN
  #--------------------------------------------------------------------#

  if (!is.null(dist_metric)) {
    plot_title = paste0("Distance metric used: ", dist_metric)
  } else {
    plot_title <- NULL
  }
  
  if (!is.null(dist_explained)) {
    x_label = paste0(pair[1], " (", round(as.numeric(dist_explained[1]), 2), "%)")
    y_label = paste0(pair[2], " (", round(as.numeric(dist_explained[2]), 2), "%)")
  } else {
    x_label = paste0(pair[1])
    y_label = paste0(pair[2])
  }

  data[[ col_name ]] <- as.factor(data[[ col_name ]])

  return(
    data %>%
      ggplot(mapping = aes(x = .data[[ pair[1] ]],
                           y = .data[[ pair[2] ]],
                           color = .data[[ col_name ]],
                           linetype = .data[[ col_name ]])) +
      geom_point(alpha = 5) +
      stat_ellipse(type = "t") +
      theme_bw() +
      theme(text=element_text(size=14),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            axis.text = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.text.x = element_text(size=12)
      ) +
      labs(title = NULL,
           subtitle = NULL,
           x = x_label,
           y = y_label
      )
  )
}
