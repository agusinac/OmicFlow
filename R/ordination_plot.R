#' Ordination plot
#'
#' @description Creates an ordination plot pre-computed principal components from \link[vegan]{wcmdscale}. 
#' This function is built into the class \link{omics} with method \code{ordination()} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#' 
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table} of Principal Components as columns and rows as loading scores.
#' @param groups A categorical variable to color the groups (e.g. \code{"treatment"}).
#' @param col_name `r lifecycle::badge("deprecated")` This argument has been renamed to `groups` for more clarity.
#' @param pair A vector of character variables indicating what dimension names (e.g. \code{c("PC1", "PC2")} or \code{c("NMDS1", "NMDS2")}).
#' @param dist_explained A vector of numeric values of the percentage dissimilarity explained for the dimension pairs (default: \code{NULL}).
#' @param dist_metric A character variable indicating what metric is used (e.g. unifrac, bray-curtis) (default: \code{NULL}).
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#' 
#' @examples 
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
#'   groups = "groups",
#'   pair = c("PC1", "PC2")
#' )
#' 
#' # Adding variance/dissimilarity explained.
#' ordination_plot(
#'   data = mock_data,
#'   groups = "groups",
#'   pair = c("PC1", "PC2"),
#'   dist_explained = c(45, 22),
#'   dist_metric = "bray-curtis"
#' )
#' @export

ordination_plot <- function(
  data, 
  col_name = lifecycle::deprecated(),
  groups = lifecycle::deprecated(),
  pair, 
  dist_explained = NULL, 
  dist_metric = NULL
  ) {

  ## lifecycle warn
  if (lifecycle::is_present(col_name)) {
    lifecycle::deprecate_warn(
      when = "1.7.0",
      what = "ordination_plot(col_name)",
      with = "ordination_plot(groups)"
    )
    if (!lifecycle::is_present(groups)) {
      groups <- col_name
    }
  }
  
  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("{.val data} must be a {.cls data.frame} or {.cls data.table}.")
  
  if (!is.character(groups) || length(groups) != 1) {
    cli::cli_abort("{.val groups} needs to contain characters with length of 1.")
  } else if (!column_exists(groups, data)) {
    cli::cli_abort("The {.val groups} column does not exist in the provided {.arg data}.")
  }

  if (!is.character(pair)) {
    cli::cli_abort("{.val pair} needs to be a characters {.cls vector}.")
  } else if (length(pair) != 2) {
    cli::cli_abort("{.val pair} needs to be a {.cls vector} of length 2.")
  }

  if (!is.null(dist_explained)) {
    if (!is.numeric(dist_explained)) {
      cli::cli_abort("{.val dist_explained} needs to be a numeric {.cls vector}.")
    } else if (length(dist_explained) != 2) {
      cli::cli_abort("{.val dist_explained} needs to be a {.cls vector} of length 2.")
    }
  }

  if (!is.null(dist_metric)) {
    if (!is.character(dist_metric) || length(dist_metric) != 1) {
      cli::cli_abort("{.val dist_metric} needs to contain characters with length of 1.")
    }
  }    

  ## MAIN
  #--------------------------------------------------------------------#

  if (!is.null(dist_metric))
    dist_metric <- paste0("Distance metric used: ", dist_metric)
  
  if (!is.null(dist_explained)) {
    x_label <- paste0(pair[1], " (", round(as.numeric(dist_explained[1]), 2), "%)")
    y_label <- paste0(pair[2], " (", round(as.numeric(dist_explained[2]), 2), "%)")
  } else {
    x_label <- paste0(pair[1])
    y_label <- paste0(pair[2])
  }

  data[[ groups ]] <- as.factor(data[[ groups ]])

  plt <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
      x = .data[[ pair[1] ]],
      y = .data[[ pair[2] ]],
      color = .data[[ groups ]],
      linetype = .data[[ groups ]]
    )
  ) +
  ggplot2::geom_point(alpha = 5) +
  ggplot2::stat_ellipse(type = "t") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size=14),
    legend.text = ggplot2::element_text(size=12),
    legend.title = ggplot2::element_text(size=14),
    axis.text = ggplot2::element_text(size=12),
    axis.text.y = ggplot2::element_text(size=12),
    axis.text.x = ggplot2::element_text(size=12)
  )

  if (length(unique(data[[groups]])) <= 8) {
    plt <- plt +
    ggplot2::scale_colour_manual(
      name = groups,
      values = colormap(data = data, groups = groups)
    )
  }
  plt <- plt +
  ggplot2::labs(
    title = dist_metric,
    subtitle = NULL,
    x = x_label,
    y = y_label
  )
  return(plt)
}
