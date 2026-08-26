#' Color map of a variable
#'
#' @description Creates an object of hexcode colors with names given a vector of characters.
#' This function is built into the \code{composition} method from the abstract class \link{omics} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#'
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table}.
#' @param groups A column name of a categorical variable in `data`.
#' @param col_name `r lifecycle::badge("deprecated")` This argument has been renamed to `groups` for more clarity.
#' @param Brewer.palID A character name that exists in \link[RColorBrewer]{brewer.pal} (Default: \code{"Set2"}).
#' @return A \link[stats]{setNames}.
#'
#' @examples 
#' library("data.table")
#' dt <- data.table(
#'   "SAMPLE_ID" = c("sample_1", "sample_2", "sample_3"),
#'   "treatment" = c("healthy", "tumor", NA)
#' )
#'
#' colors <- colormap(data = dt,
#'                    groups = "treatment")
#' @export

colormap <- function(
  data,
  col_name = lifecycle::deprecated(),
  groups,
  Brewer.palID = "Set2"
) {
  
  ## Lifecycle warn
  if (lifecycle::is_present(col_name)) {
    lifecycle::deprecate_warn(
      when = "1.7.0",
      what = "colormap(col_name)",
      with = "colormap(groups)"
    )
    if (missing(groups)) {
      groups <- col_name
    }
  }

  ## Error handling
  #--------------------------------------------------------------------#
  OPTIONS <- rownames(RColorBrewer::brewer.pal.info)

  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("{.val data} must be a {.cls data.frame} or {.cls data.table}.")

  if (!is.character(groups) || length(groups) != 1) {
    cli::cli_abort("{.val groups} needs to contain characters with length of 1.")
  } else if (!column_exists(groups, data)) {
    cli::cli_abort("The {.val groups} column does not exist in the provided data.")
  }

  if (!is.character(Brewer.palID) || length(Brewer.palID) != 1) {
    cli::cli_abort("The {.val Brewer.palID} needs to contain characters with length of 1.")
  } else if (!c(Brewer.palID %in% OPTIONS)) {
    cli::cli_abort("{.val Brewer.palID} is not a valid Brewer pal ID. \nValid options: {.val {OPTIONS}}.")
  }

  ## MAIN
  #--------------------------------------------------------------------#

  unique_groups <- unique(data[[groups]])
  len_unique_groups <- length(unique_groups)
  suppressWarnings(
    chosen_palette <- RColorBrewer::brewer.pal(len_unique_groups, Brewer.palID)
  )
  colors <- stats::setNames(chosen_palette[1:len_unique_groups], unique_groups)
  return(colors)
}
