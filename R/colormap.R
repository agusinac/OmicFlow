#' Color map of a variable
#'
#' @description Creates an object of hexcode colors with names given a vector of characters.
#' This function is built into the \code{ordination} method from the abstract class \link{omics} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#'
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table}.
#' @param col_name A column name of a categorical variable.
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
#'                    col_name = "treatment")
#' @export

colormap <- function(data,
                     col_name,
                     Brewer.palID = "Set2") {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("Data must be a data.frame or data.table.")

  if (!is.character(col_name) && length(col_name) != 1) {
    cli::cli_abort("Column name: {col_name} needs to contain characters with length of 1.")
  } else if (!column_exists(col_name, data)) {
    cli::cli_abort("The {col_name} column does not exist in the provided data.")
  }

  if (!is.character(Brewer.palID) && length(Brewer.palID) != 1)
    cli::cli_abort("The {Brewer.palID} needs to contain characters with length of 1.")

  ## MAIN
  #--------------------------------------------------------------------#

  unique_groups <- unique(data[[col_name]])
  chosen_palette <- RColorBrewer::brewer.pal(length(unique_groups), Brewer.palID)
  colors <- stats::setNames(chosen_palette, unique_groups)
  return(colors)
}
