#' Fetch palette(s)
#'
#' @description Creates an object of hexcode colors with names.
#' This function is built into the \code{ordination} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param metadata A \code{data.frame} or \code{data.table}.
#' @param col_name A column name that exists in metadata.
#' @param Brewer.palID A character name that exists in \link[RColorBrewer]{brewer.pal}, default: \code{"Set2"}.
#' @return An object of names and hexcolors
#' @seealso \link[stats]{setNames}
#'
#' @export

fetch_palette <- function(metadata, col_name, Brewer.palID="Set2") {
  # Creating color palette
  unique_groups <- unique(metadata[[col_name]])
  chosen_palette <- RColorBrewer::brewer.pal(length(unique_groups), Brewer.palID)
  colors <- stats::setNames(chosen_palette, unique_groups)
  return(colors)
}
