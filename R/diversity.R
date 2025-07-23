#' Alpha Diversity
#'
#' @description Computes the alpha diversity based on Shannon index, simpson or invsimpson.
#' Code is adapted from \link[vegan]{diversity} and uses \link[Matrix]{sparseMatrix} in triplet format over the dense matrix.
#' The code is much faster and memory efficient, while still being mathematical correct.
#' This function is built into the \code{alpha_diversity} method from the abstract class \link[OmicFlow]{omics} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metagenomics} and \link[OmicFlow]{proteomics}.
#'
#' @param x A \link[base]{matrix} or \link[Matrix]{sparseMatrix}.
#' @param method A character variable for method; shannon, simpson or invsimpson.
#' @param normalize A boolean variable for sample normalization by column sums.
#' @param base Input for \link[base]{log} to use natural logarithmic scale, log2, log10 or other.
#' @return A numeric vector with type double.
#' @seealso \link[vegan]{diversity}
#'
#' @export

diversity <- function(x,
                      method = c("shannon", "simpson", "invsimpson"),
                      normalize = TRUE,
                      base = exp(1)) {

  ## Error handling
  #--------------------------------------------------------------------#

  x <- drop(as(x, "sparseMatrix"))
  if (!is.numeric(x@x))
    cli::cli_abort("input data must be numeric")

  if (any(x@x < 0, na.rm = TRUE))
    cli::cli_abort("input data must be non-negative")

  OPTIONS <- c("shannon", "simpson", "invsimpson")
  if (!is.character(method) && length(method) != 1) {
    cli::cli_abort("{method} needs to contain characters with length of 1.")
  } else if (!method %in% OPTIONS) {
    cli::cli_abort("{method} is not a valid method. Valid options: {OPTIONS}")
  }

  ## MAIN
  #--------------------------------------------------------------------#

  if (normalize) {
    total <- rep(Matrix::colSums(x), base::diff(x@p))
    x@x <- x@x / total
  }

  if (index == "shannon") {
    x@x <- -x@x * log(x@x, base)
  } else {
    x@x <- x@x * x@x
  }
  if (length(dim(x)) > 1) {
    H <- Matrix::colSums(x, na.rm = TRUE)
  }
  if (index == "simpson") {
    H <- 1 - H
  } else if (index == "invsimpson") {
    H <- 1/H
  }
  ## check NA in data
  if (any(NAS <- is.na(total)))
    H[NAS] <- NA
  H
}
