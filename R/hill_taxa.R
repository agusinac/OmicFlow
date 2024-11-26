#' Hill numbers
#'
#' @description Computes the hill numbers for q is 0, 1 or 2.
#' Code is adapted from \link[hillR]{hill_taxa} and uses \link[Matrix]{sparseMatrix} in triplet format over the dense matrix.
#' The code is much faster and memory efficient, while still being mathematical correct.
#' This function is built into the \code{alpha_diversity} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param x A \code{matrix} or \code{sparseMatrix}.
#' @param q An integer for 0, 1 or 2, default is 0.
#' @param normalize A boolean variable for sample normalization by column sums.
#' @param base Input for \link[base]{log} to use natural logarithmic scale, log2, log10 or other.
#' @return A numeric vector with type double.
#' @seealso \link[hillR]{hill_taxa}
#'
#' @export
hill_taxa <- function (x, q = 0, normalize = TRUE, base = exp(1)) {
  x <- drop(as(x, "sparseMatrix"))
  if (!is.numeric(x@x))
    stop("input data must be numeric")
  if (any(x@x < 0, na.rm = TRUE))
    stop("input data must be non-negative")
  if (normalize) {
    total <- rep(Matrix::colSums(x), base::diff(x@p))
    x@x <- x@x / total
  }
  if (q == 0) {
    hill <- base::diff(x@p)
  } else {
    if (q == 1) {
      x@x <- -x@x * log(x@x, base)
      hill <- exp(Matrix::colSums(x, na.rm = TRUE))
    } else {
      x@x <- x@x^q
      total <- Matrix::colSums(x, na.rm = TRUE)
      hill <- total^(1/(1 - q))
    }
  }
  hill
}
