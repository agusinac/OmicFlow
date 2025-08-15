#' Sparse implementation of Hill numbers
#'
#' @description Computes the hill numbers for q is 0, 1 or 2.
#' Code is adapted from \link[hillR]{hill_taxa} and uses \link[Matrix]{sparseMatrix} in triplet format over the dense matrix.
#' The code is much faster and memory efficient, while still being mathematical correct.
#'
#' @param x A \link[base]{matrix} or \link[Matrix]{sparseMatrix}.
#' @param q A wholenumber for 0, 1 or 2, default is 0.
#' @param normalize A boolean variable for sample normalization by column sums.
#' @param base Input for \link[base]{log} to use natural logarithmic scale, log2, log10 or other.
#' @return A numeric vector with type double.
#' @seealso \link[hillR]{hill_taxa}
#' 
#' @importFrom methods as
#' 
#' @examples 
#' library("Matrix")
#' 
#' n_row <- 1000
#' n_col <- 100
#' density <- 0.2
#' num_entries <- n_row * n_col
#' num_nonzero <- round(num_entries * density)
#'
#' set.seed(123)
#' positions <- sample(num_entries, num_nonzero, replace=FALSE)
#' row_idx <- ((positions - 1) %% n_row) + 1
#' col_idx <- ((positions - 1) %/% n_row) + 1
#'
#' values <- runif(num_nonzero, min = 0, max = 1)
#' sparse_mat <- sparseMatrix(
#'   i = row_idx,
#'   j = col_idx,
#'   x = values,
#'   dims = c(n_row, n_col)
#' )
#' 
#' result <- OmicFlow::hill_taxa(
#'  x = sparse_mat,
#'  q = 2
#' )
#' @export

hill_taxa <- function (x,
                       q = 0,
                       normalize = TRUE, 
                       base = exp(1)) {

  ## Error handling
  #--------------------------------------------------------------------#

  x <- drop(as(x, "sparseMatrix"))
  if (!is.numeric(x@x))
    cli::cli_abort("input data must be numeric")

  if (!is.wholenumber(q) && q %in% c(0, 1, 2)) 
    cli::cli_abort("{q} needs to be a whole number and either 0, 1 or 2.")

  if (any(x@x < 0, na.rm = TRUE))
    cli::cli_abort("input data must be non-negative")

  ## MAIN
  #--------------------------------------------------------------------#

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
