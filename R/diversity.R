#' Sparse implementation of Alpha Diversity Metrics
#'
#' @description Computes the alpha diversity based on Shannon index, simpson or invsimpson.
#' Code is adapted from \link[vegan]{diversity} and uses \link[Matrix]{sparseMatrix} in triplet format over the dense matrix.
#' The code is much faster and memory efficient, while still being mathematical correct.
#' This function is built into the class \link{omics} with method \code{alpha_diversity()} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#'
#' @param x A \link[base]{matrix} or \link[Matrix]{sparseMatrix}.
#' @param metric A character variable for metric; shannon, simpson or invsimpson.
#' @param normalize A boolean variable for sample normalization by column sums.
#' @param base Input for \link[base]{log} to use natural logarithmic scale, log2, log10 or other.
#' @return A numeric vector with type double.
#' @seealso \link[vegan]{diversity}
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
#' # Alpha diversity is computed on column level
#' ## Transpose the sparseMatrix if required with t() from Matrix R package.
#' result <- OmicFlow::diversity(
#'   x = sparse_mat,
#'   metric = "shannon"
#' ) 
#' @export

diversity <- function(x,
                      metric = c("shannon", "simpson", "invsimpson"),
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
  if (!is.character(metric) && length(metric) != 1) {
    cli::cli_abort("{metric} needs to contain characters with length of 1.")
  } else if (!metric %in% OPTIONS) {
    cli::cli_abort("{metric} is not a valid metric. Valid options: {OPTIONS}")
  }

  ## MAIN
  #--------------------------------------------------------------------#

  total <- rep(Matrix::colSums(x), base::diff(x@p))
  if (normalize) {
    x@x <- x@x / total
  }

  if (metric == "shannon") {
    x@x <- -x@x * log(x@x, base)
  } else {
    x@x <- x@x * x@x
  }
  if (length(dim(x)) > 1) {
    H <- Matrix::colSums(x, na.rm = TRUE)
  }
  if (metric == "simpson") {
    H <- 1 - H
  } else if (metric == "invsimpson") {
    H <- 1/H
  }
  ## check NA in data
  if (any(NAS <- is.na(total)))
    H[NAS] <- NA
  H
}
