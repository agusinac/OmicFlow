# Adapted from vegan::diversity

diversity <- function(x, index = "shannon", normalize = TRUE, base = exp(1)) {
  x <- drop(as(x, "sparseMatrix"))
  if (!is.numeric(x@x))
    stop("input data must be numeric")
  if (any(x@x < 0, na.rm = TRUE))
    stop("input data must be non-negative")

  INDICES <- c("shannon", "simpson", "invsimpson")
  index <- match.arg(index, INDICES)
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
