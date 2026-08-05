#' Pairwise anosim (ANOSIM) computation
#'
#' @description Computes pairwise \link[vegan]{anosim}, given a distance matrix and a vector of labels.
#' This function is built into the class \link{omics} with method \code{ordination()} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#'
#' @inheritParams pairwise_adonis
#' @seealso \link[vegan]{anosim}
#' @return A \link[base]{data.frame} containing: \describe{
#' \item{pairs}{combinations of group comparisons}
#' \item{Df}{Degrees of Freedom}
#' \item{R2}{The R squared of the null hypothesis \eqn{H_0}}
#' \item{p.value}{The number of \eqn{F^p} higher than the null hypothesis divided by the total number of permutations.}
#' \item{p.adj}{The adjusted P-value based on the used `p.adjust.method`}
#' }
#' @examples 
#' # Create random data
#' set.seed(42)
#' mock_data <- matrix(rnorm(15 * 10), nrow = 15, ncol = 10)
#' 
#' # Create euclidean dissimilarity matrix
#' mock_dist <- dist(mock_data, method = "euclidean")
#' 
#' # Define group labels, should be equal to number of columns and rows to dist
#' mock_groups <- rep(c("A", "B", "C"), each = 5)
#' 
#' # Compute pairwise anosim
#' result <- pairwise_anosim(x = mock_dist, 
#'                           groups = mock_groups, 
#'                           p.adjust.method = "bonferroni", 
#'                           perm = 99)
#' @export

pairwise_anosim <- function(
  x,
  groups,
  metadata = NULL,
  perm_design = NULL,
  p.adjust.method = "bonferroni",
  perm = 999
  ){

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(x, "dist"))
    cli::cli_abort("{.val x} must be a {.cls dist}")

  if (is.list(groups))
    cli::cli_abort("{.val groups} must be a {.cls vector} and not a {.cls list}.")
  
  if (!is.null(metadata) && !inherits(metadata, "data.frame") && !inherits(metadata, "data.table"))
    cli::cli_abort("{.val metadata} must be a {.cls data.frame} or {.cls data.table}.")
  
  if (!is.null(perm_design) && !is.function(perm_design))
    cli::cli_abort("{.val perm_design} must be a function.")

  if (!is.character(p.adjust.method)) {
    cli::cli_abort("{.val p.adjust.method} must be a character.")
  } else if (!c(p.adjust.method %in% stats::p.adjust.methods)) {
    cli::cli_abort("{.val {p.adjust.method}} is not a valid method. \nValid options: {.val {p.adjust.methods}}.")
  }

  if (length(perm) != 1) {
    cli::cli_abort("{.val perm} must be a single whole number.")
  } else if (!is.wholenumber(perm)) {
    cli::cli_abort("{.val perm} needs to be a whole number.")
  }

  ## MAIN
  #--------------------------------------------------------------------#

  co <- utils::combn(unique(as.character(groups)), 2)
  n <- ncol(co)
  pairs <- vector(mode = "numeric", length = n)
  anosimR <- vector(mode = "numeric", length = n)
  p.value <- vector(mode = "numeric", length = n)

  for(i in 1:n){
    if(inherits(x, "dist")){
      rows_to_keep <- groups %in% co[, i]
      m <- as.matrix(x)[rows_to_keep, rows_to_keep]
    }

    # Apply permutation design
    if (!is.null(perm_design) && !is.null(metadata)) {
      sub_meta <- metadata[rows_to_keep, ]
      h1 <- perm_design(sub_meta)
      ano <- vegan::anosim(
        x = m,
        grouping = groups,
        permutations = h1
      )
    } else {
      ano <- vegan::anosim(
        x = m, 
        grouping = groups[rows_to_keep],
        permutations = perm
      )
    }
    
    pairs[i] <- paste(co[1, i],'vs',co[2, i])
    anosimR[i] <- ano$statistic
    p.value[i] <- ano$signif
  }
  p.adj <- stats::p.adjust(p.value, method = p.adjust.method)
  pairw.res <- data.frame(pairs, anosimR, p.value, p.adj)
  return(pairw.res)
}
