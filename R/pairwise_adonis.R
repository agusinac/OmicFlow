#' Pairwise adonis2 (PERMANOVA) computation
#'
#' @description Computes pairwise \link[vegan]{adonis2}, given a distance matrix and a vector of labels.
#' This function is built into the class \link{omics} with method \code{ordination()} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#'
#' @param x A \link[stats]{dist}.
#' Obtained from a dissimilarity metric, in the case of similarity metric please use \code{1-dist}.
#' @param groups A character vector (e.g. a column from a the `metadata`) to match the sample group labels.
#' @param metadata A \link[data.table]{data.table} or \link[base]{data.frame} as input to the function \code{perm_design} (default: \code{NULL}).
#' @param perm_design A function that takes the `metadata` and returns a permutation design with \link[permute]{how} (default: \code{NULL}).
#' @param p.adjust.method A character as input to adjust the p-values, see \link[stats]{p.adjust} (default: \code{"bonferroni"}).
#' @param perm A whole number to define the number of permutations in \link[vegan]{adonis2} (default: \code{999}).
#' @seealso \link[vegan]{adonis2}
#' @return A \link[base]{data.frame} containing: \describe{
#' \item{pairs}{combinations of group comparisons}
#' \item{Df}{Degrees of Freedom}
#' \item{SumsOfSqs}{The sums of squares (centroid) of the null hypothesis \eqn{H_0}}
#' \item{F.Model}{The F-test of the null hypothesis \eqn{H_0}}
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
#' # Compute pairwise adonis (PERMANOVA)
#' result <- pairwise_adonis(x = mock_dist, 
#'                           groups = mock_groups, 
#'                           p.adjust.method = "bonferroni", 
#'                           perm = 99)
#' @export

pairwise_adonis <- function(
  x,
  groups,
  metadata = NULL,
  perm_design = NULL,
  p.adjust.method = "bonferroni",
  perm = 999){

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
    cli::cli_abort("{.val {p.adjust.method}} must be a character.")
  } else if (!c(p.adjust.method %in% stats::p.adjust.methods)) {
    cli::cli_abort("{.val {p.adjust.method}} is not a valid method. \nValid options: {.val {p.adjust.methods}}.")
  }

  if (length(perm) != 1) {
    cli::cli_abort("{.val perm} must be a single whole number.")
  } else if (!is.wholenumber(perm)) {
    cli::cli_abort("{.val {perm}} needs to be a whole number.")
  }
    

  ## MAIN
  #--------------------------------------------------------------------#

  # Initialize required parameters
  co <- utils::combn(unique(as.character(groups)), 2)
  n <- ncol(co)
  pairs <- vector(mode = "numeric", length = n)
  p.value <- vector(mode = "numeric", length = n)
  Df <- vector(mode = "numeric", length = n)
  SumsOfSqs <- vector(mode = "numeric", length = n)
  F.Model <- vector(mode = "numeric", length = n)
  R2 <- vector(mode = "numeric", length = n)

  # Loops through pairs
  for(i in 1:n){
    if(inherits(x, 'dist')){
      rows_to_keep <- groups %in% co[, i]
      m <- as.matrix(x)[rows_to_keep, rows_to_keep]
    }
    tmp_m <- data.frame(Fac = groups[rows_to_keep])

    # Apply permutation design
    if (!is.null(perm_design) && !is.null(metadata)) {
      sub_meta <- metadata[rows_to_keep, ]
      h1 <- perm_design(sub_meta)
      ad <- vegan::adonis2(
        m ~ Fac,
        data = tmp_m,
        permutations = h1
      )
    } else {
      ad <- vegan::adonis2(
        m ~ Fac,
        data = tmp_m,
        permutations = perm
      )
    }

    # Saving stats
    pairs[i] <- paste(co[1, i],'vs',co[2, i])
    Df[i] <- ad$Df[1]
    SumsOfSqs[i] <- ad$SumOfSqs[1]
    F.Model[i] <- ad$F[1]
    R2[i] <- ad$R2[1]
    p.value[i] <- ad$`Pr(>F)`[1]
  }
  # Adjusts P-values and returns combined dataframe
  p.adj <- stats::p.adjust(p.value, method = p.adjust.method)
  pairw.res <- data.frame(pairs, Df, SumsOfSqs, F.Model, R2, p.value, p.adj)
  return(pairw.res)
}
