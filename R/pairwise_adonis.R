#' Pairwise adonis2 (PERMANOVA) computation
#'
#' @description Computes pairwise \link[vegan]{adonis2}, given a distance matrix and a vector of labels.
#' This function is built into the class \link{omics} with method \code{ordination()} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#'
#' @param x A distance matrix in the form of \link[stats]{dist}.
#' Obtained from a dissimilarity metric, in the case of similarity metric please use \code{1-dist}
#' @param groups A character vector (column from a table) of labels.
#' @param p.adjust.method P adjust method see \link[stats]{p.adjust}
#' @param perm  Number of permutations to compare against the null hypothesis of adonis2 (default: \code{perm=999}).
#' @seealso \link[vegan]{adonis2}
#' @return A \link[base]{data.frame} of
#'  * pairs that are used
#'  * Degrees of freedom (Df)
#'  * Sums of Squares of H_0
#'  * F.Model of H_0
#'  * R2 of H_0
#'  * p value of F^p > F
#'  * p adjusted
#' @importFrom stats p.adjust.methods
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

pairwise_adonis <- function(x,
                            groups,
                            p.adjust.method = "bonferroni",
                            perm = 999){

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(x, "dist"))
    cli::cli_abort("x must be of class dist.")

  if (!is.vector(groups))
    cli::cli_abort("{groups} must be a vector.")

  if (!c(p.adjust.method %in% p.adjust.methods))
    cli::cli_abort("Specified {p.adjust.method} is not valid. \nValid options: {p.adjust.methods}.")

  if (!is.wholenumber(perm))
    cli::cli_abort("Permutations {perm} need to be an integer.")

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
      m = as.matrix(x)[groups %in% c(as.character(co[1, i]),as.character(co[2, i])),
                      groups %in% c(as.character(co[1, i]),as.character(co[2, i]))]
    }
    # Performing adonis2 test
    tmp_m = data.frame(Fac = groups[groups %in% c(co[1,i], co[2,i])])
    ad <- vegan::adonis2(m ~ Fac,
                         data = tmp_m,
                         permutations = perm);

    # Saving stats
    pairs[i] <- paste(co[1, i],'vs',co[2, i])
    Df[i] <- ad$Df[1]
    SumsOfSqs[i] <- ad$SumOfSqs[1]
    F.Model[i] <- ad$F[1]
    R2[i] <- ad$R2[1]
    p.value[i] <- ad$`Pr(>F)`[1]
  }
  # Adjusts P-values and returns combined dataframe
  p.adj <- p.adjust(p.value, method = p.adjust.method)
  pairw.res <- data.frame(pairs, Df, SumsOfSqs, F.Model, R2, p.value, p.adj)
  return(pairw.res)
}
