#' Compute Aitchison Dissimilarity from a Sparse Matrix.
#'
#' @description
#' Calculates the Aitchison dissimilarity of compositional data represented by a \link[Matrix]{sparseMatrix} pairwise for each column.
#' The Aitchison dissimilarity between two compositions \eqn{A} and \eqn{B}, each of length \eqn{n}, is defined as the Euclidean distance of their centered log-ratio (clr) transforms:
#'
#' \eqn{ d(A, B) = \sqrt{\sum_{i}^n \left[ \log \frac{A_i}{g(x)} - \log \frac{B_i}{g(y)} \right]^2 } }
#'
#' where \eqn{g(A)} and \eqn{g(B)} are the geometric means of components of \eqn{A} and \eqn{B}.
#' When weighted is set to FALSE, counts are replaced by presence/absence data.
#'
#' @param x A \link[Matrix]{sparseMatrix} of strictly positive values.
#' @param weighted A boolean value, to use counts or presence/absence (default: TRUE).
#' @param threads A wholenumber, the number of threads to use in \link[RcppParallel]{setThreadOptions} (default: 1).
#' @return A column x column \link[stats]{dist} object.
#' @references
#' Aitchison, J. (1986). The Statistical Analysis of Compositional Data. Journal of the Royal Statistical Society. Series B (Methodological), 44(2), 139-160.
#' @examples 
#' library("OmicFlow")
#'
#' metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
#' counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
#' features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
#' tree_file <- system.file("extdata", "tree.newick", package = "OmicFlow")
#'
#' taxa <- metagenomics$new(
#'     metaData = metadata_file,
#'     countData = counts_file,
#'     featureData = features_file,
#'     treeData = tree_file
#' )
#'
#' taxa$feature_subset(Kingdom == "Bacteria")
#' taxa$normalize()
#'
#' aitchison(taxa$countData)
#' @importFrom RcppParallel setThreadOptions
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist
#' @export

aitchison <- function(x, weighted = TRUE, threads = 1) {

    ## Error handling
    #--------------------------------------------------------------------#
    if (is.vector(x))
        cli::cli_abort("Input must a matrix of class matrix or Matrix, not a vector.")

    x <- drop(as(x, "sparseMatrix"))
    if (!is.numeric(x@x))
        cli::cli_abort("Input data must be numeric.")

    if (any(x@x < 0, na.rm = TRUE))
        cli::cli_abort("Input data must be non-negative.")

    if (!is.wholenumber(threads))
        cli::cli_abort("{threads} must be a whole number.")

    ## MAIN
    #--------------------------------------------------------------------#

    if (!weighted) x@x[] <- 1

    RcppParallel::setThreadOptions(numThreads = threads)
    out <- .Call('_OmicFlow_aitchison', PACKAGE = 'OmicFlow', x)

    col_names <- colnames(x)
    if (!is.null(col_names))
        dimnames(out) <- list(col_names, col_names)

    return(as.dist(out))
}