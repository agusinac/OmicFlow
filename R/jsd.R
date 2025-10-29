#' Compute Jensen-Shannon Divergence from a Sparse Matrix.
#'
#' @description
#' Calculates the Jensen-Shannon divergence of a \link[Matrix]{sparseMatrix} pairwise for each column.
#' 
#' @details
#' The Jensen-Shannon divergence between two probability distributions \eqn{A} and \eqn{B}, each of length \eqn{n}, is defined as:
#'
#' \eqn{ d(A, B) = \frac{1}{2} D_{KL}(A \parallel M) + \frac{1}{2} D_{KL}(B \parallel M) }
#'
#' where \eqn{M = \frac{1}{2} (A + B)} is the mixture distribution,
#' and \eqn{D_{KL}} is the Kullback-Leibler divergence.
#' When weighted is set to FALSE, counts are replaced by presence/absence data.
#'
#' @param x A \link[Matrix]{sparseMatrix}.
#' @param weighted A boolean value, to use abundances (\code{weighted = TRUE}) or absence/presence (\code{weighted=FALSE}) (default: TRUE).
#' @param threads A wholenumber, the number of threads to use in \link[RcppParallel]{setThreadOptions} (default: 1).
#' @return A column x column \link[stats]{dist} object.
#' @references
#' Lin, J. (1991). Divergence measures based on the Shannon entropy. IEEE Transactions on Information Theory, 37(1), 145-151.
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
#' jsd(taxa$countData)
#' @importFrom RcppParallel setThreadOptions
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist
#' @export

jsd <- function(x, weighted = TRUE, threads = 1) {

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
    out <- .Call('_OmicFlow_jsd', PACKAGE = 'OmicFlow', x)

    col_names <- colnames(x)
    if (!is.null(col_names))
        dimnames(out) <- list(col_names, col_names)

    return(as.dist(out))
}