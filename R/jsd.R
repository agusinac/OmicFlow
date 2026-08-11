#' Compute Jensen-Shannon Divergence from a Dense or Sparse Matrix.
#'
#' @description
#' Calculates the Jensen-Shannon divergence of a Matrix pairwise for each column.
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
#' @inheritParams bray
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
#' taxa$scale(method = "tss")
#'
#' jsd(taxa$countData)
#' @export

jsd <- function(x, weighted = TRUE, threads = 1) {

    ## Error handling
    #--------------------------------------------------------------------#
    if (is.vector(x))
        cli::cli_abort("{.val x} must be a {.cls matrix}, {.cls denseMatrix} or {.cls sparseMatrix}, not a {.cls vector}.")
    
    if (inherits(x, "denseMatrix") || inherits(x, "matrix") || inherits(x, "sparseMatrix")) {
        x <- methods::as(x, "CsparseMatrix")
    } else cli::cli_abort("{.val x} isn't a {.cls matrix}, {.cls denseMatrix} or {.cls sparseMatrix}.")
    
    if (!is.numeric(x@x))
        cli::cli_abort("{.val x} must be numeric.")

    if (any(x@x < 0, na.rm = TRUE))
        cli::cli_abort("{.val x} must be non-negative.")
    
    if (!is.logical(weighted))
        cli::cli_abort("{.val weighted} needs to be either `TRUE` or `FALSE`.")

    if (length(threads) != 1) {
        cli::cli_abort("{.val threads} must be a single whole number.")
    } else if (!is.wholenumber(threads)) {
        cli::cli_abort("{.val threads} must be a whole number.")
    }

    ## MAIN
    #--------------------------------------------------------------------#

    if (!weighted) x@x[] <- 1

    RcppParallel::setThreadOptions(numThreads = threads)
    out <- .Call('_OmicFlow_jsd', PACKAGE = 'OmicFlow', x)

    col_names <- colnames(x)
    if (!is.null(col_names))
        dimnames(out) <- list(col_names, col_names)

    return(stats::as.dist(out))
}