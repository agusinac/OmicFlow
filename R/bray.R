#' Compute Bray-Curtis Dissimilarity from a Dense or Sparse Matrix.
#'
#' @description 
#' Calculates the Bray-Curtis dissimilarity of a Matrix pairwise for each column.
#' 
#' @details
#' The Bray-Curtis dissimilarity between two samples \eqn{A} and \eqn{B}, each of length \eqn{n}, is defined as:
#'
#' \eqn{d(A,B) = \frac{\sum_{i}^n |A_i - B_i|}{\sum_{i}^n (A_i + B_i)}}
#'
#' where \eqn{A_i} and \eqn{B_i} are the abundances of the \eqn{i}-th feature in sample \eqn{A} and \eqn{B}, respectively. 
#' When weighted is set to FALSE, counts are replaced by presence/absence data.
#'
#' @param x A \link[base]{matrix}, \link[Matrix]{sparseMatrix} or \link[Matrix]{Matrix}.
#' @param weighted A boolean value, to use abundances (\code{weighted = TRUE}) or absence/presence (\code{weighted=FALSE}) (default: \code{TRUE}).
#' @param threads A wholenumber, the number of threads to use in \link[RcppParallel]{setThreadOptions} (default: \code{1}).
#' @return A column x column \link[stats]{dist} object.
#' @references
#' Bray, J.R. & Curtis, J.T. (1957) An Ordination of the Upland Forest Communities of Southern Wisconsin. Ecological Monographs, 27(4), 325–349.
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
#' bray(taxa$countData)
#' @importFrom Matrix sparseMatrix
#' @export

bray <- function(x, weighted = TRUE, threads = 1) {

    ## Error handling
    #--------------------------------------------------------------------#
    if (is.vector(x))
        cli::cli_abort("{.val x} must be a {.cls matrix}, {.cls denseMatrix} or {.cls sparseMatrix}, not a {.cls vector}.")
    
    if (inherits(x, "denseMatrix") || inherits(x, "matrix") || inherits(x, "sparseMatrix")) {
        x <- methods::as(x, "CsparseMatrix")
    } else cli::cli_abort("{.val x} isn't a {.cls matrix}, {.cls denseMatrix} or {.cls sparseMatrix}.")
    
    if (!is.numeric(x@x))
        cli::cli_abort("{.val x} must be numeric.")

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
    out <- .Call('_OmicFlow_bray', PACKAGE = 'OmicFlow', x)

    col_names <- colnames(x)
    if (!is.null(col_names))
        dimnames(out) <- list(col_names, col_names)

    return(stats::as.dist(out))
}
