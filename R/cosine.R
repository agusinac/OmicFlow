#' Compute Cosine Dissimilarity from a Dense or Sparse Matrix.
#'
#' @description 
#' Calculates the cosine disimilarity of a Marix pairwise for each column.
#' 
#' @details
#' The cosine dissimilarity between two samples \eqn{A} and \eqn{B}, each of length \eqn{n}, is defined as:
#'
#' \eqn{d(A,B) = 1 - \frac{\sum_{i}^n A_i B_i}{\sqrt{\sum_{i}^n A_i^2} \sqrt{\sum_{i}^n B_i^2}} }
#'
#' where \eqn{A_i} and \eqn{B_i} are the abundances of the \eqn{i}-th feature in sample \eqn{A} and \eqn{B}, respectively.
#' When weighted is set to FALSE, counts are replaced by presence/absence data.
#'
#' @inheritParams bray
#' @return A column x column \link[stats]{dist} object.
#' @references
#' Deza, M. M., & Deza, E. (2009). Encyclopedia of Distances. Springer Science & Business Media., 308.
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
#' cosine(taxa$countData)
#' @importFrom Matrix sparseMatrix
#' @export

cosine <- function(x, weighted = TRUE, threads = 1) {

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
        cli::cli_abort("{.val {threads}} must be a whole number.")
    }

    ## MAIN
    #--------------------------------------------------------------------#

    if (!weighted) x@x[] <- 1

    RcppParallel::setThreadOptions(numThreads = threads)
    out <- .Call('_OmicFlow_cosine', PACKAGE = 'OmicFlow', x)

    col_names <- colnames(x)
    if (!is.null(col_names))
        dimnames(out) <- list(col_names, col_names)

    return(stats::as.dist(out))
}
