#' Compute Canberra Dissimilarity from a Sparse Matrix.
#'
#' @description 
#' Calculates the Canberra dissimilarity of a \link[Matrix]{sparseMatrix} pairwise for each column.
#' 
#' @details
#' The Canberra dissimilarity between two samples \eqn{A} and \eqn{B}, each of length \eqn{n}, is defined as:
#'
#' \eqn{d(A,B) = \frac{1 / NZ} \sum_{i}^n \frac{|A_i - B_i|}{|A_i| + |B_i|}}
#'
#' where \eqn{A_i} and \eqn{B_i} are the abundances of the \eqn{i}-th feature in sample \eqn{A} and \eqn{B}, respectively. NZ are the number of non-zero entries.
#' When weighted is set to FALSE, counts are replaced by presence/absence data.
#'
#' @param x A \link[Matrix]{sparseMatrix}.
#' @param weighted A boolean value, to use abundances (\code{weighted = TRUE}) or absence/presence (\code{weighted=FALSE}) (default: TRUE).
#' @param threads A wholenumber, the number of threads to use in \link[RcppParallel]{setThreadOptions} (default: 1).
#' @return A column x column \link[stats]{dist} object.
#' @references
#' Lance, G.N. & Williams, W.T. (1967) Mixed-data classificatory programs. I. Agglomerative systems. Australian Computer Journal, 1(1), 15-20.
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
#' canberra(taxa$countData)
#' @importFrom RcppParallel setThreadOptions
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist
#' @export

canberra <- function(x, weighted = TRUE, threads = 1) {

    ## Error handling
    #--------------------------------------------------------------------#
    if (is.vector(x))
        cli::cli_abort("Input must a matrix of class matrix or Matrix, not a vector.")

    x <- drop(as(x, "sparseMatrix"))
    if (!is.numeric(x@x))
        cli::cli_abort("Input data must be numeric.")

    if (!is.wholenumber(threads))
        cli::cli_abort("{threads} must be a whole number.")

    ## MAIN
    #--------------------------------------------------------------------#

    if (!weighted) x@x[] <- 1

    RcppParallel::setThreadOptions(numThreads = threads)
    out <- .Call('_OmicFlow_canberra', PACKAGE = 'OmicFlow', x)

    col_names <- colnames(x)
    if (!is.null(col_names))
        dimnames(out) <- list(col_names, col_names)

    return(as.dist(out))
}