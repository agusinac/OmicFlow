#' Compute UniFrac Dissimilarity from a Sparse Matrix.
#'
#' @description Calculates the UniFrac dissimilarity between samples based on phylogenetic branch lengths and abundance or presence/absence data.
#' 
#' @details
#' The UniFrac distance between two samples \eqn{A} and \eqn{B}, with phylogenetic tree edges \eqn{i = 1 \ldots n} of lengths \eqn{L_i}, is computed differently depending on the \code{weighted} and \code{normalized} flags. 
#' When \code{weighted = FALSE}, input counts are first converted to presence/absence data.
#' \describe{
#'  \item{Weighted UniFrac (\code{normalized = FALSE} and \code{weighted = TRUE}):}{
#'      \eqn{d(A,B) = \frac{\sum_{i}^n L_i |A_i - B_i|}{\sum_{i}^n L_i (A_i + B_i)}}
#'  }
#'  \item{Normalized Weighted UniFrac (\code{normalized = TRUE} and \code{weighted = TRUE}):}{
#'      \eqn{d(A,B) = \sum_{i}^n L_i |A_i - B_i|}
#'  }
#'  \item{Unweighted UniFrac (\code{weighted = FALSE}, unweighted is always normalized):}{
#'      \eqn{d(A,B) = \frac{\sum_{i}^n L_i |A_i - B_i|}{\sum_{i}^n L_i \max(A_i, B_i)}}
#'  }
#'}
#' @param x A \link[Matrix]{sparseMatrix} of strictly positive counts or presence/absence data.
#' @param tree A `phylo` class tree.
#' @param weighted A boolean value, to use abundances (\code{weighted = TRUE}) or absence/presence (\code{weighted=FALSE}) (default: TRUE).
#' @param normalized A boolean value, whether to normalize weighted UniFrac distances to be between 0 and 1 (default: TRUE). Unweighted UniFrac is always normalized.
#' @param threads A wholenumber, the number of threads to use in \link[RcppParallel]{setThreadOptions} (default: 1).
#' @return A column x column \link[stats]{dist} object.
#' @references
#' Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method for comparing microbial communities. Applied and Environmental Microbiology, 71(12), 8228–8235.
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
#' # Weighted UniFrac
#' unifrac(x = taxa$countData, tree = taxa$treeData, weighted=TRUE, normalized=FALSE)
#' 
#' # Weighted Normalized UniFrac
#' unifrac(x = taxa$countData, tree = taxa$treeData, weighted=TRUE, normalized=TRUE)
#' 
#' # Unweighted UniFrac
#' unifrac(x = taxa$countData, tree = taxa$treeData, weighted=FALSE)
#' @importFrom RcppParallel setThreadOptions
#' @importFrom Matrix sparseMatrix
#' @importFrom stats as.dist
#' @export

unifrac <- function(x, tree, weighted = TRUE, normalized = TRUE, threads = 1) {

    ## Error handling
    #--------------------------------------------------------------------#
    if (is.vector(x))
        cli::cli_abort("Input must a matrix of class matrix or Matrix, not a vector.")

    if (!inherits(tree, "phylo"))
        cli::cli_abort("Tree must be of class `phylo`.")

    x <- drop(as(x, "sparseMatrix"))
    if (!is.numeric(x@x))
        cli::cli_abort("Input data must be numeric.")

    if (any(x@x < 0, na.rm = TRUE))
        cli::cli_abort("Input data must be non-negative.")

    if (!is.wholenumber(threads))
        cli::cli_abort("{threads} must be a whole number.")

    ## MAIN
    #--------------------------------------------------------------------#
    RcppParallel::setThreadOptions(numThreads = threads)

    out <- .Call(
        '_OmicFlow_unifrac', PACKAGE = 'OmicFlow', 
        x, tree$edge-1, tree$edge.length, weighted, normalized
        )

    col_names <- colnames(x)
    if (!is.null(col_names))
        dimnames(out) <- list(col_names, col_names)

    return(as.dist(out))
}