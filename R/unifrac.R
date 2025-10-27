#' Compute UniFrac Dissimilarity from a Sparse Matrix.
#'
#' @description Calculates the UniFrac dissimilarity between samples based on phylogenetic branch lengths and abundance or presence/absence data.
#' The UniFrac distance between two samples \eqn{A} and \eqn{B}, with phylogenetic tree edges \eqn{i = 1 \ldots n} of lengths \eqn{L_i}, is computed differently depending on the \code{weighted} and \code{normalized} flags:
#' Weighted UniFrac:
#' \eqn{d(A,B) = \frac{\sum_{i}^n L_i |A_i - B_i|}{\sum_{i}^n L_i (A_i + B_i)}}
#' Normalized Weighted UniFrac:
#' \eqn{d(A,B) = \sum_{i}^n L_i |A_i - B_i|}
#' where \eqn{A_i} and \eqn{B_i} are the abundance weights for branch \eqn{i} computed by propagating tip abundances to internal nodes.
#'
#' Unweighted UniFrac (always normalized): 
#' \eqn{d(A,B) = \frac{\sum_{i}^n L_i |A'_i - B'_i|}{\sum_{i}^n L_i \max(A'_i, B'_i)}}
#'
#' where \eqn{A'_i} and \eqn{B'_i} are binary indicators (presence/absence) on branches \eqn{i}, derived by propagating tip presence to internal nodes.
#' When \code{weighted = FALSE}, input counts are first converted to presence/absence before distance computation.
#'
#' @param x A \link[Matrix]{sparseMatrix} of strictly positive counts or presence/absence data.
#' @param tree A `phylo` class tree.
#' @param weighted Logical indicating whether to compute weighted (abundance) or unweighted (presence/absence) UniFrac (default: TRUE).
#' @param normalized Logical indicating whether to normalize weighted UniFrac distances to be between 0 and 1 (default: TRUE). Unweighted UniFrac is always normalized.
#' @param threads Integer number of threads to use for parallel computation (default: 1).
#' @return A symmetric matrix of pairwise UniFrac distances between columns of \code{x}.
#' @references
#' Lozupone, C., Hamady, M., Kelley, S. T., & Knight, R. (2007). Quantitative and qualitative beta diversity measures lead to different insights into factors that structure microbial communities. Applied and Environmental Microbiology, 73(5), 1576–1585.
#'
#' Lozupone, C., & Knight, R. (2005). UniFrac: a new phylogenetic method for comparing microbial communities. Applied and Environmental Microbiology, 71(12), 8228–8235.
#'
#' Lozupone, C., & Knight, R. (2008). Species divergence and the measurement of microbial diversity. FEMS Microbiology Reviews, 32(4), 557–578.
#'
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