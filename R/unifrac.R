#' Compute UniFrac Dissimilarity from a Dense or Sparse Matrix.
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
#' @inheritParams bray
#' @param tree A `phylo` class tree.
#' @param normalize A boolean value, whether to normalize weighted UniFrac distances to be between 0 and 1 (default: \code{TRUE}). Unweighted UniFrac is always normalized.
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
#' taxa$scale(method = "tss")
#'
#' # Weighted UniFrac
#' unifrac(x = taxa$countData, tree = taxa$treeData, weighted=TRUE, normalize=FALSE)
#' 
#' # Weighted Normalized UniFrac
#' unifrac(x = taxa$countData, tree = taxa$treeData, weighted=TRUE, normalize=TRUE)
#' 
#' # Unweighted UniFrac
#' unifrac(x = taxa$countData, tree = taxa$treeData, weighted=FALSE)
#' @importFrom Matrix sparseMatrix
#' @export

unifrac <- function(x, tree, weighted = TRUE, normalize = TRUE, threads = 1) {

    ## Error handling
    #--------------------------------------------------------------------#
    if (!inherits(tree, "phylo"))
        cli::cli_abort("{.val tree} must be a {.cls phylo}.")

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

    if (!is.logical(normalize))
        cli::cli_abort("{.val normalize} needs to be either `TRUE` or `FALSE`.")

    if (length(threads) != 1) {
        cli::cli_abort("{.val threads} must be a single whole number.")
    } else if (!is.wholenumber(threads)) {
        cli::cli_abort("{.val threads} must be a whole number.")
    }

    ## MAIN
    #--------------------------------------------------------------------#
    RcppParallel::setThreadOptions(numThreads = threads)

    ## Check if rows are aligned by tip.labels

    out <- .Call(
        '_OmicFlow_unifrac', PACKAGE = 'OmicFlow', 
        x, tree$edge-1, tree$edge.length, weighted, normalize
        )

    col_names <- colnames(x)
    if (!is.null(col_names))
        dimnames(out) <- list(col_names, col_names)

    return(stats::as.dist(out))
}