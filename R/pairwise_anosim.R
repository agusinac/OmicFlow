#' Pairwise anosim (ANOSIM) computation
#'
#' @description Computes pairwise anosim from vegan, given a distance matrix and a vector of labels.
#' This function is built into the \code{ordination} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param x  A distance matrix in the form of \link[stats]{dist}.
#' Obtained from a dissimilarity metric, in the case of similarity metric please use \code{1-dist}
#' @param groups  A vector (column from a table) of labels.
#' @param p.adjust.m P adjust method see \link[stats]{p.adjust}
#' @param perm  Number of permutations to compare against the null hypothesis of anosim default: \code{perm=999}.
#' @seealso \link[vegan]{anosim}
#' @return A \code{data.frame} of
#'  * pairs that are used
#'  * R2 of H_0
#'  * p value of F^p > F
#'  * p adjusted
#'
#' @export

pairwise_anosim <- function(x, groups, p.adjust.m = "bonferroni", perm = 999){
  co <- combn(unique(as.character(groups)), 2)
  n <- ncol(co)
  pairs <- vector(mode = "numeric", length = n)
  anosimR <- vector(mode = "numeric", length = n)
  p.value <- vector(mode = "numeric", length = n)

  for(i in 1:n){
    if(inherits(x, "dist")){
      m <- as.matrix(x)[groups %in% co[, i], groups %in% co[, i]]
    }

    ano <- vegan::anosim(m, groups[groups %in% co[, i]], permutations = perm)
    pairs[i] <- paste0(co[1,i], ".vs.", co[2,i])
    anosimR[i] <- ano$statistic
    p.value[i] <- ano$signif
  }
  p.adj <- p.adjust(p.value, method = p.adjust.m)
  pairw.res <- data.frame(pairs, anosimR, p.value, p.adj)
  class(pairw.res) <- c("pwanosim", "data.frame")
  return(pairw.res)
}
