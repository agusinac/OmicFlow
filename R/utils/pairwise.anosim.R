#' Pairwise multilevel comparison using anosim
#'
#' @param x Data matrix or data frame in which rows are samples and columns are response variable(s),
#' or a dissimilarity object or a symmetric square matrix of dissimilarities
#' @param grouping Factor for grouping observations.
#' @param sim.function Function used to calculate the similarity matrix,
#' one of 'daisy' or 'vegdist' default is 'vegdist'. Ignored if x is a distance matrix.
#' @param sim.method Similarity method from daisy or vegdist, default is 'bray'. Ignored if x is a distance matrix.
#' @param p.adjust.m The p.value correction method, one of the methods supported by p.adjust(),
#' default is 'bonferroni'.
#' @param perm permutations
#'
#' @return Table with the pairwise factors, anosim R, p.value and adjusted p.value.
#' @export
#' @author ZH Gai
#' @importFrom stats p.adjust
#' @importFrom utils combn
#' @importFrom vegan anosim vegdist
#' @importFrom cluster daisy
#' @examples
#' anosim.pairwise(iris[, 1:4], iris$Species)
pairwise.anosim <- function(x, grouping,
                            sim.function = "vegdist",
                            sim.method = "bray",
                            p.adjust.m = "bonferroni",
                            perm = 999){
  co <- combn(unique(as.character(grouping)), 2) # pairwise grouping
  n <- ncol(co)
  pairs <- vector(mode = "numeric", length = n)
  anosimR <- vector(mode = "numeric", length = n)
  p.value <- vector(mode = "numeric", length = n)
  
  for(i in 1:n){
    if(inherits(x, "dist")){
      m <- as.matrix(x)[grouping %in% co[, i], grouping %in% co[, i]]
    }else if (sim.function == "daisy") {
      m <- cluster::daisy(x[grouping %in% co[, i], ], metric = sim.method)
    }else{
      m <- vegan::vegdist(x[grouping %in% co[, i], ], method = sim.method)
    }
    
    ano <- vegan::anosim(m, grouping[grouping %in% co[, i]], permutations = perm)
    pairs[i] <- paste0(co[1,i], ".vs.", co[2,i])
    anosimR[i] <- ano$statistic
    p.value[i] <- ano$signif
  }
  p.adj <- p.adjust(p.value, method = p.adjust.m)
  pairw.res <- data.frame(pairs, anosimR, p.value, p.adj)
  class(pairw.res) <- c("pwanosim", "data.frame")
  return(pairw.res)
}