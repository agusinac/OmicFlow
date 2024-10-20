#' Pairwise adonis2 (PERMANOVA) computation
#'
#' @description Computes pairwise adonis2 from vegan, given a distance matrix and a vector of labels.
#' This function is built into the \code{ordination} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param x  A distance matrix in the form of \link[stats]{dist}.
#' Obtained from a dissimilarity metric, in the case of similarity metric please use \code{1-dist}
#' @param groups  A vector (column from a table) of labels.
#' @param p.adjust.m P adjust method see \link[stats]{p.adjust}
#' @param perm  Number of permutations to compare against the null hypothesis of adonis2 default: \code{perm=999}.
#' @seealso \link[vegan]{adonis2}
#' @return A \code{data.frame} of
#'  * pairs that are used
#'  * Degrees of freedom (Df)
#'  * Sums of Squares of H_0
#'  * F.Model of H_0
#'  * R2 of H_0
#'  * p value of F^p > F
#'  * p adjusted
#'
#' @export

pairwise_adonis <- function(x, groups, p.adjust.m = "bonferroni", perm = 999){
  # Initialize required parameters
  co <- combn(unique(as.character(groups)), 2)
  n <- ncol(co)
  pairs <- vector(mode = "numeric", length = n)
  p.value <- vector(mode = "numeric", length = n)
  Df <- vector(mode = "numeric", length = n)
  SumsOfSqs <- vector(mode = "numeric", length = n)
  F.Model <- vector(mode = "numeric", length = n)
  R2 <- vector(mode = "numeric", length = n)

  # Loops through pairs
  for(i in 1:n){
    if(inherits(x, 'dist')){
      m = as.matrix(x)[groups %in% c(as.character(co[1, i]),as.character(co[2, i])),
                      groups %in% c(as.character(co[1, i]),as.character(co[2, i]))]
    }
    # Performing adonis2 test
    tmp_m = data.frame(Fac = groups[groups %in% c(co[1,i], co[2,i])])
    ad <- vegan::adonis2(m ~ Fac,
                         data = tmp_m,
                         permutations = perm);

    # Saving stats
    pairs[i] <- paste(co[1, i],'vs',co[2, i])
    Df[i] <- ad$Df[1]
    SumsOfSqs[i] <- ad$SumOfSqs[1]
    F.Model[i] <- ad$F[1]
    R2[i] <- ad$R2[1]
    p.value[i] <- ad$`Pr(>F)`[1]
  }
  # Adjusts P-values and returns combined dataframe
  p.adj <- p.adjust(p.value, method = p.adjust.m)
  pairw.res <- data.frame(pairs, Df, SumsOfSqs, F.Model, R2, p.value, p.adj)
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
}
