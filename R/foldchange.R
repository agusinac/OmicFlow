#' Computes non-paired Log2(A) - Log2(B) Fold Change
#'
#' @description Computes non-paired Log2(A) - Log2(B) Fold Change. The function can handle multiple categorical variables in A and B. It is possible to use multiple cores to speed up computation time.
#' This function is built into the \code{differential_feature_expression} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param dt A \code{data.table}.
#' @param sample.id A column name of a categorical variable containing sample IDs.
#' @param condition_A A vector of categorical characters, it is possible to specify multiple labels.
#' @param condition_B A vector of categorical characters, it is possible to specify multiple labels.
#' @param condition_labels A column name of a categorical variable in where \code{condtion_A} and \code{condition_B} can be found.
#' @param feature_rank A character variable of the feature level (e.g. "Genus" in taxonomy).
#' @param cpus An integer specifying the number of cores to use for embarrassing parallelism, see \link[parallel]{makeCluster}.
#' @return A list of \code{data.table} for \link{volcano_plot} and \link{ViolinBoxPlot}.
#'
#' @export


foldchange <- function(dt, sample.id, condition_A, condition_B, condition_labels, feature_rank, unique.ids, paired=FALSE) {
  # Creates tmp data table
  tmp_dt <- data.table::copy(dt)

  # subset feature labels before removing them
  feature_labels <- tmp_dt[[ feature_rank ]]
  tmp_dt <- tmp_dt[, .SD, .SDcols = !c(feature_rank)]

  # Create data.tables for results
  foldchange_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(foldchange_dt) <- feature_rank

  # Computing for multiple conditions
  for (i in seq_along(condition_A)) {
    # Subset by condition_A value
    dt_A <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[grepl(condition_A[i], condition_labels)]]
    dt_B <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[grepl(condition_B[i], condition_labels)]]

    ## In case of paired samples
    if (paired) {
      pairs <- find_pairs(dt_A, dt_B, unique.ids)
      dt_A <- dt_A[, .SD, .SDcols = pairs$colnames_A]
      dt_B <- dt_B[, .SD, .SDcols = pairs$colnames_B]
    }

    result <- base::rowMeans(dt_A) - base::rowMeans(dt_B)

    foldchange_dt <- cbind(foldchange_dt, result)
    colnames(foldchange_dt) <- c(feature_rank, paste0("Log2FC_", i))

    # Compute pvalues with wilcox test
    mat_A <- as.matrix(dt_A)
    mat_B <- as.matrix(dt_B)
    for (k in seq_along(feature_labels)) {
      # save p-values in data.table
      foldchange_dt[k, (paste0("pvalue_", i)) := stats::wilcox.test(mat_A[k, ], mat_B[k, ],
                                                                  correct = TRUE,
                                                                  paired = paired)$p.value]
    }
  }

  return(foldchange_dt)
}
