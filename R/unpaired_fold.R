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


unpaired_fold <- function(dt, sample.id, condition_A, condition_B, condition_labels, feature_rank, cpus = 8) {
  # Creates tmp data table
  tmp_dt <- data.table::copy(dt)

  # subset feature labels before removing them
  feature_labels <- tmp_dt[[ feature_rank ]]
  tmp_dt <- tmp_dt[, .SD, .SDcols = !c(feature_rank)]

  # Create data.tables for results
  unpaired_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(unpaired_dt) <- feature_rank

  # Register Parallel backend
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)

  # Computing for multiple conditions
  for (i in seq_along(condition_A)) {
    # Subset by condition_A value
    dt_A <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[grepl(condition_A[i], condition_labels)]]
    dt_B <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[grepl(condition_B[i], condition_labels)]]

    # Improve with foreach and parallelize it!
    # Create cross-wise combinations
    # combinations <- data.table::data.table(expand.grid(colnames(dt_A), colnames(dt_B)))

    # Perform subtraction cross wise
    # results <- foreach(j = seq_along(combinations$Var1), .combine = cbind, .packages = 'data.table') %dopar% {
    #   col_A <- combinations$Var1[j]
    #   col_B <- combinations$Var2[j]
    #
    #   dt_A[, ..col_A] - dt_B[, ..col_B]
    # }
    result <- base::rowMeans(dt_A) - base::rowMeans(dt_B)

    unpaired_dt <- cbind(unpaired_dt, result)
    colnames(unpaired_dt) <- c(feature_rank, paste0("Log2FC_", i))

    # Compute pvalues with wilcox test
    mat_A <- as.matrix(dt_A)
    mat_B <- as.matrix(dt_B)
    for (k in seq_along(feature_labels)) {
      # save p-values in data.table
      unpaired_dt[k, (paste0("pvalue_", i)) := stats::wilcox.test(mat_A[k, ], mat_B[k, ], correct = TRUE)$p.value]
    }
  }
  # Stop the cluster
  parallel::stopCluster(cl)

  return(unpaired_dt)
}
