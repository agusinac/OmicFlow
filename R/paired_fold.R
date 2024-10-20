#' Computes paired Log2(A) - Log2(B) Fold Change given a paired.id
#'
#' @description Computes paired Log2(A) - Log2(B) Fold Change of paired samples. The function can handle multiple categorical variables in A and B. It is possible to use multiple cores to speed up computation time.
#' This function is built into the \code{differential_feature_expression} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param dt A \code{data.table}.
#' @param sample.id A column name of a categorical variable containing sample IDs.
#' @param paired.id A column name of a categorical variable containing paired sample IDs.
#' @param unique.id A list of unique character variables from the paired.id column.
#' @param condition_A A vector of categorical characters, it is possible to specify multiple labels.
#' @param condition_B A vector of categorical characters, it is possible to specify multiple labels.
#' @param condition_labels A column name of a categorical variable in where \code{condtion_A} and \code{condition_B} can be found.
#' @param feature_rank A character variable of the feature level (e.g. "Genus" in taxonomy).
#' @param cpus An integer specifying the number of cores to use for embarrassing parallelism, see \link[parallel]{makeCluster}.
#' @return A list of \code{data.table} for \link{volcano_plot} and \link{ViolinBoxPlot}.
#'
#' @export

paired_fold <- function(dt, sample.id, paired.id, unique.id, condition_A, condition_B, feature_rank, condition_labels, cpus = 8) {
  # tmp data.table
  tmp_dt <- data.table::copy(dt)

  feature_labels <- tmp_dt[[ feature_rank ]]
  paired_dt <- data.table::data.table(feature_rank = feature_labels)
  colnames(paired_dt) <- feature_rank

  # Register Parallel backend
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)

  results <- foreach(id = seq_along(unique.id), .combine = rbind, .packages = 'data.table') %dopar% {
    # extract samples, takes only matching pairs
    pair <- colnames(dt)[grepl(unique.id[id], colnames(dt))]
    if (length(pair) == 2) {
      sample_pair <- tmp_dt[, .SD, .SDcols = pair]
    } else {
      sample_pair <- NULL
    }

    # Computes 2foldChange for matching pairs
    if (!is.null(sample_pair)) {
      # initialize empty table
      dt_diff <- data.table::data.table()

      for (i in seq_along(condition_A)) {
        # Subset column pair by condition A and B
        dt_A <- sample_pair[, .SD, .SDcols = colnames(sample_pair)[grepl(condition_A[i], condition_labels)][1]]
        dt_B <- sample_pair[, .SD, .SDcols = colnames(sample_pair)[grepl(condition_B[i], condition_labels)][1]]

        # subtraction
        dt_diff[, (paste0("diff_", i)) := dt_A - dt_B]
        dt_diff[, (feature_rank) := feature_labels]
      }
      dt_diff[, (sample.id) := unique.id[id]]
    }
  }
  # Wraps and ends Parallel backend
  parallel::stopCluster(cl)
  paired_dt <- na.omit(rbind(paired_dt, results, fill = TRUE))

  return(paired_dt)
}
