#' Computes Log2(A) - Log2(B) Fold Change of (non-) paired data.
#'
#' @description Computes (non-)paired Log2(A) - Log2(B) Fold Change.
#' This function is built into the \code{differential_feature_expression} method from the abstract class \link[OmicFlow]{omics} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metagenomics} and \link[OmicFlow]{proteomics}.
#'
#' @param data A \link[data.table]{data.table}.
#' @param sample.id A column name of a categorical variable containing sample IDs.
#' @param feature_rank A character variable of the feature level (e.g. "Genus" in taxonomy).
#' @param condition_A A vector of categorical characters, it is possible to specify multiple labels.
#' @param condition_B A vector of categorical characters, it is possible to specify multiple labels.
#' @param condition_labels A vector characters wherein `condition_A` and `condition_B` are present.
#' @return A \link[data.table]{data.table}
#'
#' @export

foldchange <- function(data,
                       sample.id,
                       feature_rank,
                       condition_A,
                       condition_B,
                       condition_labels,
                       paired = FALSE) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") || !inherits(data, "data.table"))
    cli::cli_abort("data must be a data.frame or data.table.")

  if (!is.character(sample.id) && length(values) != 1) {
    cli::cli_abort("Column name: {sample.id} needs to contain characters with length of 1.")
  } else if (!column_exists(sample.id, data)) {
    cli::cli_abort("The {sample.id} column does not exist in the provided data.")
  }

  if (!is.character(feature_rank) && length(values) != 1) {
    cli::cli_abort("Column name: {feature_rank} needs to contain characters with length of 1.")
  } else if (!column_exists(feature_rank, data)) {
    cli::cli_abort("The {feature_rank} column does not exist in the provided data.")
  }

  if (!is.vector(condition_labels))
    cli::cli_abort("{condition_labels} needs to be a vector.")

  ## MAIN
  #--------------------------------------------------------------------#

  # Creates tmp data table
  tmp_dt <- data.table::copy(data)

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

    result <- base::rowMeans(dt_A) - base::rowMeans(dt_B)

    foldchange_dt <- cbind(foldchange_dt, result)
    colnames(foldchange_dt) <- c(feature_rank, paste0("Log2FC_", i))

    # Compute pvalues with wilcox test
    mat_A <- as.matrix(dt_A)
    mat_B <- as.matrix(dt_B)
    for (k in seq_along(feature_labels)) {
      # save p-values in data.table
      foldchange_dt[
        k, (paste0("pvalue_", i)) := stats::wilcox.test(
          mat_A[k, ], mat_B[k, ],
          correct = TRUE,
          paired = paired
          )$p.value
        ]
    }
  }

  return(foldchange_dt)
}
