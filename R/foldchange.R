#' Computes Log2(A) - Log2(B) Fold Change of (non-) paired data.
#'
#' @description Computes (non-)paired Log2(A) - Log2(B) Fold Change.
#' This function is built into the class \link{omics} with method \code{DFE()} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}. The function handles zero's, and doesn't return +/- infinites.
#'
#' @param data A \link[data.table]{data.table}.
#' @param feature_rank A character variable of the feature level (e.g. "Genus" in taxonomy).
#' @param condition_A A vector of categorical characters, it is possible to specify multiple labels.
#' @param condition_B A vector of categorical characters, it is possible to specify multiple labels.
#' @param condition_labels A vector character wherein \code{condition_A} and \code{condition_B} are present.
#' @param paired A Boolean value to perform paired or non-paired test, see \link[stats]{wilcox.test}.
#' @return A \link[data.table]{data.table}
#' @examples 
#' #-------------------------#
#' ##      NON-PAIRED       ##
#' #-------------------------#
#' library(data.table)
#'
#' # Define parameters and variables
#' sample_ids <- c("S1", "S2", "S3", "S4", "S5", "S6")
#' groups <- c("A", "A", "B", "B", "C", "C")
#' feature_ids <- c("Feature1", "Feature2", "Feature3")
#'
#' # Simulated abundance matrix (features x samples)
#' abundances <- matrix(
#'   c(
#'     100, 120, 110, 55, 60, 65,
#'     50, 65, 60, 130, 120, 125,
#'     80, 85, 90, 80, 85, 90
#'   ),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(feature_ids, sample_ids)
#' )
#'
#' # Convert to a data.table
#' mock_data <- OmicFlow::matrix_to_dtable(abundances)
#' mock_data$Genus <- feature_ids
#'
#' # It uses exact matching and multiple conditions are allowed.
#' res <- foldchange(
#'   data = mock_data,
#'   feature_rank = "Genus",
#'   condition_A = c("A", "B"),
#'   condition_B = c("B", "C"),
#'   condition_labels = groups,
#'   paired = FALSE
#' )
#' print(res)
#' 
#' #---------------------#
#' ##      PAIRED       ##
#' #---------------------#
#' library(data.table)
#' 
#' # In the paired case both conditions A and B must be of the same length!
#' # We re-use the above mock_data and only change group labels
#' groups <- c("A", "A", "B", "B", "A", "B")
#' 
#' res <- foldchange(
#'   data = mock_data,
#'   feature_rank = "Genus",
#'   condition_A = c("A"),
#'   condition_B = c("B"),
#'   condition_labels = groups,
#'   paired = TRUE
#' )
#' print(res)
#' 
#' @export

foldchange <- function(data,
                       feature_rank,
                       condition_A,
                       condition_B,
                       condition_labels,
                       paired = FALSE) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("Data must be a {.cls data.frame} or {.cls data.table}.")

  if (!is.character(feature_rank) && length(feature_rank) != 1) {
    cli::cli_abort("{.val {feature_rank}} needs to contain characters with length of 1.")
  } else if (!column_exists(feature_rank, data)) {
    cli::cli_abort("The {.val {feature_rank}} column does not exist in the provided {.arg data}.")
  }

  if (!is.vector(condition_labels))
    cli::cli_abort("{.val {condition_labels}} needs to be {.cls vector}.")

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
    ## TODO: convert dt_A to a sparseMatrix
    dt_A <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[condition_labels %in% condition_A[i]]]
    dt_B <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[condition_labels %in% condition_B[i]]]

    # Convert to dense matrix
    mat_A <- as.matrix(dt_A)
    mat_B <- as.matrix(dt_B)

    # Feature means per condition
    row_means_A <- Matrix::rowMeans(mat_A)
    row_means_B <- Matrix::rowMeans(mat_B)

    # Empty vector
    result <- numeric(length(row_means_A))
    max_val <- base::max(mat_A)

    # Find zero's to prevent Inf
    both_zero <- row_means_A == 0 & row_means_B == 0
    row_means_A_zero <- row_means_A == 0 & row_means_B != 0
    row_means_B_zero <- row_means_A != 0 & row_means_B == 0
    both_non_zero <- row_means_A != 0 & row_means_B != 0

    # Compute log2 fold change
    result[both_zero] <- 0
    result[row_means_A_zero] <- row_means_A[row_means_A_zero] - log2(row_means_B[row_means_A_zero])
    result[row_means_B_zero] <- log2(row_means_A[row_means_B_zero]) - row_means_B[row_means_B_zero]
    result[both_non_zero] <- log2(row_means_A[both_non_zero]) - log2(row_means_B[both_non_zero])

    # Reverse flipped values with zero's based on max_val
    if (max_val < 1.0) {
      result[row_means_A_zero] <- result[row_means_A_zero] * -1
      result[row_means_B_zero] <- result[row_means_B_zero] * -1
    }

    # Combines to final foldchange data table
    foldchange_dt <- cbind(foldchange_dt, result)
    colnames(foldchange_dt)[grepl("result", colnames(foldchange_dt))] <- paste0("Log2FC_", i)

    # Compute pvalues with wilcox test
    for (k in seq_along(feature_labels)) {
      # save p-values in data.table
      suppressWarnings(
        foldchange_dt[
          k, (paste0("pvalue_", i)) := stats::wilcox.test(
            mat_A[k, ], mat_B[k, ],
            correct = TRUE,
            paired = paired
            )$p.value
          ]
      )
    }
  }

  return(foldchange_dt)
}
