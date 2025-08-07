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
#' # Load required library
#' library(data.table)
#' 
#' # Define parameters and variables
#' sample_ids <- c("S1_A", "S2_A", "S3_A", "S4_B", "S5_B", "S6_B")
#' feature_ids <- c("Feature1", "Feature2", "Feature3")
#' 
#' # Simulated abundance matrix (features x samples)
#' abundances <- matrix(
#'   c(
#'     # Feature1 (e.g. GenusA)
#'     100, 120, 110,  55, 60, 65,
#'     # Feature2 (e.g. GenusB)
#'     50, 65, 60,    130, 120, 125,
#'     # Feature3 (e.g. GenusC)
#'     80, 85, 90,     80, 85, 90
#'   ),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(feature_ids, sample_ids)
#' )
#' 
#' # A wide table with columns as samples, rows as features
#' # And an additional column as the feature_rank, a column for feature comparison.
#' mock_data <- data.table(
#'   Genus = feature_ids,  # feature_rank column (e.g. "Genus")
#'   S1_A = abundances[ , 1],
#'   S2_A = abundances[ , 2],
#'   S3_A = abundances[ , 3],
#'   S4_B = abundances[ , 4],
#'   S5_B = abundances[ , 5],
#'   S6_B = abundances[ , 6]
#' )
#' print(mock_data)
#' 
#' # It uses substring matching, and multiple conditions can be used
#' res <- foldchange(
#'   data = mock_data,
#'   feature_rank = "Genus",
#'   condition_A = c("_A", "_B"),
#'   condition_B = c("_B", "_A"),
#' 
#'   # This can also be a column wherein, conditions A and B are present
#'   condition_labels = sample_ids, 
#'   paired = FALSE
#' )
#' print(res)
#' 
#' #---------------------#
#' ##      PAIRED       ##
#' #---------------------#
#' library(data.table)
#' 
#' # Define paired sample ids for 3 pairs:
#' paired_ids <- paste0("Pair", 1:3)
#' 
#' # Features:
#' feature_ids <- c("Feature1", "Feature2", "Feature3")
#' 
#' # Simulate abundances for each paired sample:
#' # For each pair, we have two samples: condition A and condition B.
#' # Make sure the length of condition A and condition B are the same!
#' 
#' # Construct the data.table with features as rows
#' mock_data_paired <- data.table(
#'   Genus = feature_ids,
#'   Pair1_A = c(100, 50, 80),
#'   Pair1_B = c(60, 130, 75),
#'   Pair2_A = c(120, 65, 85),
#'   Pair2_B = c(60, 120, 90),
#'   Pair3_A = c(110, 60, 90),
#'   Pair3_B = c(65, 125, 85) 
#' )
#' 
#' res <- foldchange(
#'   data = mock_data_paired,
#'   feature_rank = "Genus",
#'   condition_A = c("_A", "_B"),
#'   condition_B = c("_B", "_A"),
#' 
#'   # This can also be a column wherein, conditions A and B are present
#'   condition_labels = names(mock_data_paired)[-1], 
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
    cli::cli_abort("data must be a data.frame or data.table.")

  if (!is.character(feature_rank) && length(feature_rank) != 1) {
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
    ## TODO: convert dt_A to a sparseMatrix
    dt_A <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[grepl(condition_A[i], condition_labels)]]
    dt_B <- tmp_dt[, .SD, .SDcols = colnames(tmp_dt)[grepl(condition_B[i], condition_labels)]]

    # Convert to dense matrix
    mat_A <- as.matrix(dt_A)
    mat_B <- as.matrix(dt_B)

    # Feature means per condition
    row_means_A <- Matrix::rowMeans(mat_A)
    row_means_B <- Matrix::rowMeans(mat_B)

    # Empty vector
    result <- numeric(length(row_means_A))

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

    # Combines to final foldchange data table
    foldchange_dt <- cbind(foldchange_dt, result)
    colnames(foldchange_dt)[grepl("result", colnames(foldchange_dt))] <- paste0("Log2FC_", i)

    # Compute pvalues with wilcox test
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
