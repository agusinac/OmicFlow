#' Converting a Matrix to data.table
#'
#' @description Wrapper function that converts a sparseMatrix to data.table
#'
#' @param x A \link[base]{matrix}, \link[Matrix]{sparseMatrix} or \link[Matrix]{Matrix}.
#' @return A \link[data.table]{data.table} class.
#' @noRd
matrix_to_dtable <- function(x) {
  if (inherits(x, "denseMatrix") || inherits(x, "matrix") || inherits(x, "sparseMatrix")) {
      return(data.table::data.table(as.matrix(x)))
  } else cli::cli_abort("Input isn't a {.cls matrix}, {.cls denseMatrix} or {.cls sparseMatrix}.")
}
#' Checks if column exists in table
#'
#' @description Mainly used within \link{omics} and other functions to check if given column name does exist in the table and is not completely empty (containing NAs).
#'
#' @param column A character of length 1.
#' @param table A \link[data.table]{data.table} or \link[base]{data.frame}.
#' @return A boolean value.
#' @noRd
column_exists <- function(column, table) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!is.character(column) && length(column) != 1)
    cli::cli_abort("{.val {column}} needs to contain characters with length of 1.")

  if (!inherits(table, "data.frame") && !inherits(table, "data.table"))
    cli::cli_abort("{.arg table} must be a {.cls data.frame} or {.cls data.table}.")

  ## MAIN
  #--------------------------------------------------------------------#

  valid_columns <- column[column %in% colnames(table)]

  if (length(valid_columns) == 0) {
    return(FALSE)
  }

  # For each existing column, check if it's *not entirely NA*
  columns_empty <- all(sapply(valid_columns, function(col) {
    any(!is.na(table[[col]]))
  }))

  return (length(valid_columns) == length(column) && columns_empty)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.character(x)) {
    return(FALSE)
  } else {
    abs(x - round(x)) < tol
  }
}

#' Helper function in `omics$autoFlow` to combine conditions 
#' 
#' @noRd
combine_conditions <- function(condition1, condition2) {
  if (!is.null(condition1) && !is.null(condition2)) {
    if (!inherits(condition1, "data.frame") && !inherits(condition1, "data.table"))
      cli::cli_abort("{.arg condition1} must be a {.cls data.frame} or {.cls data.table}.")

    if (!inherits(condition2, "data.frame") && !inherits(condition2, "data.table"))
      cli::cli_abort("{.arg condition2} must be a {.cls data.frame} or {.cls data.table}.")
  }

  # Combine to strings for easy comparison
  cond1_str <- paste(
    pmin(condition1$group1, condition1$group2),
    pmax(condition1$group1, condition1$group2), 
    sep = "_")

  cond2_str <- paste(
    pmin(condition2$group1, condition2$group2),
    pmax(condition2$group1, condition2$group2), 
    sep = "_")

  # Extend new unique pairs
  new_pairs_idx <- !cond2_str %in% cond1_str

  if (any(new_pairs_idx)) {
    new_rows <- condition2[new_pairs_idx, ]
    updated_conditions <- rbind(condition1, new_rows)
  } else {
    updated_conditions <- condition1
  }

  return(updated_conditions)
}