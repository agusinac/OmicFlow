#' Converting a Matrix to data.table
#'
#' @description Wrapper function that converts a sparseMatrix to data.table
#'
#' @param x A \link[base]{matrix}, \link[Matrix]{sparseMatrix} or \link[Matrix]{Matrix}.
#' @return A \link[data.table]{data.table} class.
#' @export
matrix_to_dtable <- function(x) {
  if (inherits(x, "denseMatrix") || inherits(x, "matrix") || inherits(x, "sparseMatrix")) {
      return(data.table::data.table(as.matrix(x)))
  } else cli::cli_abort("{.val x} isn't a {.cls matrix}, {.cls denseMatrix} or {.cls sparseMatrix}.")
}
#' Checks if column exists in table
#'
#' @description Mainly used within \link{omics} and other functions to check if given column name(s) exist in the table and is not completely empty (containing NAs).
#'
#' @param column A character vector.
#' @param table A \link[data.table]{data.table} or \link[base]{data.frame}.
#' @return A boolean value.
#' @export
column_exists <- function(column, table) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!is.character(column))
    cli::cli_abort("{.val column} needs to contain characters.")

  if (!inherits(table, "data.frame") && !inherits(table, "data.table"))
    cli::cli_abort("{.val table} must be a {.cls data.frame} or {.cls data.table}.")

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

#' Check if number is integer
#'
#' @description Checks if a given value `x` is a wholenumber, so it should be 1, 4 and not a float.
#' \link{is.integer} also accepts floats.
#' 
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.character(x)) {
    return(FALSE)
  } else {
    abs(x - round(x)) < tol
  }
}

#' Check if input is a valid color
#'
#' @description Checks if a given value `x` contains red, green and blue channels.
#' 
#' @noRd
is.color <- function(x) {
  return(
    all(sapply(x, function(X) {
      tryCatch(is.matrix(grDevices::col2rgb(X)), 
              error = function(e) FALSE)
      }))
  )
}

#' Helper function in `omics$autoFlow` to combine conditions 
#' 
#' @noRd
combine_conditions <- function(condition1, condition2) {
  if (!inherits(condition1, "data.frame") && !inherits(condition1, "data.table"))
    cli::cli_abort("{.val condition1} must be a {.cls data.frame} or {.cls data.table}.")

  if (!inherits(condition2, "data.frame") && !inherits(condition2, "data.table"))
    cli::cli_abort("{.val condition2} must be a {.cls data.frame} or {.cls data.table}.")

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