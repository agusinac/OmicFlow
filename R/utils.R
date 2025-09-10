#' Converting a sparse matrix to data.table
#'
#' @description Wrapper function that converts a sparseMatrix to data.table
#'
#' @param sparsemat A \link[Matrix]{sparseMatrix} class.
#' @return A \link[data.table]{data.table} class.
#' @export
sparse_to_dtable <- function(sparsemat) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(sparsemat, "sparseMatrix"))
    cli::cli_abort("sparsemat must be a sparseMatrix.")

  ## MAIN
  #--------------------------------------------------------------------#

  return(data.table::data.table(as.matrix(sparsemat)))
}

#' Loads a rarefied alpha diversity table from Qiime2
#'
#' @description Parses a QIIME2 table of rarefied data into a data.table as input to \link{diversity_plot}
#'
#' @param filepath A character value, filename or filepath to existing file.
#' @return A \link[data.table]{data.table}.
#' @export
read_rarefraction_qiime <- function(filepath) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!file.exists(filepath))
    cli::cli_abort("{filepath} does not exist.")

  ## MAIN
  #--------------------------------------------------------------------#

  df_shannon <- data.table::fread(filepath)

  # Pivot into long table
  shannon_long <- data.table::melt(data = df_shannon,
                                   measure.vars = colnames(df_shannon)[grepl("depth-", colnames(df_shannon))],
                                   variable.name = "iters",
                                   variable.factor = FALSE,
                                   value.name = "alpha_div")
  # Corrects colnames
  colnames(shannon_long) <- c("SAMPLE_ID", "iters", "alpha_div")

  return(shannon_long)
}

#' Checks if column exists in table
#'
#' @description Mainly used within \link{omics} and other functions to check if given column name does exist in the table and is not completely empty (containing NAs).
#'
#' @param column A character of length 1.
#' @param table A \link[data.table]{data.table} or \link[base]{data.frame}.
#' @return A boolean value.
#' @export
column_exists <- function(column, table) {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!is.character(column) && length(column) != 1)
    cli::cli_abort("{column} needs to contain characters with length of 1.")

  if (!inherits(table, "data.frame") && !inherits(table, "data.table"))
    cli::cli_abort("table must be a data.frame or data.table.")

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

combine_conditions <- function(condition1, condition2) {
  if (!is.null(condition1) && !is.null(condition2)) {
    if (!inherits(condition1, "data.frame") && !inherits(condition1, "data.table"))
      cli::cli_abort("condition1 must be a data.frame or data.table.")

    if (!inherits(condition2, "data.frame") && !inherits(condition2, "data.table"))
      cli::cli_abort("condition2 must be a data.frame or data.table.")
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

  # Find which condition2 are NOT already in condition1
  new_pairs_idx <- !cond2_str %in% cond1_str

  if (any(new_pairs_idx)) {
    # There are new pairs in condition2 not in condition1;
    # append only the new ones
    new_rows <- condition2[new_pairs_idx, ]
    updated_conditions <- rbind(condition1, new_rows)
  } else {
    # All pairs in condition2 are already included
    updated_conditions <- condition1
  }

  return(updated_conditions)
}