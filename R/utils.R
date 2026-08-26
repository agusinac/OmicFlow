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

#' Perform pairwise or omnibus statistical tests between groups
#'
#' Performs statistical comparisons between groups using various tests from the
#' \code{matrixTests} package. Supports location tests (t-test, Wilcoxon, Kruskal-Wallis,
#' Van der Waerden) and scale tests (Levene, Bartlett, Brown-Forsythe, Fligner-Killeen).
#' 
#' @param data A \link[data.table]{data.table} or \link[base]{data.frame}.
#' @param x_col A character string specifying the column with numeric values in \code{data}
#' @param g_col A character string specifying the column with group labels in \code{data}
#' @param test A character string specifying the statistical test:
#'   \describe{
#'     \item{\strong{Pairwise tests:}}{
#'       \itemize{
#'         \item \code{"t.test"} - Welch t-test (default, unequal variances assumed)
#'         \item \code{"t.equalvar"} - Equal variance t-test (Student's t-test)
#'         \item \code{"wilcox"} - Wilcoxon rank-sum/Mann-Whitney U test
#'       }
#'     }
#'   }
#' @param paired A boolean value whether to use paired tests (default: \code{FALSE}).
#'   For \code{test = "t.test"}, uses \link[matrixTests]{row_t_paired} when \code{TRUE}
#'   or \link[matrixTests]{row_t_welch} when \code{FALSE}.
#'   For \code{test = "wilcox"}, uses \link[matrixTests]{row_wilcoxon_paired} when \code{TRUE}
#'   or \link[matrixTests]{row_wilcoxon_twosample} when \code{FALSE}.
#' @param p.adjust.method A character string to specify the p-adjust method to use in `stats::p.adjust` (default: \code{"fdr"}).
#' @param ... Extra arguments passed to the underlying \code{matrixTests} function.
#' @return A data.table with pairwise comparison results including: \describe{
#'   \item{group1, group2}{The two groups being compared}
#'   \item{obs.x, obs.y, obs.tot}{Number of observations in each group and total}
#'   \item{statistic}{Test statistic}
#'   \item{pvalue}{Raw p-value}
#'   \item{p.adj}{P-value adjusted for multiple testing}
#'   \item{y.position}{Y position for plotting significance brackets}
#'   \item{xmin, xmax}{X positions for plotting significance brackets}
#'   \item{Additional columns}{Test-specific statistics (e.g., conf.low, conf.high for t-tests)}
#' }
#' 
#' For omnibus tests (≥2 groups), returns a single-row result with overall test statistics.
#' @examples
#' \dontrun{
#' # Pairwise Welch t-test (2 groups)
#' pairwise_test(data, x_col = "value", g_col = "group", test = "t.test")
#'
#' # Pairwise Wilcoxon (2 groups)
#' pairwise_test(data, x_col = "value", g_col = "group", test = "wilcox")
#' }
#' @noRd
pairwise_test <- function(
  data, 
  x_col, 
  g_col, 
  test = "wilcox", 
  paired = FALSE, 
  p.adjust.method = "fdr", ...) {
  
  ## Error handling
  #--------------------------------------------------------------------#
  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("{.val data} must be a {.cls data.frame} or {.cls data.table}.")
  
  if (!is.character(x_col) || length(x_col) != 1) {
    cli::cli_abort("{.val x_col} needs to contain characters with length of 1.")
  } else if (!column_exists(x_col, data)) {
    cli::cli_abort("The {.val x_col} column does not exist in the provided {.arg data}.")
  }
  if (!is.character(g_col) || length(g_col) != 1) {
    cli::cli_abort("{.val g_col} needs to contain characters with length of 1.")
  } else if (!column_exists(g_col, data)) {
    cli::cli_abort("The {.val g_col} column does not exist in the provided {.arg data}.")
  }

  # Define valid methods by category
  OPTIONS <- c("t.test", "t.equalvar", "wilcox")
  if (!test %in% OPTIONS) {
    cli::cli_abort(
      "{.val {test}} is not a valid test.\nValid options: {.val {OPTIONS}}."
    )
  }

  if (!is.logical(paired))
    cli::cli_abort("{.val paired} needs to be either `TRUE` or `FALSE`.")

  if (!is.character(p.adjust.method)) {
    cli::cli_abort("{.val p.adjust.method} must be a character.")
  } else if (!c(p.adjust.method %in% stats::p.adjust.methods)) {
    cli::cli_abort("{.val {p.adjust.method}} is not a valid method. \nValid options: {.val {p.adjust.methods}}.")
  }

  ## MAIN
  #--------------------------------------------------------------------#
  data_tmp <- data.table::copy(data)
  pvalue <- group1 <- group2 <- NULL

  # Initialize required parameters
  unique_groups <- unique(as.character(data_tmp[[ g_col ]]))
  n_groups <- length(unique_groups)

  if (n_groups < 2) {
    cli::cli_abort("At least 2 groups are required for test {.val {test}}.")
  }

  # Fetch right test
  get_test <- function(test, paired) {
    if (paired) {
      switch(test,
        "t.test" = matrixTests::row_t_paired,
        "t.equalvar" = matrixTests::row_t_paired,
        "wilcox" = matrixTests::row_wilcoxon_paired,
        cli::cli_abort("Invalid paired test.")
      )
    } else {
      switch(test,
        "t.test" = matrixTests::row_t_welch,
        "t.equalvar" = matrixTests::row_t_equalvar,
        "wilcox" = matrixTests::row_wilcoxon_twosample,
        cli::cli_abort("Invalid test.")
      )
    }
  }

  test_fun <- get_test(test, paired)

  # Initialize pair combinations
  co <- utils::combn(unique_groups, 2)
  n <- ncol(co)
  out_list <- list()

  ## Grouped max values for Y positioning
  groups_max <- data_tmp[, max(.SD[[x_col]], na.rm = TRUE), by = g_col]

  # Loop through pairs
  for (i in seq_len(n)) {
    pair_1 <- co[1, i]
    pair_2 <- co[2, i]

    X <- data_tmp[[x_col]][data_tmp[[g_col]] %in% pair_1]
    Y <- data_tmp[[x_col]][data_tmp[[g_col]] %in% pair_2]

    # Run the test
    out <- test_fun(x = X, y = Y)

    # Saving stats
    out[["group1"]] <- as.character(pair_1)
    out[["group2"]] <- as.character(pair_2)
    out[["y.position"]] <- max(groups_max[groups_max[[g_col]] %in% c(pair_1, pair_2), ]$V1) * 1.01

    out_list[[i]] <- out
  }
  # Combine pairwise subsets, adjust p-value, set new order
  pairw.res <- data.table::rbindlist(out_list)
  pairw.res[, "p.adj" := data.table::fifelse(pvalue == 1, pvalue,
                                             stats::p.adjust(pvalue, method = p.adjust.method))]

  # Define column order based on method
  base_cols <- c("group1", "group2", "obs.x", "obs.y", "obs.tot", "statistic",
                 "pvalue", "p.adj", "alternative", "y.position")

  # Add method-specific columns
  if (test %in% c("t.test", "t.equalvar")) {
    extra_cols <- c("conf.low", "conf.high", "mean.x", "mean.y", "mean.diff")
  } else if (test == "wilcox") {
    extra_cols <- c("location.null", "exact", "corrected")
  }

  if (paired) {
    base_cols[5] <- "obs.paired"
  }

  col_order <- c(base_cols, extra_cols)
  col_order <- intersect(col_order, names(pairw.res))

  data.table::setcolorder(x = pairw.res, neworder = col_order)

  ## Adding X positions
  pairw.res[, "xmin" := as.numeric(as.factor(group1))]
  pairw.res[, "xmax" := as.numeric(as.factor(group2)) + 1]

  return(pairw.res)
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