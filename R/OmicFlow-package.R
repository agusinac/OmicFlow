#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table := .SD
#' @importFrom Rcpp sourceCpp
#' @useDynLib OmicFlow, .registration = TRUE
#' @importFrom utils globalVariables
## usethis namespace: end

utils::globalVariables(c(".", "group_col"))
NULL
