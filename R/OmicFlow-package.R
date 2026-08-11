#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table := .SD
#' @importFrom ggplot2 .data
#' @importFrom Matrix sparseMatrix
#' @useDynLib OmicFlow, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom RcppParallel setThreadOptions
#' @importFrom utils globalVariables
## usethis namespace: end

utils::globalVariables(c(".", "group_col"))
NULL
