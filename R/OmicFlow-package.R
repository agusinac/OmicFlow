#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom magrittr %>%
#' @importFrom data.table := .SD
#' @importFrom stats na.omit p.adjust
#' @useDynLib OmicFlow, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils globalVariables
## usethis namespace: end

utils::globalVariables(c(".", "group_col"))
NULL
