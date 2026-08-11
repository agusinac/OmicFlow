#' @noRd
#' When a function is called via `pkg::fun` within R6 class, 
#' then devtools will not see this and deliver a note: `Namespace in Imports field not declared from: matrixTests`
#' This is a placeholder to prevent devtools note appearing.
.check_package_dependencies <- function() {
  rhdf5::h5read
  data.table::fread
  Matrix::Matrix
  ggplot2::ggplot
  matrixTests::row_levene
  matrixStats::rowMeans2
  R6::R6Class
  ape::read.tree
  yyjsonr::read_json_file
  jsonvalidate::json_validate
  patchwork::plot_layout
  invisible(NULL)
}