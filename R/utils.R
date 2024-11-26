sparse_to_dtable <- function(sparsemat) {
  return(data.table::data.table(as.matrix(sparsemat)))
}
