sparse_to_dtable <- function(sparsemat) {
  return(data.table::data.table(as.matrix(sparsemat)))
}

read_tsv_matrix <- function(filename) {
  distmat <- data.table::fread(filename,
                               header = TRUE)

  rownames <- distmat[[1]]
  colnames <- names(distmat)[-1]
  values <- as.matrix(distmat[, -1, with = FALSE])

  sparse_matrix <- sparseMatrix(
    i = rep(1:nrow(values), ncol(values)),
    j = rep(1:ncol(values), each = nrow(values)),
    x = as.vector(values)
  )
  rownames(sparse_matrix) <- rownames
  colnames(sparse_matrix) <- colnames

  return(sparse_matrix)
}

read_rarefraction_qiime <- function(filename) {
  df_shannon <- data.table::fread(filename)

  # Pivot into long table
  shannon_long <- data.table::melt(data = df_shannon,
                                   measure.vars = colnames(df_shannon)[grepl("depth-", colnames(df_shannon))],
                                   variable.name = "iters",
                                   variable.factor = FALSE,
                                   value.name = "alpha_div")
  # Corrects colnames
  colnames(shannon_long) <- c("SAMPLE-ID", "iters", "alpha_div")

  return(shannon_long)
}
