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

read_sparseTable <- function(filename) {
  # Read text file, supports csv, excel and tsv formats
  dt <- data.table::fread(filename)
  dt[, (names(dt)) := lapply(.SD, function(x) {
    x <- gsub("\\s+", "", x)                      # Removes spaces between strings
    x <- gsub("^[A-Za-z]*", "", x)                # Removes letters
    })]

  # Convert to matrix format
  mat_1 <- as.matrix(dt,
                     rownames = rownames(dt),
                     colnames = colnames(dt))

  # Change character values to numeric
  mat_2 <- matrix(data = as.numeric(mat_1),
                  ncol = ncol(dt))
  colnames(mat_2) <- colnames(dt)

  mat_2[is.na(mat_2) | mat_2 == ""] <- 0          # Empty strings from cleaning step

  # Return sparseMatrix
  return(as(mat_2, "sparseMatrix"))
}

find_pairs <- function(dt_A, dt_B, unique.ids) {
  result <- data.table::rbindlist(lapply(unique.ids, function(id) {

    A_matches <- colnames(dt_A)[stringr::str_detect(colnames(dt_A), id)]
    B_matches <- colnames(dt_B)[stringr::str_detect(colnames(dt_B), id)]

    if (length(A_matches) == 1 && length(B_matches) == 1) {
      data.table::data.table(colnames_A = A_matches, colnames_B = B_matches, id = id)
    } else {
      NULL
    }

  }))
  result[, id := NULL]

  return(result)
}
