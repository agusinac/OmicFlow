## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`diversity_plot()` -- Argument checks", {

})

test_that("`diversity_plot()` -- Behavioral checks", { 
 
})

test_that("Testing `diversity_plot` error-handling", {  
  # Testing Shannon
  n_row <- 1000
  n_col <- 100
  density <- 0.2
  num_entries <- n_row * n_col
  num_nonzero <- round(num_entries * density)

  set.seed(123)
  positions <- sample(num_entries, num_nonzero, replace=FALSE)
  row_idx <- ((positions - 1) %% n_row) + 1
  col_idx <- ((positions - 1) %/% n_row) + 1

  values <- runif(num_nonzero, min = 0, max = 1)
  sparse_mat <- Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = values,
    dims = c(n_row, n_col)
  )

  div <- OmicFlow::diversity(
  x = sparse_mat,
  metric = "shannon"
  )

  dt <- data.table::data.table(
  "shannon" = div,
  "treatment" = c(rep("healthy", n_col / 2), rep("tumor", n_col / 2)),
  "sex" = c(rep("male", n_col / 4), rep("female", n_col / 4))
  )

  colors <- OmicFlow::colormap(dt, "treatment")

  expect_error(OmicFlow::diversity_plot(
    data = as.matrix(dt),
    values = "shannon",
    col_name = "treatment",
    palette = colors,
    method = "shannon"
  ))
  expect_error(OmicFlow::diversity_plot(
    data = dt,
    values = "val",
    col_name = "treatment",
    palette = colors,
    method = "shannon"
  ))
  expect_error(OmicFlow::diversity_plot(
    data = dt,
    values = "shannon",
    col_name = "treat",
    palette = colors,
    method = "shannon"
  ))
  expect_error(OmicFlow::diversity_plot(
    data = dt,
    values = "shannon",
    col_name = "treatment",
    palette = colors,
    method = 5
  ))
  expect_no_error(OmicFlow::diversity_plot(
    data = dt,
    values = "shannon",
    col_name = "treatment",
    palette = colors,
    method = "shannon",
    group_by = "sex"
  ))
})
