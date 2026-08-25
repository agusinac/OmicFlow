## Load example data
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

colors <- OmicFlow::colormap(dt, "sex")

test_that("`boxjitter_test()` -- Argument checks", {
  expect_snapshot(boxjitter_test(data = list()), error = TRUE)
  expect_snapshot(boxjitter_test(data = matrix()), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt, values = 1), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt, values = c("1", "2")), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt, values = "nonexisting"), error = TRUE)

  expect_snapshot(boxjitter_test(data = dt, values = "shannon", col_name = 1), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt, values = "shannon", col_name = c("1", "2")), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt,, values = "shannon", col_name = "nonexisting"), error = TRUE)

  expect_snapshot(boxjitter_test(data = dt, values = "shannon", col_name = "sex", group_by = 1), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt, values = "shannon", col_name = "sex", group_by = c("1", "2")), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt,, values = "shannon", col_name = "sex", group_by = "nonexisting"), error = TRUE)

  expect_snapshot(boxjitter_test(data = dt, values = "shannon", col_name = "sex", palette = c(1, 2, 3)), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt, values = "shannon", col_name = "sex", palette = c("foo1", "foo2")), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt, values = "shannon", col_name = "sex", palette = colors, method = 1), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt, values = "shannon", col_name = "sex", palette = colors, paired = "FALSE"), error = TRUE)
  expect_snapshot(boxjitter_test(data = dt, values = "shannon", col_name = "sex", palette = colors, p.adjust.method = 5), error = TRUE)
})

test_that("`boxjitter_test()` -- Behavioral checks", { 
  ## Testing default settings
  expect_no_error(res <- boxjitter_test(
    data = dt,
    values = "shannon",
    col_name = "sex",
    palette = colors
  ))
  expect_equal(names(res), c("plot", "stats"))
  expect_s3_class(res$plot, "ggplot")

  ## Testing with paired
  expect_no_error(res <- boxjitter_test(
    data = dt,
    values = "shannon",
    col_name = "sex",
    palette = colors,
    paired = TRUE
  ))
  expect_equal(names(res), c("plot", "stats"))
  expect_s3_class(res$plot, "ggplot")

  ## Testing with `group_by`
  expect_no_error(res <- boxjitter_test(
    data = dt,
    values = "shannon",
    col_name = "sex",
    group_by = "treatment",
    palette = colors
  ))
  expect_equal(names(res), c("plot", "stats"))
  expect_s3_class(res$plot, "ggplot")
})