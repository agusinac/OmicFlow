## Load example data
set.seed(123)
mock_data <- data.frame(
  SampleID = paste0("Sample", 1:10),
  PC1 = rnorm(10, mean = 0, sd = 1),
  PC2 = rnorm(10, mean = 0, sd = 1),
  groups = rep(c("Group1", "Group2"), each = 5)
)

test_that("`ordination_plot()` -- Argument checks", {
  expect_snapshot(ordination_plot(data = matrix()), error = TRUE)
  expect_snapshot(ordination_plot(data = mock_data, groups = 1), error = TRUE)
  expect_snapshot(ordination_plot(data = mock_data, groups = c("1", "2")), error = TRUE)
  expect_snapshot(ordination_plot(data = mock_data, groups = "nothing"), error = TRUE)

  expect_snapshot(ordination_plot(data = mock_data, groups = "groups", pair = 1), error = TRUE)
  expect_snapshot(ordination_plot(data = mock_data, groups = "groups", pair = c("PC1", "PC2", "PC3")), error = TRUE)

  expect_snapshot(ordination_plot(data = mock_data, groups = "groups", pair = c("PC1", "PC2"), dist_explained = c(0.2)), error = TRUE)
  expect_snapshot(ordination_plot(data = mock_data, groups = "groups", pair = c("PC1", "PC2"), dist_explained = c("0.2", "0.6")), error = TRUE)
  expect_snapshot(ordination_plot(data = mock_data, groups = "groups", pair = c("PC1", "PC2"), dist_explained = c(0.2, 0.5, 0.3)), error = TRUE)

  expect_snapshot(ordination_plot(data = mock_data, groups = "groups", pair = c("PC1", "PC2"), dist_metric = 1), error = TRUE)
  expect_snapshot(ordination_plot(data = mock_data, groups = "groups", pair = c("PC1", "PC2"), dist_metric = c("1", "2")), error = TRUE)
})

test_that("`ordination_plot()` -- Behavioral checks", { 
  ## Testing default settings
  expect_no_error(p <- ordination_plot(data = mock_data, groups = "groups", pair = c("PC1", "PC2")))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(res@plot@labels$x, "PC1")
  expect_equal(res@plot@labels$y, "PC2")
  expect_equal(res@plot@labels$colour, "groups")
  expect_equal(res@plot@data[1, 2:3], mock_data[1, 2:3])

  ## Testing with label `dist_metric`
  expect_no_error(p <- ordination_plot(data = mock_data, groups = "groups", pair = c("PC1", "PC2"), dist_metric = "bray"))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(res@plot@labels$x, "PC1")
  expect_equal(res@plot@labels$y, "PC2")
  expect_equal(res@plot@labels$title, "Distance metric used: bray")
  expect_equal(res@plot@labels$colour, "groups")
  expect_equal(res@plot@data[1, 2:3], mock_data[1, 2:3])

  ## Testing with label `dist_explained`
  expect_no_error(p <- ordination_plot(data = mock_data, groups = "groups", pair = c("PC1", "PC2"), dist_explained = c(60, 20)))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(res@plot@labels$x, "PC1 (60%)")
  expect_equal(res@plot@labels$y, "PC2 (20%)")
  expect_equal(res@plot@labels$colour, "groups")
  expect_equal(res@plot@data[1, 2:3], mock_data[1, 2:3])
})