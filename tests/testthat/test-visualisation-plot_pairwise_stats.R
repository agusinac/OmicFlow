## Load example data
set.seed(42)
mock_data <- matrix(rnorm(15 * 10), nrow = 15, ncol = 10)
mock_dist <- dist(mock_data, method = "euclidean")
mock_groups <- rep(c("A", "B", "C"), each = 5)

# Compute pairwise adonis
adonis_res <- pairwise_adonis(
  x = mock_dist,
  groups = mock_groups, 
  p.adjust.method = "bonferroni", 
  perm = 99)

test_that("`plot_pairwise_stats()` -- Argument checks", {
  expect_snapshot(plot_pairwise_stats(data = matrix()), error = TRUE)
  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = 1), error = TRUE)
  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = c("1", "2")), error = TRUE)
  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = "nothing"), error = TRUE)

  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = 1), error = TRUE)
  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = c("1", "2")), error = TRUE)
  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "nothing"), error = TRUE)

  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs", label_col = 1), error = TRUE)
  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs", label_col = c("1", "2")), error = TRUE)
  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs", label_col = "nothing"), error = TRUE)

  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs", label_col = "p.adj", y_axis_title = 1), error = TRUE)
  expect_snapshot(plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs", label_col = "p.adj", plot_title = 1), error = TRUE)
})

test_that("`plot_pairwise_stats()` -- Behavioral checks", { 
  ## Testing default settings
  expect_no_error(p <- plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs", label_col = "p.adj"))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(res@plot@labels$label, "p.adj")
  expect_equal(res@plot@data[1, ], adonis_res[1, ])

  ## Testing addition of labels
  title_label <- "ANOSIM"
  y_label <- "ANOSIM R statistic"
  expect_no_error(p <- plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs", label_col = "p.adj", y_axis_title = y_label, plot_title = title_label))
  expect_s3_class(p, "ggplot")

  res <- ggplot2::ggplot_build(p)
  expect_equal(res@plot@labels$label, "p.adj")
  expect_equal(res@plot@labels$title, title_label)
  expect_equal(res@plot@labels$y, y_label)
  expect_equal(res@plot@data[1, ], adonis_res[1, ])
})
