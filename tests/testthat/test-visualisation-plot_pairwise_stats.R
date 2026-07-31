## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`plot_pairwise_stats()` -- Argument checks", {

})

test_that("`plot_pairwise_stats()` -- Behavioral checks", { 
 
})

test_that("Testing `plot_pairwise_stats` error-handling", {
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


  expect_error(plot_pairwise_stats(
    data = as.matrix(adonis_res),
    group_col = "pairs",
    stats_col = "F.Model",
    label_col = "p.adj",
    y_axis_title = "Pseudo F test statistic",
    plot_title = "PERMANOVA"
  ))

  expect_error(plot_pairwise_stats(
    data = adonis_res,
    group_col = "pair",
    stats_col = "F.Model",
    label_col = "p.adj",
    y_axis_title = "Pseudo F test statistic",
    plot_title = "PERMANOVA"
  ))

  expect_error(plot_pairwise_stats(
    data = adonis_res,
    group_col = "pairs",
    stats_col = "mod",
    label_col = "p.adj",
    y_axis_title = "Pseudo F test statistic",
    plot_title = "PERMANOVA"
  ))

  expect_error(plot_pairwise_stats(
    data = adonis_res,
    group_col = "pairs",
    stats_col = "F.Model",
    label_col = "p_adj",
    y_axis_title = "Pseudo F test statistic",
    plot_title = "PERMANOVA"
  ))

  expect_error(plot_pairwise_stats(
    data = adonis_res,
    group_col = "pairs",
    stats_col = "F.Model",
    label_col = "p.adj",
    y_axis_title = 5,
    plot_title = "PERMANOVA"
  ))

  expect_error(plot_pairwise_stats(
    data = adonis_res,
    group_col = "pairs",
    stats_col = "F.Model",
    label_col = "p.adj",
    y_axis_title = "Pseudo F test statistic",
    plot_title = 10
  ))
})
