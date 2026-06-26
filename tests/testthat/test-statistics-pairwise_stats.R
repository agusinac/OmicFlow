test_that("Testing pairwise_adonis & pairwise_anosim", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )

  # Example test data
  set.seed(42)
  mock_data <- matrix(rnorm(15 * 10), nrow = 15, ncol = 10)
  mock_dist <- dist(mock_data, method = "euclidean")
  mock_groups <- rep(c("A", "B", "C"), each = 5)

  expect_error(pairwise_adonis(x = mock_dist, groups = mock_groups, p.adjust.method = "test", perm = 99))
  expect_error(pairwise_adonis(x = as.matrix(mock_dist), groups = mock_groups))
  expect_error(pairwise_adonis(x = mock_dist, groups = list(mock_groups)))
  expect_error(pairwise_adonis(x = mock_dist, groups = mock_groups, perm = 5.5))

  expect_error(pairwise_anosim(x = mock_dist, groups = mock_groups, p.adjust.method = "test", perm = 99))
  expect_error(pairwise_anosim(x = as.matrix(mock_dist), groups = mock_groups))
  expect_error(pairwise_anosim(x = mock_dist, groups = list(mock_groups)))
  expect_error(pairwise_anosim(x = mock_dist, groups = mock_groups, perm = 5.5))
})