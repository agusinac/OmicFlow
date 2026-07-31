## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`ordination_plot()` -- Argument checks", {

})

test_that("`ordination_plot()` -- Behavioral checks", { 
 
})

test_that("Testing `ordination_plot` error-handling", {
  set.seed(123)
  mock_data <- data.frame(
    SampleID = paste0("Sample", 1:10),
    PC1 = rnorm(10, mean = 0, sd = 1),
    PC2 = rnorm(10, mean = 0, sd = 1),
    groups = rep(c("Group1", "Group2"), each = 5)
  )

  expect_error(ordination_plot(
    data = as.matrix(mock_data),
    col_name = "groups",
    pair = c("PC1", "PC2"),
    dist_explained = c(45, 22),
    dist_metric = "bray-curtis"
  ))
  expect_error(ordination_plot(
    data = mock_data,
    col_name = "group",
    pair = c("PC1", "PC2"),
    dist_explained = c(45, 22),
    dist_metric = "bray-curtis"
  ))
  expect_error(ordination_plot(
    data = mock_data,
    col_name = "groups",
    pair = c("PC1"),
    dist_explained = c(45, 22),
    dist_metric = "bray-curtis"
  ))
  expect_error(ordination_plot(
    data = mock_data,
    col_name = "groups",
    pair = c("PC1", "PC2"),
    dist_explained = c(45),
    dist_metric = "bray-curtis"
  ))
})
