## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`colormap()` -- Argument checks", {
  expect_snapshot(colormap(data = as.matrix(taxa$metaData)), error = TRUE)
  expect_snapshot(colormap(data = data.frame(), col_name = "nonexisting"), error = TRUE)

  expect_snapshot(colormap(data = data.frame(), col_name = c("col1", "col2")), error = TRUE)
  expect_snapshot(colormap(data = taxa$metaData, col_name = "CONTRAST_sex", Brewer.palID = 2), error = TRUE)
  expect_snapshot(colormap(data = taxa$metaData, col_name = "CONTRAST_sex", Brewer.palID = "colSet"), error = TRUE)
})

test_that("`colormap()` -- Behavioral checks", {
  col1 <- colormap(data = taxa$metaData, col_name = "CONTRAST_sex")

  expect_equal(names(col1), unique(taxa$metaData$CONTRAST_sex))
  expect_equal(rownames(col2rgb(col1)), c("red", "green", "blue"))

  col2 <- colormap(data = taxa$metaData, col_name = "SAMPLE_ID", Brewer.palID = "Oranges")

  expect_equal(names(col2), unique(taxa$metaData$SAMPLE_ID))
  expect_equal(rownames(col2rgb(col2)), c("red", "green", "blue"))
})