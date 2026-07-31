## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`colormap()` -- Argument checks", {

})

test_that("`colormap()` -- Behavioral checks", { 
 
})

test_that("Testing colormap", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  expect_error(colormap(data = as.matrix(taxa$metaData), col_name = "CONTRAST_sex"))
  expect_error(colormap(data = taxa$metaData))
  expect_error(colormap(data = taxa$metaData, col_name = "sex"))
  expect_error(colormap(data = taxa$metaData, col_name = "CONTRAST_sex", Brewer.palID = 2))
  expect_error(colormap(data = taxa$metaData, col_name = "CONTRAST_sex", Brewer.palID = "colSet"))
})