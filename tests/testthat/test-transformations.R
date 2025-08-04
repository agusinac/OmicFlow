test_that("Testing transformations of data", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  # Perform log transformation
  taxa$transform(log2)
  expect_snapshot(taxa)
  
  # In between reset
  taxa$reset()
  
  # Perform sqrt transformation
  taxa$transform(sqrt)
  expect_snapshot(taxa)
})