test_that("Testing transformations of data", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  # Perform log transformation
  taxa$scale(transform = log2)
  expect_snapshot(as.vector(taxa$countData[, 1]))
  
  # Perform sqrt transformation
  taxa$reset()
  taxa$scale(transform = sqrt)
  expect_snapshot(as.vector(taxa$countData[, 1]))
  
  # Perform clr standardisation
  taxa$reset()
  taxa$scale(method = "clr")
  expect_snapshot(as.vector(taxa$countData[, 1]))
  
  # Perform tss normalisation
  taxa$reset()
  taxa$scale(method = "tss")
  expect_snapshot(as.vector(taxa$countData[, 1]))
})