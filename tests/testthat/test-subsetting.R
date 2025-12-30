test_that("Testing subsetting of data", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  # Perform feature subset
  taxa$feature_subset(Kingdom == "Bacteria")
  expect_snapshot(taxa)
  
  # In between reset
  taxa$reset()
  
  # Perform metadata subset
  taxa$sample_subset(treatment == "tumor")
  expect_snapshot(taxa)

  # active binding subset
  taxa$reset()
  taxa$featureData <- taxa$featureData[1:200, ]
  expect_snapshot(taxa)

  taxa$reset()
  taxa$countData <- taxa$countData[1:50, ]
  expect_snapshot(taxa)
})