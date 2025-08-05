test_that("Testing Metagenomics reading and writing of BIOM files", {
  taxa_hdf5 <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  taxa_ref <- metagenomics$new(
    countData = taxa_hdf5$countData,
    metaData = taxa_hdf5$metaData,
    treeData = taxa_hdf5$treeData,
    featureData = taxa_hdf5$featureData
  )
  
  expect_snapshot(taxa_hdf5)
  expect_snapshot(taxa_ref)
})