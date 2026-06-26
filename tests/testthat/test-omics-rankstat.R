test_that("Testing Rankstat", {
  # Loading hdf5 format
  taxa_hdf5 <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  expect_no_error(taxa_hdf5$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species")))
  expect_no_error(taxa_hdf5$rankstat(feature_ranks = c("Kingdom", "Phylum", "Family", "Genus", "Species"), unique=TRUE))
})