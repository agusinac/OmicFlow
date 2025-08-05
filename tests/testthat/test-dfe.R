test_that("Testing Log2 Foldchanges", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  suppressWarnings(
    dfe <- taxa$DFE(
      feature_rank = "Genus",
      feature_filter = c("uncultured"),
      paired = FALSE,
      condition.group = "CONTRAST_sex",
      condition_A = c("male"),
      condition_B = c("female")
    ))
  
  expect_snapshot(dfe$data)
})