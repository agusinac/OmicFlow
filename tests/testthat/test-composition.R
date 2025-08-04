test_that("Testing Compositional data", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  res <- taxa$composition(
    feature_rank = "Genus",
    feature_filter = c("uncultured"),
    col_name = "CONTRAST_sex",
    feature_top = 10
  )
  
  expect_snapshot(res$data)
  expect_snapshot(res$palette)
})