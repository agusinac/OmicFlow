test_that("Testing Log2 Foldchanges", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  suppressWarnings(
    dfe_ungrouped <- taxa$foldchange(
      feature_rank = "Genus",
      feature_filter = c("uncultured"),
      paired = FALSE,
      normalize = TRUE,
      condition.group = "CONTRAST_sex",
      condition_A = c("male"),
      condition_B = c("female")
    ))

  suppressWarnings(
    dfe_grouped <- taxa$foldchange(
      feature_rank = "Genus",
      feature_filter = c("uncultured"),
      paired = FALSE,
      normalize = TRUE,
      group_by = "CONTRAST_sex",
      condition.group = "treatment",
      condition_A = c("tumor"),
      condition_B = c("healthy")
    ))

  skip_if(grepl("devel", R.version$status)) # Due to numerical differences in p-value
  expect_snapshot(dfe_ungrouped$data)
  expect_snapshot(dfe_grouped$data)
})