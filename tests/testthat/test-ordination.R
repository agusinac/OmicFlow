test_that("Testing UniFrac ordination", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  res <- taxa$ordination(
    metric = "unifrac",
    method = "pcoa",
    group_by = "CONTRAST_sex",
    weighted = TRUE,
    normalize = TRUE,
    cpus = 1
  )
  
  expect_snapshot(res$anova_data)
  expect_snapshot(res$dist)
  expect_snapshot(res$pcs)
})