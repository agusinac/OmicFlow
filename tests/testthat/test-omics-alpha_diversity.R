test_that("Testing Alpha diversity", {
  taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
  )
  
  ## Testing Shannon
  res_shannon <- suppressWarnings(
    taxa$alpha_diversity(
      metric = "shannon",
      col_name = "CONTRAST_sex",
    ))

  expect_snapshot(res_shannon$data)
  expect_snapshot(res_shannon$stats)

  # with stratification
  res_grouped_shannon <- suppressWarnings(
    taxa$alpha_diversity(
      metric = "shannon",
      col_name = "CONTRAST_sex",
      group_by = "treatment"
    ))

  expect_snapshot(res_grouped_shannon$data)
  expect_snapshot(res_grouped_shannon$stats)
  
  ## Testing inverse simpson
  res_invsimpson <- suppressWarnings(
    taxa$alpha_diversity(
      metric = "invsimpson",
      col_name = "CONTRAST_sex",
    ))
  
  expect_snapshot(res_invsimpson$data)
  expect_snapshot(res_invsimpson$stats)
  
  ## Testing simpson
  res_simpson <- suppressWarnings(
    taxa$alpha_diversity(
      metric = "simpson",
      col_name = "CONTRAST_sex",
    ))
  
  expect_snapshot(res_simpson$data)
  expect_snapshot(res_simpson$stats)
  
})
