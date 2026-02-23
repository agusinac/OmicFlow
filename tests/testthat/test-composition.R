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

  # Testing composition_plot
  expect_no_error(composition_plot(
    data = res$data, 
    palette = res$palette, 
    feature_rank = "Genus", 
    group_by = "CONTRAST_sex", 
    title_name = "Awesome title"
  ))
  expect_error(composition_plot(data = res$data, palette = res$palette, feature_rank = "Species"))
  expect_error(composition_plot(data = res$data, palette = res$palette, feature_rank = "Genus", group_by = "treatment"))
  expect_error(composition_plot(data = res$data, palette = res$palette, feature_rank = "Genus", group_by = 1))
  expect_error(composition_plot(data = res$data, palette = c(1,2,3,4), feature_rank = "Genus"))
  expect_error(composition_plot(data = as.matrix(res$data), palette = res$palette, feature_rank = "Genus"))
})