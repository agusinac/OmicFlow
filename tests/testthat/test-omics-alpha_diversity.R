## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`omics$alpha_diversity()` -- Argument checks", {
  expect_error(taxa$alpha_diversity(), 'argument "col_name" is missing, with no default')
  expect_snapshot(taxa$alpha_diversity(col_name = 1), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(col_name = "1"), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(col_name = c("1", "2")), error = TRUE)
  
  expect_snapshot(taxa$alpha_diversity(col_name = "CONTRAST_sex", group_by = 1), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(col_name = "CONTRAST_sex", group_by = "1"), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(col_name = "CONTRAST_sex", group_by = c("1", "2")), error = TRUE)

  expect_snapshot(taxa$alpha_diversity(col_name = "CONTRAST_sex", evenness = "FALSE"), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(col_name = "CONTRAST_sex", paired = "FALSE"), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(col_name = "CONTRAST_sex", p.adjust.method = "nothing"), error = TRUE)
})

test_that("`omics$alpha_diversity()` -- Behavioral checks", {
  ## Testing with single `col_name`
  adiv <- suppressWarnings(taxa$alpha_diversity(col_name = "CONTRAST_sex"))
  expect_snapshot(adiv$data)
  expect_equal(round(adiv$stats$p, 2), 0.67)
  expect_s3_class(adiv$plot, "ggplot")
  expect_snapshot(taxa)

  ## Testing with `group_by`
  adiv <- suppressWarnings(taxa$alpha_diversity(col_name = "CONTRAST_sex", group_by = "treatment"))
  expect_snapshot(adiv$data)
  expect_equal(round(adiv$stats$p, 2), c(1, 1))
  expect_s3_class(adiv$plot, "ggplot")
})
