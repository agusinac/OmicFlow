## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`omics$alpha_diversity()` -- Argument checks", {
  expect_error(taxa$alpha_diversity(), 'argument "groups" is missing, with no default')
  expect_snapshot(taxa$alpha_diversity(groups = 1), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(groups = "1"), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(groups = c("1", "2")), error = TRUE)
  
  expect_snapshot(taxa$alpha_diversity(groups = "CONTRAST_sex", split_by = 1), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(groups = "CONTRAST_sex", split_by = "1"), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(groups = "CONTRAST_sex", split_by = c("1", "2")), error = TRUE)

  expect_snapshot(taxa$alpha_diversity(groups = "CONTRAST_sex", evenness = "FALSE"), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(groups = "CONTRAST_sex", paired = "FALSE"), error = TRUE)
  expect_snapshot(taxa$alpha_diversity(groups = "CONTRAST_sex", p.adjust.method = "nothing"), error = TRUE)
})

test_that("`omics$alpha_diversity()` -- Behavioral checks", {
  ## Testing with single `groups`
  adiv <- suppressWarnings(taxa$alpha_diversity(groups = "CONTRAST_sex"))
  expect_snapshot(adiv$data)
  expect_equal(round(adiv$stats$pvalue, 2), 0.67)
  expect_s3_class(adiv$plot, "ggplot")
  expect_snapshot(taxa)

  ## Testing with `split_by`
  adiv <- suppressWarnings(taxa$alpha_diversity(groups = "CONTRAST_sex", split_by = "treatment"))
  expect_snapshot(adiv$data)
  expect_equal(round(adiv$stats$pvalue, 2), c(1, 1))
  expect_s3_class(adiv$plot, "ggplot")
})
