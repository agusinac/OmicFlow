## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`omics$rankstat()` -- Argument checks", {
  expect_snapshot(taxa$rankstat(feature_ranks = list()), error = TRUE)
  expect_snapshot(taxa$rankstat(feature_ranks = data.frame()), error = TRUE)
  expect_snapshot(taxa$rankstat(feature_ranks = c(1,2,4)), error = TRUE)
  expect_snapshot(taxa$rankstat(feature_ranks = c("Kingdom", "Genus"), unique = 1), error = TRUE)
  expect_snapshot(taxa$rankstat(unique = "1"), error = TRUE)
})

test_that("`omics$rankstat()` -- Behavioral checks", {
  taxa_rank_5 <- c("Kingdom", "Phylum", "Family", "Genus", "Species")
  taxa_rank_3 <- c("Family", "Genus", "Species")

  expect_no_error(p <- taxa$rankstat(feature_ranks = taxa_rank_5))
  data <- ggplot_build(p)
  expect_s3_class(p, "ggplot")
  expect_equal(length(data@data[[1]]$x), length(taxa_rank_5))
  expect_no_error(taxa$rankstat(feature_ranks = taxa_rank_5, unique = TRUE))

  expect_no_error(p <- taxa$rankstat(feature_ranks = taxa_rank_3))
  expect_s3_class(p, "ggplot")
  data <- ggplot_build(p)
  expect_equal(length(data@data[[1]]$x), length(taxa_rank_3))
  expect_no_error(taxa$rankstat(feature_ranks = taxa_rank_3, unique = TRUE))
})