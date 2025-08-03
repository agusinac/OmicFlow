library("OmicFlow")

taxa <- metataxonomics$new(metaData = "input/metataxonomics/metadata.tsv",
                           biomData = "input/metataxonomics/biom_with_taxonomy.biom",
                           treeData = "input/metataxonomics/rooted_tree.newick")

test_that("metataxonomics data wrangling functions", {
  # Checking class types to be correct
  testthat::expect_equal(class(taxa), c("metataxonomics", "tools", "R6"))
  testthat::expect_equal(class(taxa$featureData), c("data.table", "data.frame"))
  testthat::expect_equal(class(taxa$metaData), c("data.table", "data.frame"))
  testthat::expect_equal(class(taxa$countData), c("data.table", "data.frame"))
  testthat::expect_equal(class(taxa$treeData), c("phylo"))

  # Checking the correct dimensions
  testthat::expect_equal(dim(taxa$featureData), c(242, 8))
  testthat::expect_equal(dim(taxa$metaData), c(4, 6))
  testthat::expect_equal(dim(taxa$countData), c(242, 4))
  testthat::expect_equal(length(taxa$treeData$tip.label), 242)
  testthat::expect_equal(taxa$treeData$Nnode, 241)

  # Checking data wrangling utils
  # Feature subset
  taxa$feature_subset(Domain == "Bacteria")
  testthat::expect_equal(dim(taxa$featureData), c(185, 8))
  testthat::expect_equal(dim(taxa$metaData), c(4, 6))
  testthat::expect_equal(dim(taxa$countData), c(185, 4))
  testthat::expect_equal(length(taxa$treeData$tip.label), 185)
  testthat::expect_equal(taxa$treeData$Nnode, 184)

  # Sample subset
  taxa$sample_subset(BMI %in% c(22, 23))
  testthat::expect_equal(dim(taxa$featureData), c(133, 8))
  testthat::expect_equal(dim(taxa$metaData), c(3, 6))
  testthat::expect_equal(dim(taxa$countData), c(133, 3))
  testthat::expect_equal(length(taxa$treeData$tip.label), 133)
  testthat::expect_equal(taxa$treeData$Nnode, 132)

  # Feature agglomerate
  taxa$feature_glom(feature_rank = "Genus")
  testthat::expect_equal(dim(taxa$featureData), c(42, 8))
  testthat::expect_equal(dim(taxa$metaData), c(3, 6))
  testthat::expect_equal(dim(taxa$countData), c(42, 3))
  testthat::expect_equal(length(taxa$treeData$tip.label), 42)
  testthat::expect_equal(taxa$treeData$Nnode, 41)

  # column Sums prior to transformation
  testthat::expect_equal(as.numeric(colSums(taxa$countData)), c(964, 629, 707))

  # Transformation: Custom function
  taxa$transform(function(x) x / sum(x))
  testthat::expect_equal(as.numeric(colSums(taxa$countData)), c(1, 1, 1))

  # Transformation: inbuilt function
  taxa$transform(log2)
  testthat::expect_equal(as.numeric(colSums(taxa$countData)), c(-120.55790, -84.23311, -120.82340))

  # Reset all above changes
  taxa$reset()
  testthat::expect_equal(dim(taxa$featureData), c(242, 8))
  testthat::expect_equal(dim(taxa$metaData), c(4, 6))
  testthat::expect_equal(dim(taxa$countData), c(242, 4))
  testthat::expect_equal(length(taxa$treeData$tip.label), 242)
  testthat::expect_equal(taxa$treeData$Nnode, 241)
})
