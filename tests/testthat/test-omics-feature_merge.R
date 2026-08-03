## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`omics$feature_merge()` -- Argument checks", {
  expect_error(taxa$feature_merge(), 'argument "feature_rank" is missing, with no default')
  expect_snapshot(taxa$feature_merge(feature_rank = 1), error = TRUE)
  expect_snapshot(taxa$feature_merge(feature_rank = list()), error = TRUE)
  expect_snapshot(taxa$feature_merge(feature_rank = "1"), error = TRUE)
  expect_snapshot(taxa$feature_merge(feature_rank = c("1", "Kingdom")), error = TRUE)

  expect_snapshot(taxa$feature_merge(feature_rank = "Kingdom", feature_filter = list()), error = TRUE)
  expect_snapshot(taxa$feature_merge(feature_rank = "Kingdom", feature_filter = 1), error = TRUE)
})

test_that("`omics$feature_merge()` -- Behavioral checks", { 
    ## Testing single rank
    taxa$reset()
    features_original <- data.table::copy(taxa$featureData)
    counts_original <- taxa$countData
    
    taxa$feature_merge(
      feature_rank = "Genus",
      feature_filter = c("uncultured")
    )
    expect_snapshot(taxa)
    expect_equal(any(grepl("uncultured", taxa$featureData$Genus)), FALSE)

    expect_equal(
      Matrix::colSums(counts_original[features_original[Genus == "Woesearchaeales"]$FEATURE_ID, ]),
      taxa$countData[taxa$featureData[Genus == "Woesearchaeales"]$FEATURE_ID, ] 
    )

    ## Testing multiple ranks
    taxa$reset()

    taxa$feature_merge(feature_rank = c("Kingdom", "Phylum", "Genus"))
    expect_snapshot(taxa)

    expect_equal(
      Matrix::colSums(counts_original[features_original[Genus == "Woesearchaeales"]$FEATURE_ID, ]),
      taxa$countData[taxa$featureData[Genus == "Woesearchaeales"]$FEATURE_ID, ] 
    )
})