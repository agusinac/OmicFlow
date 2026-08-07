## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`omics$composition()` -- Argument checks", {
  ## feature_rank and feature_filter are already tested in `omics-feature_merge`
  expect_snapshot(taxa$composition(feature_rank = "Genus", col_name = "1"), error = TRUE)
  expect_snapshot(taxa$composition(feature_rank = "Genus", col_name = "nonexisting"), error = TRUE)

  expect_snapshot(taxa$composition(feature_rank = "Genus", feature_top = "10"), error = TRUE)
  expect_snapshot(taxa$composition(feature_rank = "Genus", feature_top = c(10, 15)), error = TRUE)
  expect_snapshot(taxa$composition(feature_rank = "Genus", feature_top = 16), error = TRUE)
})

test_that("`omics$composition()` -- Behavioral checks", { 
    ## Testing a single column with top 10    
    taxa$reset()
    taxa$scale(method = "tss")
    res <- taxa$composition(
      feature_rank = "Genus",
      col_name = "CONTRAST_sex",
      feature_top = 10
    )
    expect_snapshot(res$data)
    expect_snapshot(res$palette)
    expect_equal(length(res$palette), 11)

    ## Testing a single column with top 15
    taxa$reset()
    taxa$scale(method = "tss")
    res <- taxa$composition(
      feature_rank = "Genus",
      col_name = "CONTRAST_sex",
      feature_top = 15
    )
    expect_snapshot(res$data)
    expect_snapshot(res$palette)
    expect_equal(length(res$palette), 16)

    ## Testing without col_name with top 10
    taxa$reset()
    taxa$scale(method = "tss")
    res <- taxa$composition(
      feature_rank = "Genus",
      feature_top = 10
    )
    expect_snapshot(res$data)
    expect_snapshot(res$palette)
    expect_equal(length(res$palette), 11)

    ## Testing a single column with top 15
    taxa$reset()
    taxa$scale(method = "tss")
    res <- taxa$composition(
      feature_rank = "Genus",
      feature_top = 15
    )
    expect_snapshot(res$data)
    expect_snapshot(res$palette)
    expect_equal(length(res$palette), 16)    
})