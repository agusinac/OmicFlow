## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`matrix_to_dtable()` -- Argument checks", {
    expect_snapshot(matrix_to_dtable(x = list()), error = TRUE)
    expect_snapshot(matrix_to_dtable(x = data.frame()), error = TRUE)
})

test_that("`matrix_to_dtable()` -- Behavioral checks", {
    dt <- matrix_to_dtable(taxa$countData)
    expect_equal(class(dt), c("data.table", "data.frame"))
    expect_snapshot(dt)
    expect_equal(colnames(dt), taxa$metaData$SAMPLE_ID)
})