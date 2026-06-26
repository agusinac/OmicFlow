test_that("Testing exported utils functions", {
    # Load test data
    taxa <- metagenomics$new(
        biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
        metaData = "input/metagenomics/metadata.tsv",
        treeData = "input/metagenomics/rooted_tree.newick"
    )

    # Testing `column_exists`
    expect_true(column_exists("CONTRAST_sex", taxa$metaData))
    expect_false(column_exists("features", taxa$metaData))
})