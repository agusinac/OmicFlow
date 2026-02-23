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

    # Testing `matrix_to_dtable`
    expect_no_error(matrix_to_dtable(taxa$countData))
    expect_error(matrix_to_dtable(dt_sparse))
})

test_that("Testing non-exported utils functions", {
    condition1 <- data.frame(
        group1 = c("A", "A", "B"),
        group2 = c("B", "C", "C")
    )

    condition2 <- data.frame(
        group1 = c("B", "D", "C", "E"),
        group2 = c("A", "B", "D", "A")
    )

    expect_no_error(combine_conditions(condition1, condition2))
    expect_error(combine_conditions(as.matrix(condition1), condition2))
})