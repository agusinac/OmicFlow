test_that("Testing utils functions", {
    ## Testing column_exists
    taxa <- metagenomics$new(
        biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
        metaData = "input/metagenomics/metadata.tsv",
        treeData = "input/metagenomics/rooted_tree.newick"
    )

    expect_true(column_exists("CONTRAST_sex", taxa$metaData))
    expect_false(column_exists("features", taxa$metaData))

    ## Testing sparse_to_dtable
    sparsemat <- Matrix::sparseMatrix(
        i = c(1, 3, 4, 5), 
        j = c(2, 1, 4, 3), 
        x = c(10, 20, 30, 40), 
        dims = c(5, 4)
        )

    dt_sparse <- matrix_to_dtable(sparsemat)
    expect_snapshot(dt_sparse)
})