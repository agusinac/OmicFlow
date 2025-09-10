test_that("Testing utils functions", {
    ## Testing read_rarefraction_qiime

    mock_df <- data.frame(
        SAMPLE_ID = c("SampleA", "SampleB", "SampleC"),
        `depth-1` = c(4.01, 3.95, NA),
        `depth-2` = c(4.22, 4.06, 4.32),
        `depth-3` = c(4.31, NA, 4.41),
        check.names = FALSE
    )
    data.table::fwrite(mock_df, "mock_rarefaction.txt", sep = "\t")
    result <- read_rarefraction_qiime("mock_rarefaction.txt")

    # check & cleanup
    expect_snapshot(result)
    file.remove("mock_rarefaction.txt")

    ## Testing column_exists
    expect_true(column_exists("iters", result))
    expect_false(column_exists("features", result))

    ## Testing sparse_to_dtable
    sparsemat <- Matrix::sparseMatrix(
        i = c(1, 3, 4, 5), 
        j = c(2, 1, 4, 3), 
        x = c(10, 20, 30, 40), 
        dims = c(5, 4)
        )

    dt_sparse <- sparse_to_dtable(sparsemat)
    expect_snapshot(dt_sparse)
})