df <- data.frame(A = c(5,2,2), B = c(2,2,NA), C = c(NA,NA,NA))
dt <- data.table::data.table(A = c(5,2,2), B = c(2,2,NA), C = c(NA,NA,NA))

test_that("`column_exists()` -- Argument checks", {
    expect_snapshot(column_exists(column = 1), error = TRUE)
    expect_snapshot(column_exists(column = "1", table = matrix()), error = TRUE)
})

test_that("`column_exists()` -- Behavioral checks", {
    ## Testing on data.frame
    expect_equal(column_exists(column = "A", table = df), TRUE)
    expect_equal(column_exists(column = c("A", "B"), table = df), TRUE)
    expect_equal(column_exists(column = "nothing", table = df), FALSE)
    expect_equal(column_exists(column = "C", table = df), FALSE)
    
    ## Testing on data.table
    expect_equal(column_exists(column = "A", table = dt), TRUE)
    expect_equal(column_exists(column = c("A", "B"), table = dt), TRUE)
    expect_equal(column_exists(column = "nothing", table = dt), FALSE)
    expect_equal(column_exists(column = "C", table = dt), FALSE)
 
})