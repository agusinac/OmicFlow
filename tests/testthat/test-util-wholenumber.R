test_that("`is.wholenumber()` -- Behavioral checks", {
    expect_equal(is.wholenumber(5.1), FALSE)
    expect_equal(is.wholenumber(5), TRUE)
})