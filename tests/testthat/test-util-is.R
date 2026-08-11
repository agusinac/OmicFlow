test_that("`is.wholenumber()` -- Behavioral checks", {
    expect_equal(is.wholenumber(5.1), FALSE)
    expect_equal(is.wholenumber(5), TRUE)
})

test_that("`is.color()` -- Behavioral checks", {
    expect_equal(is.color("foo"), FALSE)
    expect_equal(is.color("black"), TRUE)

    colors <- stats::setNames(c("#F54927", "#55DE1B"), c("t1", "t2"))
    expect_equal(is.color(colors), TRUE)
})