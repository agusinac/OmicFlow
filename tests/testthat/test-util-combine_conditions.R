test_that("`combine_conditions()` -- Argument checks", {
  expect_snapshot(combine_conditions(condition1 = matrix()), error = TRUE)
  expect_snapshot(combine_conditions(condition1 = list()), error = TRUE)
  
  expect_snapshot(combine_conditions(condition2 = matrix()), error = TRUE)
  expect_snapshot(combine_conditions(condition2 = list()), error = TRUE)

  expect_snapshot(combine_conditions(condition1 = matrix(), condition2 = data.frame()), error = TRUE)
  expect_snapshot(combine_conditions(condition1 = data.frame(), condition2 = matrix()), error = TRUE)
})

test_that("`combine_conditions()` -- Behavioral checks", {
  condition1 <- data.table::data.table(
    group1 = c("A", "B", "C"),
    group2 = c("B", "C", "D"),
    value  = c(1, 2, 3)
  )

  condition2 <- data.table::data.table(
    group1 = c("B", "C", "D", "E"),
    group2 = c("A", "D", "E", "F"),
    value  = c(10, 20, 30, 40)
  )

  out <- combine_conditions(condition1, condition2)

  expect_s3_class(out, "data.table")
  expect_equal(nrow(out), 5)

  expected <- rbind(
    condition1,
    condition2[3:4]   # D_E and E_F
  )

  expect_equal(out, expected)

  # Check that reversed duplicate pairs are treated as duplicates
  expect_false(any(paste(pmin(out$group1, out$group2), pmax(out$group1, out$group2), sep = "_") %in%
                     c("A_B", "B_C", "C_D", "D_E", "E_F") == FALSE))
})