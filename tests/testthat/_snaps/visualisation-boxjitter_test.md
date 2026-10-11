# `boxjitter_test()` -- Argument checks

    Code
      boxjitter_test(data = list())
    Condition
      Error in `boxjitter_test()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      boxjitter_test(data = matrix())
    Condition
      Error in `boxjitter_test()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      boxjitter_test(data = dt, values = 1)
    Condition
      Error in `boxjitter_test()`:
      ! "values" needs to contain characters with length of 1.

---

    Code
      boxjitter_test(data = dt, values = c("1", "2"))
    Condition
      Error in `boxjitter_test()`:
      ! "values" needs to contain characters with length of 1.

---

    Code
      boxjitter_test(data = dt, values = "nonexisting")
    Condition
      Error in `boxjitter_test()`:
      ! The "values" column does not exist in the provided `data`.

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = 1)
    Condition
      Error in `boxjitter_test()`:
      ! "groups" needs to contain characters with length of 1.

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = c("1", "2"))
    Condition
      Error in `boxjitter_test()`:
      ! "groups" needs to contain characters with length of 1.

---

    Code
      boxjitter_test(data = dt, , values = "shannon", groups = "nonexisting")
    Condition
      Error in `boxjitter_test()`:
      ! The "groups" column does not exist in the provided `data`.

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = "sex", split_by = 1)
    Condition
      Error in `boxjitter_test()`:
      ! "split_by" needs to contain characters with length of 1.

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = "sex", split_by = c("1",
        "2"))
    Condition
      Error in `boxjitter_test()`:
      ! "split_by" needs to contain characters with length of 1.

---

    Code
      boxjitter_test(data = dt, , values = "shannon", groups = "sex", split_by = "nonexisting")
    Condition
      Error in `boxjitter_test()`:
      ! The "split_by" column does not exist in the provided `data`.

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = "sex", palette = c(1, 2,
        3))
    Condition
      Error in `boxjitter_test()`:
      ! "palette" needs to contain characters.

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = "sex", palette = c(
        "foo1", "foo2"))
    Condition
      Error in `boxjitter_test()`:
      ! "palette" contains invalid colors.

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = "sex", palette = colors,
        test = 1)
    Condition
      Error in `pairwise_test()`:
      ! 1 is not a valid test. Valid options: "t.test", "t.equalvar", and "wilcox".

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = "sex", palette = colors,
        test = "none")
    Condition
      Error in `pairwise_test()`:
      ! "none" is not a valid test. Valid options: "t.test", "t.equalvar", and "wilcox".

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = "sex", palette = colors,
        paired = "FALSE")
    Condition
      Error in `boxjitter_test()`:
      ! "paired" needs to be either `TRUE` or `FALSE`.

---

    Code
      boxjitter_test(data = dt, values = "shannon", groups = "sex", palette = colors,
        p.adjust.method = 5)
    Condition
      Error in `boxjitter_test()`:
      ! 5 is not a valid option.  Valid options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none"

