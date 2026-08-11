# `diversity_plot()` -- Argument checks

    Code
      diversity_plot(data = list())
    Condition
      Error in `diversity_plot()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      diversity_plot(data = matrix())
    Condition
      Error in `diversity_plot()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      diversity_plot(data = dt, values = 1)
    Condition
      Error in `diversity_plot()`:
      ! "values" needs to contain characters with length of 1.

---

    Code
      diversity_plot(data = dt, values = c("1", "2"))
    Condition
      Error in `diversity_plot()`:
      ! "values" needs to contain characters with length of 1.

---

    Code
      diversity_plot(data = dt, values = "nonexisting")
    Condition
      Error in `diversity_plot()`:
      ! The "values" column does not exist in the provided `data`.

---

    Code
      diversity_plot(data = dt, values = "shannon", col_name = 1)
    Condition
      Error in `diversity_plot()`:
      ! "col_name" needs to contain characters with length of 1.

---

    Code
      diversity_plot(data = dt, values = "shannon", col_name = c("1", "2"))
    Condition
      Error in `diversity_plot()`:
      ! "col_name" needs to contain characters with length of 1.

---

    Code
      diversity_plot(data = dt, , values = "shannon", col_name = "nonexisting")
    Condition
      Error in `diversity_plot()`:
      ! The "col_name" column does not exist in the provided `data`.

---

    Code
      diversity_plot(data = dt, values = "shannon", col_name = "sex", group_by = 1)
    Condition
      Error in `diversity_plot()`:
      ! "group_by" needs to contain characters with length of 1.

---

    Code
      diversity_plot(data = dt, values = "shannon", col_name = "sex", group_by = c(
        "1", "2"))
    Condition
      Error in `diversity_plot()`:
      ! "group_by" needs to contain characters with length of 1.

---

    Code
      diversity_plot(data = dt, , values = "shannon", col_name = "sex", group_by = "nonexisting")
    Condition
      Error in `diversity_plot()`:
      ! The "group_by" column does not exist in the provided `data`.

---

    Code
      diversity_plot(data = dt, values = "shannon", col_name = "sex", palette = c(1,
        2, 3))
    Condition
      Error in `diversity_plot()`:
      ! "palette" needs to contain characters.

---

    Code
      diversity_plot(data = dt, values = "shannon", col_name = "sex", palette = c(
        "foo1", "foo2"))
    Condition
      Error in `diversity_plot()`:
      ! "palette" contains invalid colors.

---

    Code
      diversity_plot(data = dt, values = "shannon", col_name = "sex", palette = colors,
        method = 1)
    Condition
      Error in `diversity_plot()`:
      ! "method" needs to be a character <vector>.

---

    Code
      diversity_plot(data = dt, values = "shannon", col_name = "sex", palette = colors,
        paired = "FALSE")
    Condition
      Error in `diversity_plot()`:
      ! "paired" needs to be either `TRUE` or `FALSE`.

---

    Code
      diversity_plot(data = dt, values = "shannon", col_name = "sex", palette = colors,
        p.adjust.method = 5)
    Condition
      Error in `diversity_plot()`:
      ! 5 is not a valid option.  Valid options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none"

