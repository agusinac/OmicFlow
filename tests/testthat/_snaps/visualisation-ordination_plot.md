# `ordination_plot()` -- Argument checks

    Code
      ordination_plot(data = matrix())
    Condition
      Error in `ordination_plot()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      ordination_plot(data = mock_data, col_name = 1)
    Condition
      Error in `ordination_plot()`:
      ! "col_name" needs to contain characters with length of 1.

---

    Code
      ordination_plot(data = mock_data, col_name = c("1", "2"))
    Condition
      Error in `ordination_plot()`:
      ! "col_name" needs to contain characters with length of 1.

---

    Code
      ordination_plot(data = mock_data, col_name = "nothing")
    Condition
      Error in `ordination_plot()`:
      ! The "col_name" column does not exist in the provided `data`.

---

    Code
      ordination_plot(data = mock_data, col_name = "groups", pair = 1)
    Condition
      Error in `ordination_plot()`:
      ! "pair" needs to be a characters <vector>.

---

    Code
      ordination_plot(data = mock_data, col_name = "groups", pair = c("PC1", "PC2",
        "PC3"))
    Condition
      Error in `ordination_plot()`:
      ! "pair" needs to be a <vector> of length 2.

---

    Code
      ordination_plot(data = mock_data, col_name = "groups", pair = c("PC1", "PC2"),
      dist_explained = c(0.2))
    Condition
      Error in `ordination_plot()`:
      ! "dist_explained" needs to be a <vector> of length 2.

---

    Code
      ordination_plot(data = mock_data, col_name = "groups", pair = c("PC1", "PC2"),
      dist_explained = c("0.2", "0.6"))
    Condition
      Error in `ordination_plot()`:
      ! "dist_explained" needs to be a numeric <vector>.

---

    Code
      ordination_plot(data = mock_data, col_name = "groups", pair = c("PC1", "PC2"),
      dist_explained = c(0.2, 0.5, 0.3))
    Condition
      Error in `ordination_plot()`:
      ! "dist_explained" needs to be a <vector> of length 2.

---

    Code
      ordination_plot(data = mock_data, col_name = "groups", pair = c("PC1", "PC2"),
      dist_metric = 1)
    Condition
      Error in `ordination_plot()`:
      ! "dist_metric" needs to contain characters with length of 1.

---

    Code
      ordination_plot(data = mock_data, col_name = "groups", pair = c("PC1", "PC2"),
      dist_metric = c("1", "2"))
    Condition
      Error in `ordination_plot()`:
      ! "dist_metric" needs to contain characters with length of 1.

