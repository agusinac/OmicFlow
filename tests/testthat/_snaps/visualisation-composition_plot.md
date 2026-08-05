# `composition_plot()` -- Argument checks

    Code
      composition_plot(data = list())
    Condition
      Error in `composition_plot()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      composition_plot(data = matrix())
    Condition
      Error in `composition_plot()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      composition_plot(data = comp$data, palette = list())
    Condition
      Error in `composition_plot()`:
      ! "palette" needs to contain characters.

---

    Code
      composition_plot(data = comp$data, palette = comp$palette, feature_rank = 1)
    Condition
      Error in `composition_plot()`:
      ! "feature_rank" needs to contain characters with length of 1.

---

    Code
      composition_plot(data = comp$data, palette = comp$palette, feature_rank = c("1",
        "2"))
    Condition
      Error in `composition_plot()`:
      ! "feature_rank" needs to contain characters with length of 1.

---

    Code
      composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus",
      title_name = 1)
    Condition
      Error in `composition_plot()`:
      ! "title_name" must be a character and of length 1

---

    Code
      composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus",
      title_name = c("1", "2"))
    Condition
      Error in `composition_plot()`:
      ! "title_name" must be a character and of length 1

---

    Code
      composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus",
      group_by = 1)
    Condition
      Error in `composition_plot()`:
      ! "group_by" must be a character and of length 1

---

    Code
      composition_plot(data = comp$data, palette = comp$palette, feature_rank = "Genus",
      group_by = c("1", "2"))
    Condition
      Error in `composition_plot()`:
      ! "group_by" must be a character and of length 1

