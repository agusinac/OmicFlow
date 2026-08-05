# `plot_pairwise_stats()` -- Argument checks

    Code
      plot_pairwise_stats(data = matrix())
    Condition
      Error in `plot_pairwise_stats()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = 1)
    Condition
      Error in `plot_pairwise_stats()`:
      ! "stats_col" needs to contain characters with length of 1.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = c("1", "2"))
    Condition
      Error in `plot_pairwise_stats()`:
      ! "stats_col" needs to contain characters with length of 1.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = "nothing")
    Condition
      Error in `plot_pairwise_stats()`:
      ! The "stats_col" column does not exist in the provided `data`.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = 1)
    Condition
      Error in `plot_pairwise_stats()`:
      ! "group_col" needs to contain characters with length of 1.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = c("1",
        "2"))
    Condition
      Error in `plot_pairwise_stats()`:
      ! "group_col" needs to contain characters with length of 1.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "nothing")
    Condition
      Error in `plot_pairwise_stats()`:
      ! The "group_col" column does not exist in the provided `data`.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs",
        label_col = 1)
    Condition
      Error in `plot_pairwise_stats()`:
      ! "label_col" needs to contain characters with length of 1.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs",
        label_col = c("1", "2"))
    Condition
      Error in `plot_pairwise_stats()`:
      ! The "label_col" column does not exist in the provided `data`.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs",
        label_col = "nothing")
    Condition
      Error in `plot_pairwise_stats()`:
      ! The "label_col" column does not exist in the provided `data`.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs",
        label_col = "p.adj", y_axis_title = 1)
    Condition
      Error in `plot_pairwise_stats()`:
      ! "y_axis_title" needs to contain characters.

---

    Code
      plot_pairwise_stats(data = adonis_res, stats_col = "F.Model", group_col = "pairs",
        label_col = "p.adj", plot_title = 1)
    Condition
      Error in `plot_pairwise_stats()`:
      ! "plot_title" needs to contain characters.

