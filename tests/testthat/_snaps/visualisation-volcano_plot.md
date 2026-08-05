# `volcano_plot()` -- Argument checks

    Code
      volcano_plot(data = data.frame())
    Condition
      Error in `volcano_plot()`:
      ! `data` must be a <data.table>.

---

    Code
      volcano_plot(data = mock_data, logfold_col = 1)
    Condition
      Error in `volcano_plot()`:
      ! "logfold_col" needs to contain characters with length of 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = c("1", "2"))
    Condition
      Error in `volcano_plot()`:
      ! "logfold_col" needs to contain characters with length of 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "nothing")
    Condition
      Error in `volcano_plot()`:
      ! The "logfold_col" column does not exist in the provided `data`.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = 1)
    Condition
      Error in `volcano_plot()`:
      ! "pvalue_col" needs to contain characters with length of 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = c("1", "2"))
    Condition
      Error in `volcano_plot()`:
      ! "pvalue_col" needs to contain characters with length of 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "nothing")
    Condition
      Error in `volcano_plot()`:
      ! The "pvalue_col" column does not exist in the provided `data`.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = 1)
    Condition
      Error in `volcano_plot()`:
      ! "abundance_col" needs to contain characters with length of 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = c("1", "2"))
    Condition
      Error in `volcano_plot()`:
      ! "abundance_col" needs to contain characters with length of 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "nothing")
    Condition
      Error in `volcano_plot()`:
      ! The "abundance_col" column does not exist in the provided `data`.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = 1)
    Condition
      Error in `volcano_plot()`:
      ! "feature_rank" needs to contain characters with length of 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = c("1", "2"))
    Condition
      Error in `volcano_plot()`:
      ! "feature_rank" needs to contain characters with length of 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "nothing")
    Condition
      Error in `volcano_plot()`:
      ! The "feature_rank" column does not exist in the provided `data`.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", pvalue.threshold = "1")
    Condition
      Error in `volcano_plot()`:
      ! "pvalue.threshold" need to be a single numeric value.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", pvalue.threshold = c(1,
          2))
    Condition
      Error in `volcano_plot()`:
      ! "pvalue.threshold" need to be a single numeric value.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", logfold.threshold = "1")
    Condition
      Error in `volcano_plot()`:
      ! "logfold.threshold" need to be a single numeric value.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", logfold.threshold = c(1,
          2))
    Condition
      Error in `volcano_plot()`:
      ! "logfold.threshold" need to be a single numeric value.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", abundance.threshold = "1")
    Condition
      Error in `volcano_plot()`:
      ! "abundance.threshold" need to be a single numeric value.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", abundance.threshold = c(
          1, 2))
    Condition
      Error in `volcano_plot()`:
      ! "abundance.threshold" need to be a single numeric value.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", label_A = 1)
    Condition
      Error in `volcano_plot()`:
      ! "label_A" needs to contain characters of length 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", label_A = c("A1", "A2"))
    Condition
      Error in `volcano_plot()`:
      ! "label_A" needs to contain characters of length 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", label_B = 1)
    Condition
      Error in `volcano_plot()`:
      ! "label_B" needs to contain characters of length 1.

---

    Code
      volcano_plot(data = mock_data, logfold_col = "log2FC", pvalue_col = "pvalue",
        abundance_col = "rel_abun", feature_rank = "Feature", label_B = c("A1", "A2"))
    Condition
      Error in `volcano_plot()`:
      ! "label_B" needs to contain characters of length 1.

