# `hill_taxa()` -- Argument checks

    Code
      hill_taxa(x = data.frame())
    Condition
      Error in `hill_taxa()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      hill_taxa(x = c(2, 1, 1))
    Condition
      Error in `hill_taxa()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      hill_taxa(x = test$countData, normalize = "FALSE")
    Condition
      Error in `hill_taxa()`:
      ! "normalize" needs to be either `TRUE` or `FALSE`.

---

    Code
      hill_taxa(x = test$countData, base = "1")
    Condition
      Error in `hill_taxa()`:
      ! "base" needs to be a <numeric> type with length of 1.

---

    Code
      hill_taxa(x = test$countData, base = c(1, 2))
    Condition
      Error in `hill_taxa()`:
      ! "base" needs to be a <numeric> type with length of 1.

---

    Code
      hill_taxa(x = test$countData, q = "1")
    Condition
      Error in `hill_taxa()`:
      ! "q" is not a whole number.

---

    Code
      hill_taxa(x = test$countData, q = c(0, 1))
    Condition
      Error in `hill_taxa()`:
      ! "q" needs to be a single whole number.

---

    Code
      hill_taxa(x = test$countData)
    Condition
      Error in `hill_taxa()`:
      ! "x" must be non-negative.

# `hill_taxa()` -- Behavioral checks

    Code
      res
    Output
       [1] 50 50 50 50 50 50 50 50 50 50

---

    Code
      res
    Output
      Sample_01 Sample_02 Sample_03 Sample_04 Sample_05 Sample_06 Sample_07 Sample_08 
       49.54390  49.47553  49.21509  49.36563  49.45759  49.44173  49.66269  49.45646 
      Sample_09 Sample_10 
       49.36807  49.38007 

---

    Code
      res
    Output
      Sample_01 Sample_02 Sample_03 Sample_04 Sample_05 Sample_06 Sample_07 Sample_08 
       49.10609  48.96610  48.47464  48.77214  48.94309  48.91131  49.33188  48.91635 
      Sample_09 Sample_10 
       48.77275  48.78022 

---

    Code
      res
    Output
       [1] 50 50 50 50 50 50 50 50 50 50

---

    Code
      res
    Output
       [1] 50 50 50 50 50 50 50 50 50 50

---

    Code
      res
    Output
       [1] 50 50 50 50 50 50 50 50 50 50

