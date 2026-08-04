# `diversity()` -- Argument checks

    Code
      diversity(x = data.frame())
    Condition
      Error in `diversity()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      diversity(x = c(2, 1, 1))
    Condition
      Error in `diversity()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      diversity(x = test$countData, normalize = "FALSE")
    Condition
      Error in `diversity()`:
      ! "normalize" needs to be either `TRUE` or `FALSE`.

---

    Code
      diversity(x = test$countData, base = "1")
    Condition
      Error in `diversity()`:
      ! "base" needs to be a <numeric> type with length of 1.

---

    Code
      diversity(x = test$countData, base = c(1, 2))
    Condition
      Error in `diversity()`:
      ! "base" needs to be a <numeric> type with length of 1.

---

    Code
      diversity(x = test$countData, metric = 1)
    Condition
      Error in `diversity()`:
      ! 1 needs to contain characters with length of 1.

---

    Code
      diversity(x = test$countData, metric = c("shannon", "simpson"))
    Condition
      Error in `diversity()`:
      ! "shannon" and "simpson" needs to contain characters with length of 1.

---

    Code
      diversity(x = test$countData)
    Condition
      Error in `diversity()`:
      ! "x" must be non-negative.

# `diversity()` -- Behavioral checks

    Code
      res
    Output
      Sample_01 Sample_02 Sample_03 Sample_04 Sample_05 Sample_06 Sample_07 Sample_08 
       3.902859  3.901478  3.896200  3.899254  3.901116  3.900795  3.905254  3.901093 
      Sample_09 Sample_10 
       3.899304  3.899547 

---

    Code
      res
    Output
      Sample_01 Sample_02 Sample_03 Sample_04 Sample_05 Sample_06 Sample_07 Sample_08 
      0.9796359 0.9795777 0.9793707 0.9794965 0.9795681 0.9795548 0.9797291 0.9795569 
      Sample_09 Sample_10 
      0.9794967 0.9794999 

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
      Sample_01 Sample_02 Sample_03 Sample_04 Sample_05 Sample_06 Sample_07 Sample_08 
       5.630635  5.628643  5.621029  5.625435  5.628120  5.627657  5.634091  5.628087 
      Sample_09 Sample_10 
       5.625506  5.625857 

---

    Code
      res
    Output
      Sample_01 Sample_02 Sample_03 Sample_04 Sample_05 Sample_06 Sample_07 Sample_08 
       1.694990  1.694390  1.692098  1.693425  1.694233  1.694094  1.696030  1.694223 
      Sample_09 Sample_10 
       1.693446  1.693552 

---

    Code
      res
    Output
      Sample_01 Sample_02 Sample_03 Sample_04 Sample_05 Sample_06 Sample_07 Sample_08 
      -4710.670 -4559.788 -4480.961 -4632.565 -4599.439 -4619.721 -4857.336 -4399.274 
      Sample_09 Sample_10 
      -4663.782 -4573.771 

