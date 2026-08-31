# `omics$autoFlow()` -- Argument checks

    Code
      taxa$autoFlow(feature_contrast = 1)
    Condition
      Error in `taxa$autoFlow()`:
      ! "feature_contrast" needs to contain characters.

---

    Code
      taxa$autoFlow(feature_contrast = c("2", "3"))
    Condition
      Error in `taxa$autoFlow()`:
      ! "2" and "3" does not exist in featureData!

---

    Code
      taxa$autoFlow(distmat = "path/nonexisting")
    Condition
      Error in `taxa$autoFlow()`:
      ! Error in "distmat":
      x Input must be an existing "filepath", non-empty <matrix>, <Matrix>, <data.frame> or <data.table>.

---

    Code
      taxa$autoFlow(distmat = matrix(0, 0, 0))
    Condition
      Error in `taxa$autoFlow()`:
      ! Error in "distmat":
      x Input must be an existing "filepath", non-empty <matrix>, <Matrix>, <data.frame> or <data.table>.

---

    Code
      taxa$autoFlow(report = "FALSE")
    Condition
      Error in `taxa$autoFlow()`:
      ! "report" needs to be either `TRUE` or `FALSE`.

---

    Code
      taxa$autoFlow(filename = c("path/1", "path/2"))
    Condition
      Error in `taxa$autoFlow()`:
      ! "filename" needs to be a character with a length of 1

---

    Code
      taxa$autoFlow()
    Condition
      Error in `taxa$autoFlow()`:
      ! No columns with prefix "CONTRAST" found.. Did you forgot to add a prefix?

