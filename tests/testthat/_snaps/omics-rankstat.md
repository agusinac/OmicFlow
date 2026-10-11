# `omics$rankstat()` -- Argument checks

    Code
      taxa$rankstat(feature_ranks = list())
    Condition
      Error in `taxa$rankstat()`:
      ! "feature_ranks" needs to contains characters.

---

    Code
      taxa$rankstat(feature_ranks = data.frame())
    Condition
      Error in `taxa$rankstat()`:
      ! "feature_ranks" needs to contains characters.

---

    Code
      taxa$rankstat(feature_ranks = c(1, 2, 4))
    Condition
      Error in `taxa$rankstat()`:
      ! "feature_ranks" needs to contains characters.

---

    Code
      taxa$rankstat(feature_ranks = c("Kingdom", "Genus"), unique = 1)
    Condition
      Error in `taxa$rankstat()`:
      ! "unique" needs to be either `TRUE` or `FALSE`.

---

    Code
      taxa$rankstat(unique = "1")
    Condition
      Error in `taxa$rankstat()`:
      ! argument "feature_ranks" is missing, with no default

