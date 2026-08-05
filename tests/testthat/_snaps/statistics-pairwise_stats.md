# `pairwise_adonis()` -- Argument checks

    Code
      pairwise_adonis(x = c(1, 2, 3))
    Condition
      Error in `pairwise_adonis()`:
      ! "x" must be a <dist>

---

    Code
      pairwise_adonis(x = data.frame())
    Condition
      Error in `pairwise_adonis()`:
      ! "x" must be a <dist>

---

    Code
      pairwise_adonis(x = distmat, groups = list("4", "2"))
    Condition
      Error in `pairwise_adonis()`:
      ! "groups" must be a <vector> and not a <list>.

---

    Code
      pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment,
      metadata = matrix())
    Condition
      Error in `pairwise_adonis()`:
      ! "metadata" must be a <data.frame> or <data.table>.

---

    Code
      pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment,
      perm_design = matrix())
    Condition
      Error in `pairwise_adonis()`:
      ! "perm_design" must be a function.

---

    Code
      pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment,
      p.adjust.method = 1)
    Condition
      Error in `pairwise_adonis()`:
      ! "p.adjust.method" must be a character.

---

    Code
      pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = c(
        1, 5))
    Condition
      Error in `pairwise_adonis()`:
      ! "perm" must be a single whole number.

---

    Code
      pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = "1")
    Condition
      Error in `pairwise_adonis()`:
      ! "perm" needs to be a whole number.

---

    Code
      pairwise_adonis(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = 5.2)
    Condition
      Error in `pairwise_adonis()`:
      ! "perm" needs to be a whole number.

# `pairwise_anosim()` -- Argument checks

    Code
      pairwise_anosim(x = c(1, 2, 3))
    Condition
      Error in `pairwise_anosim()`:
      ! "x" must be a <dist>

---

    Code
      pairwise_anosim(x = data.frame())
    Condition
      Error in `pairwise_anosim()`:
      ! "x" must be a <dist>

---

    Code
      pairwise_anosim(x = distmat, groups = list("4", "2"))
    Condition
      Error in `pairwise_anosim()`:
      ! "groups" must be a <vector> and not a <list>.

---

    Code
      pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment,
      metadata = matrix())
    Condition
      Error in `pairwise_anosim()`:
      ! "metadata" must be a <data.frame> or <data.table>.

---

    Code
      pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment,
      perm_design = matrix())
    Condition
      Error in `pairwise_anosim()`:
      ! "perm_design" must be a function.

---

    Code
      pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment,
      p.adjust.method = 1)
    Condition
      Error in `pairwise_anosim()`:
      ! "p.adjust.method" must be a character.

---

    Code
      pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = c(
        1, 5))
    Condition
      Error in `pairwise_anosim()`:
      ! "perm" must be a single whole number.

---

    Code
      pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = "1")
    Condition
      Error in `pairwise_anosim()`:
      ! "perm" needs to be a whole number.

---

    Code
      pairwise_anosim(x = distmat, groups = test$metaData$CONTRAST_treatment, perm = 5.2)
    Condition
      Error in `pairwise_anosim()`:
      ! "perm" needs to be a whole number.

