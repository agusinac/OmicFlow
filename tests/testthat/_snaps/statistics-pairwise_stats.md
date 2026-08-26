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

# `pairwise_test()` -- Argument checks

    Code
      pairwise_test(data = list())
    Condition
      Error in `pairwise_test()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      pairwise_test(data = matrix())
    Condition
      Error in `pairwise_test()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      pairwise_test(data = dt, x_col = 1)
    Condition
      Error in `pairwise_test()`:
      ! "x_col" needs to contain characters with length of 1.

---

    Code
      pairwise_test(data = dt, x_col = c("1", "2"))
    Condition
      Error in `pairwise_test()`:
      ! "x_col" needs to contain characters with length of 1.

---

    Code
      pairwise_test(data = dt, x_col = "nonexisting")
    Condition
      Error in `pairwise_test()`:
      ! The "x_col" column does not exist in the provided `data`.

---

    Code
      pairwise_test(data = dt, x_col = "shannon", g_col = 1)
    Condition
      Error in `pairwise_test()`:
      ! "g_col" needs to contain characters with length of 1.

---

    Code
      pairwise_test(data = dt, x_col = "shannon", g_col = c("1", "2"))
    Condition
      Error in `pairwise_test()`:
      ! "g_col" needs to contain characters with length of 1.

---

    Code
      pairwise_test(data = dt, x_col = "shannon", g_col = "nonexisting")
    Condition
      Error in `pairwise_test()`:
      ! The "g_col" column does not exist in the provided `data`.

---

    Code
      pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", paired = 1)
    Condition
      Error in `pairwise_test()`:
      ! "paired" needs to be either `TRUE` or `FALSE`.

---

    Code
      pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", paired = "FALSE")
    Condition
      Error in `pairwise_test()`:
      ! "paired" needs to be either `TRUE` or `FALSE`.

---

    Code
      pairwise_test(data = dt, x_col = "shannon", g_col = "treatment",
        p.adjust.method = 1)
    Condition
      Error in `pairwise_test()`:
      ! "p.adjust.method" must be a character.

---

    Code
      pairwise_test(data = dt, x_col = "shannon", g_col = "treatment",
        p.adjust.method = "nonexisting")
    Condition
      Error in `pairwise_test()`:
      ! "nonexisting" is not a valid method.  Valid options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none".

---

    Code
      pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = 1)
    Condition
      Error in `pairwise_test()`:
      ! 1 is not a valid test. Valid options: "t.test", "t.equalvar", and "wilcox".

---

    Code
      pairwise_test(data = dt, x_col = "shannon", g_col = "treatment", test = "nonexisting")
    Condition
      Error in `pairwise_test()`:
      ! "nonexisting" is not a valid test. Valid options: "t.test", "t.equalvar", and "wilcox".

