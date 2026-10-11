# `omics$distance()` -- Argument checks

    Code
      test$distance()
    Condition
      Error in `test$distance()`:
      ! "metric" must be specified!

---

    Code
      test$distance(metric = "nothing")
    Condition
      Error in `test$distance()`:
      ! "nothing" is not a valid "metric".  Valid options: "bray", "jaccard", "cosine", "manhattan", "jsd", "canberra", "unifrac", "euclidean", and "aitchison"

---

    Code
      test$distance(metric = c("n1", "n2"))
    Condition
      Error in `test$distance()`:
      ! "metric" needs to be a character with a length of 1

---

    Code
      test$distance(metric = 1)
    Condition
      Error in `test$distance()`:
      ! "metric" needs to be a character with a length of 1

---

    Code
      test$distance(metric = "unifrac")
    Condition
      Error in `test$distance()`:
      ! The specified "metric" is invalid since no treeData is supplied.

---

    Code
      test$distance(metric = "bray", threads = "1")
    Condition
      Error in `test$distance()`:
      ! "threads" need to be a whole number!

---

    Code
      test$distance(metric = "bray", threads = 50.2)
    Condition
      Error in `test$distance()`:
      ! "threads" need to be a whole number!

# `omics$ordination()` -- Argument checks

    Code
      test$ordination()
    Condition
      Error in `test$ordination()`:
      ! "group_by" must be specified!

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", metric = "nothing")
    Condition
      Error in `self$distance()`:
      ! "nothing" is not a valid "metric".  Valid options: "bray", "jaccard", "cosine", "manhattan", "jsd", "canberra", "unifrac", "euclidean", and "aitchison"

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", metric = c("n1", "n2"))
    Condition
      Error in `self$distance()`:
      ! "metric" needs to be a character with a length of 1

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", metric = 1)
    Condition
      Error in `self$distance()`:
      ! "metric" needs to be a character with a length of 1

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", method = "nothing")
    Condition
      Error in `test$ordination()`:
      ! "nothing" is not a valid "method".  Valid options: "pcoa" and "nmds"

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", method = c("n1", "n2"))
    Condition
      Error in `test$ordination()`:
      ! "method" needs to be a character with a length of 1

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", method = 1)
    Condition
      Error in `test$ordination()`:
      ! "method" needs to be a character with a length of 1

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", perm = "999")
    Condition
      Error in `test$ordination()`:
      ! "perm" need to be a whole number.

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", perm = 50.2)
    Condition
      Error in `test$ordination()`:
      ! "perm" need to be a whole number.

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", threads = "999")
    Condition
      Error in `self$distance()`:
      ! "threads" need to be a whole number!

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", threads = 50.2)
    Condition
      Error in `self$distance()`:
      ! "threads" need to be a whole number!

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", perm_design = list())
    Condition
      Error in `test$ordination()`:
      ! "perm_design" must be a function.

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", perm_design = function(x)
        print(x))
    Output
          SAMPLE_ID SAMPLEPAIR_ID CONTRAST_treatment
             <char>        <char>             <char>
       1: Sample_01          S001              tumor
       2: Sample_02          S002              tumor
       3: Sample_03          S003              tumor
       4: Sample_04          S004              tumor
       5: Sample_05          S005              tumor
       6: Sample_06          S001            control
       7: Sample_07          S002            control
       8: Sample_08          S010            control
       9: Sample_09          S004            control
      10: Sample_10          S012            control
    Condition
      Error in `Math.data.frame()`:
      ! non-numeric-alike variable(s) in data frame: SAMPLE_ID, SAMPLEPAIR_ID, CONTRAST_treatment

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", distmat = list())
    Condition
      Error in `test$ordination()`:
      ! "distmat" need to be <Matrix> or <dist>

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", distmat = dist(c(2, 1, 2)))
    Condition
      Error in `test$ordination()`:
      ! None "SAMPLE_ID" from metaData match the "distmat" colnames!

