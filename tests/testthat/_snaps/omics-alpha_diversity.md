# `omics$alpha_diversity()` -- Argument checks

    Code
      taxa$alpha_diversity(col_name = 1)
    Condition
      Error in `taxa$alpha_diversity()`:
      ! 1 must be a character and of length 1

---

    Code
      taxa$alpha_diversity(col_name = "1")
    Condition
      Error in `taxa$alpha_diversity()`:
      ! The specified "1" does not exist in the metaData.

---

    Code
      taxa$alpha_diversity(col_name = c("1", "2"))
    Condition
      Error in `taxa$alpha_diversity()`:
      ! "1" and "2" must be a character and of length 1

---

    Code
      taxa$alpha_diversity(col_name = "CONTRAST_sex", group_by = 1)
    Condition
      Error in `taxa$alpha_diversity()`:
      ! 1 must be a character and of length 1

---

    Code
      taxa$alpha_diversity(col_name = "CONTRAST_sex", group_by = "1")
    Condition
      Error in `taxa$alpha_diversity()`:
      ! The specified "1" does not exist in the metaData.

---

    Code
      taxa$alpha_diversity(col_name = "CONTRAST_sex", group_by = c("1", "2"))
    Condition
      Error in `taxa$alpha_diversity()`:
      ! "1" and "2" must be a character and of length 1

---

    Code
      taxa$alpha_diversity(col_name = "CONTRAST_sex", evenness = "FALSE")
    Condition
      Error in `taxa$alpha_diversity()`:
      ! "evenness" can only be a `TRUE` or `FALSE`.

---

    Code
      taxa$alpha_diversity(col_name = "CONTRAST_sex", paired = "FALSE")
    Condition
      Error in `taxa$alpha_diversity()`:
      ! "paired" can only be a `TRUE` or `FALSE`.

---

    Code
      taxa$alpha_diversity(col_name = "CONTRAST_sex", p.adjust.method = "nothing")
    Condition
      Error in `taxa$alpha_diversity()`:
      ! Specified "nothing" is not valid.  Valid options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", and "none"

# `omics$alpha_diversity()` -- Behavioral checks

    Code
      adiv$data
    Output
               V1 CONTRAST_sex
            <num>       <char>
      1: 3.403898         male
      2: 3.776849       female
      3: 3.682609       female
      4: 3.686005         male

---

    Code
      taxa
    Message
      
      -- <metagenomics> object 
      metaData: 9 variables x 4 samples
      countData: 4 samples x 242 features
      featureData: 7 attributes x 242 features
      treeData: 242 tips x 241 nodes

---

    Code
      adiv$data
    Output
               V1 CONTRAST_sex group_col
            <num>       <char>    <char>
      1: 3.403898         male     tumor
      2: 3.776849       female     tumor
      3: 3.682609       female   healthy
      4: 3.686005         male   healthy

