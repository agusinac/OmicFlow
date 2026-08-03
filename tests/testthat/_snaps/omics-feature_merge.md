# `omics$feature_merge()` -- Argument checks

    Code
      taxa$feature_merge(feature_rank = 1)
    Condition
      Error in `taxa$feature_merge()`:
      ! 1 needs to be a character or vector containing characters

---

    Code
      taxa$feature_merge(feature_rank = list())
    Condition
      Error in `taxa$feature_merge()`:
      ! needs to be a character or vector containing characters

---

    Code
      taxa$feature_merge(feature_rank = "1")
    Condition
      Error in `taxa$feature_merge()`:
      ! "1" does not exist in featureData!

---

    Code
      taxa$feature_merge(feature_rank = c("1", "Kingdom"))
    Condition
      Error in `taxa$feature_merge()`:
      ! "1" and "Kingdom" does not exist in featureData!

---

    Code
      taxa$feature_merge(feature_rank = "Kingdom", feature_filter = list())
    Condition
      Error in `taxa$feature_merge()`:
      ! needs to be a character or vector containing characters

---

    Code
      taxa$feature_merge(feature_rank = "Kingdom", feature_filter = 1)
    Condition
      Error in `taxa$feature_merge()`:
      ! 1 needs to be a character or vector containing characters

# `omics$feature_merge()` -- Behavioral checks

    Code
      taxa
    Message
      
      -- <metagenomics> object 
      metaData: 9 variables x 4 samples
      countData: 4 samples x 63 features
      featureData: 7 attributes x 63 features
      treeData: 63 tips x 62 nodes

---

    Code
      taxa
    Message
      
      -- <metagenomics> object 
      metaData: 9 variables x 4 samples
      countData: 4 samples x 76 features
      featureData: 7 attributes x 76 features
      treeData: 76 tips x 75 nodes

