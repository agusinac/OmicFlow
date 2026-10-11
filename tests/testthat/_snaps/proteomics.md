# `proteomics` -- Argument checks

    Code
      proteomics$new()
    Condition
      Error in `super$initialize()`:
      ! metaData cannot be empty, please provide a <data.frame>, <data.table> or "filepath"

---

    Code
      proteomics$new(countData = counts_with_rownames_file)
    Condition
      Error in `super$initialize()`:
      ! metaData cannot be empty, please provide a <data.frame>, <data.table> or "filepath"

---

    Code
      proteomics$new(metaData = data.frame())
    Condition
      Error in `super$initialize()`:
      ! Error in metaData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      proteomics$new(metaData = data.table::data.table())
    Condition
      Error in `super$initialize()`:
      ! Error in metaData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      proteomics$new(metaData = metadata_file, featureData = data.frame())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `super$initialize()`:
      ! Error in featureData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      proteomics$new(metaData = metadata_file, featureData = data.table::data.table())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `super$initialize()`:
      ! Error in featureData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      proteomics$new(metaData = metadata_file, countData = data.frame())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `super$initialize()`:
      ! Error in countData:
      x "filepath", <data.frame> or <data.table> cannot be empty.

---

    Code
      proteomics$new(metaData = metadata_file, countData = data.table::data.table())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `super$initialize()`:
      ! Error in countData:
      x "filepath", <data.frame> or <data.table> cannot be empty.

---

    Code
      proteomics$new(metaData = metadata_file, countData = matrix(0))
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v countData is loaded.
      i Created placeholder featureData.
    Condition
      Error in `private$sync()`:
      ! None SAMPLE_IDs are matching, check if "SAMPLE_ID" are matching the colnames in countData!

---

    Code
      proteomics$new(metaData = metadata_file, treeData = data.frame())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! countData cannot be empty.. did you forgot to specify the countData in `proteomics$new()` ?

---

    Code
      proteomics$new(metaData = metadata_file, treeData = ape::rtree(50))
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! countData cannot be empty.. did you forgot to specify the countData in `proteomics$new()` ?

# `proteomics` -- Behavioral checks

    Code
      test <- proteomics$new(metaData = metadata_file, countData = counts_with_rownames_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v countData is loaded.
      i Created placeholder featureData.
      
      -- <proteomics> object 
      metaData: 3 variables x 10 samples
      countData: 10 samples x 50 features
      featureData: 0 attributes x 50 features

---

    Code
      test <- proteomics$new(metaData = metadata_file, countData = counts_without_rownames_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v countData is loaded.
      i Created placeholder featureData.
      
      -- <proteomics> object 
      metaData: 3 variables x 10 samples
      countData: 10 samples x 50 features
      featureData: 0 attributes x 50 features

---

    Code
      test <- proteomics$new(metaData = metadata_file, countData = counts_with_rownames_file,
        treeData = tree_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v countData is loaded.
      i Created placeholder featureData.
      v treeData is loaded.
      
      -- <proteomics> object 
      metaData: 3 variables x 10 samples
      countData: 10 samples x 50 features
      featureData: 0 attributes x 50 features
      treeData: 50 tips x 49 nodes

---

    Code
      prot_ref <- proteomics$new(countData = test$countData, metaData = test$metaData,
      treeData = test$treeData, featureData = test$featureData)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.
      v treeData is loaded.
      
      -- <proteomics> object 
      metaData: 3 variables x 10 samples
      countData: 10 samples x 50 features
      featureData: 0 attributes x 50 features
      treeData: 50 tips x 49 nodes

