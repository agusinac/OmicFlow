# `omics` -- Argument checks

    Code
      omics$new()
    Condition
      Error in `initialize()`:
      ! metaData cannot be empty, please provide a <data.frame>, <data.table> or "filepath"

---

    Code
      omics$new(featureData = features_file)
    Condition
      Error in `initialize()`:
      ! metaData cannot be empty, please provide a <data.frame>, <data.table> or "filepath"

---

    Code
      omics$new(countData = counts_sparse_file)
    Condition
      Error in `initialize()`:
      ! metaData cannot be empty, please provide a <data.frame>, <data.table> or "filepath"

---

    Code
      omics$new(metaData = data.frame())
    Condition
      Error in `initialize()`:
      ! Error in metaData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      omics$new(metaData = data.table::data.table())
    Condition
      Error in `initialize()`:
      ! Error in metaData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      omics$new(metaData = metadata_file, featureData = data.frame())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! Error in featureData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      omics$new(metaData = metadata_file, featureData = data.table::data.table())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! Error in featureData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      omics$new(metaData = metadata_file, countData = data.frame())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! Error in countData:
      x Input must be an existing "filepath", non-empty <matrix> or <Matrix>.

---

    Code
      omics$new(metaData = metadata_file, countData = data.table::data.table())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! Error in countData:
      x Input must be an existing "filepath", non-empty <matrix> or <Matrix>.

---

    Code
      omics$new(metaData = metadata_file, countData = matrix(0))
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v countData is loaded.
      ! Created placeholder featureData.
    Condition
      Error in `private$sync()`:
      ! None SAMPLE_IDs are matching, check if "SAMPLE_ID" are matching the colnames in countData!

# `omics` -- Behavioral checks

    Code
      omics$new(metaData = metadata_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      
      -- <omics> object 
      metaData: 9 variables x 4 samples

---

    Code
      omics$new(metaData = metadata_file, featureData = features_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      
      -- <omics> object 
      metaData: 9 variables x 4 samples
      featureData: 7 attributes x 242 features

---

    Code
      test <- omics$new(metaData = metadata_file, countData = counts_sparse_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v countData is loaded.
      ! Created placeholder featureData.

---

    Code
      test <- omics$new(metaData = metadata_file, countData = counts_sparse_with_rownames_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v countData is loaded.
      ! Created placeholder featureData.

---

    Code
      sparse <- omics$new(metaData = metadata_file, featureData = features_file,
        countData = counts_sparse_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.

---

    Code
      dense <- omics$new(metaData = metadata_file, featureData = features_file,
        countData = counts_dense_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.

---

    Code
      ends_with_featureData
    Message
      
      -- <omics> object 
      metaData: 1 variables x 5 samples
      countData: 5 samples x 100 features
      featureData: 0 attributes x 100 features

---

    Code
      ends_with_countData
    Message
      
      -- <omics> object 
      metaData: 1 variables x 5 samples
      countData: 5 samples x 100 features
      featureData: 0 attributes x 100 features

