# `metagenomics` -- Argument checks

    Code
      metagenomics$new()
    Condition
      Error in `super$initialize()`:
      ! metaData cannot be empty, please provide a <data.frame>, <data.table> or "filepath"

---

    Code
      metagenomics$new(biomData = "nonexisting.biom")
    Condition
      Error in `super$initialize()`:
      ! metaData cannot be empty, please provide a <data.frame>, <data.table> or "filepath"

---

    Code
      metagenomics$new(featureData = features_file)
    Condition
      Error in `super$initialize()`:
      ! metaData cannot be empty, please provide a <data.frame>, <data.table> or "filepath"

---

    Code
      metagenomics$new(countData = counts_sparse_file)
    Condition
      Error in `super$initialize()`:
      ! metaData cannot be empty, please provide a <data.frame>, <data.table> or "filepath"

---

    Code
      metagenomics$new(metaData = data.frame())
    Condition
      Error in `super$initialize()`:
      ! Error in metaData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      metagenomics$new(metaData = data.table::data.table())
    Condition
      Error in `super$initialize()`:
      ! Error in metaData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      metagenomics$new(metaData = metadata_file, featureData = data.frame())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `super$initialize()`:
      ! Error in featureData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      metagenomics$new(metaData = metadata_file, featureData = data.table::data.table())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `super$initialize()`:
      ! Error in featureData:
      x Input must be an existing "filepath", non-empty <data.frame> or <data.table>.

---

    Code
      metagenomics$new(metaData = metadata_file, biomData = "nonexisting.biom")
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! biomData doesn't exist, please provide an existing "filepath"

---

    Code
      metagenomics$new(metaData = metadata_file, biomData = metadata_file)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! biomData could not be loaded. Not a valid JSON or HDF5 format!

---

    Code
      metagenomics$new(metaData = metadata_file, countData = data.frame())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `super$initialize()`:
      ! Error in countData:
      x "filepath", <data.frame> or <data.table> cannot be empty.

---

    Code
      metagenomics$new(metaData = metadata_file, countData = data.table::data.table())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `super$initialize()`:
      ! Error in countData:
      x "filepath", <data.frame> or <data.table> cannot be empty.

---

    Code
      metagenomics$new(metaData = metadata_file, countData = matrix(0))
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
      metagenomics$new(metaData = metadata_file, treeData = data.frame())
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! countData cannot be empty.. did you forgot to specify the countData or biomData in `metagenomics$new()` ?

---

    Code
      metagenomics$new(metaData = metadata_file, treeData = ape::rtree(50))
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
    Condition
      Error in `initialize()`:
      ! countData cannot be empty.. did you forgot to specify the countData or biomData in `metagenomics$new()` ?

---

    Code
      metagenomics$new(metaData = metadata_file, biomData = biom_hdf5, treeData = ape::rtree(
        50))
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.
      v treeData is loaded.
      i Final steps .. cleaning & creating back-up
    Condition
      Error in `private$sync()`:
      ! None FEATURE_IDs are matching, check if "FEATURE_ID" matches the tip labels in treeData!

# `metagenomics` -- Behavioral checks

    Code
      test <- metagenomics$new(metaData = metadata_file, biomData = biom_hdf5)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.
      i Final steps .. cleaning & creating back-up
      
      -- <metagenomics> object 
      metaData: 9 variables x 4 samples
      countData: 4 samples x 242 features
      featureData: 7 attributes x 242 features

---

    Code
      test <- metagenomics$new(biomData = biom_json, metaData = data.table::data.table(
        SAMPLE_ID = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5",
          "Sample6")))
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.
      i Final steps .. cleaning & creating back-up
      
      -- <metagenomics> object 
      metaData: 1 variables x 6 samples
      countData: 6 samples x 5 features
      featureData: 7 attributes x 5 features

---

    Code
      test <- metagenomics$new(metaData = metadata_file, biomData = biom_hdf5,
        treeData = tree)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.
      v treeData is loaded.
      i Final steps .. cleaning & creating back-up
      
      -- <metagenomics> object 
      metaData: 9 variables x 4 samples
      countData: 4 samples x 242 features
      featureData: 7 attributes x 242 features
      treeData: 242 tips x 241 nodes

---

    Code
      taxa_ref <- metagenomics$new(countData = test$countData, metaData = test$
        metaData, treeData = test$treeData, featureData = test$featureData)
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.
      v treeData is loaded.
      i Final steps .. cleaning & creating back-up
      
      -- <metagenomics> object 
      metaData: 9 variables x 4 samples
      countData: 4 samples x 242 features
      featureData: 7 attributes x 242 features
      treeData: 242 tips x 241 nodes

---

    Code
      test <- metagenomics$new(metaData = metadata_file, biomData = biom_hdf5,
        feature_names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
          "Species", "variants"))
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.
      i Final steps .. cleaning & creating back-up
      
      -- <metagenomics> object 
      metaData: 9 variables x 4 samples
      countData: 4 samples x 242 features
      featureData: 7 attributes x 242 features

---

    Code
      test <- metagenomics$new(metaData = metadata_file, biomData = biom_hdf5,
        feature_names = c("Kingdom", "Phylum", "Class"))
    Message
      v metaData template passed the JSON validation.
      i Checking for duplicated identifiers ..
      v featureData is loaded.
      v countData is loaded.
      i Final steps .. cleaning & creating back-up
      ! The number of columns in featureData is greater than "feature_names". 
      Please check and rename the featureData columns by yourself!
      
      -- <metagenomics> object 
      metaData: 9 variables x 4 samples
      countData: 4 samples x 242 features
      featureData: 7 attributes x 242 features

