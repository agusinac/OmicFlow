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
    Condition
      Error:
      ! object 'metadata_file' not found

---

    Code
      metagenomics$new(metaData = metadata_file, featureData = data.table::data.table())
    Condition
      Error:
      ! object 'metadata_file' not found

---

    Code
      metagenomics$new(metaData = metadata_file, biomData = "nonexisting.biom")
    Condition
      Error:
      ! object 'metadata_file' not found

---

    Code
      metagenomics$new(metaData = metadata_file, biomData = metadata_file)
    Condition
      Error:
      ! object 'metadata_file' not found

---

    Code
      metagenomics$new(metaData = metadata_file, countData = data.frame())
    Condition
      Error:
      ! object 'metadata_file' not found

---

    Code
      metagenomics$new(metaData = metadata_file, countData = data.table::data.table())
    Condition
      Error:
      ! object 'metadata_file' not found

---

    Code
      metagenomics$new(metaData = metadata_file, countData = matrix(0))
    Condition
      Error:
      ! object 'metadata_file' not found

---

    Code
      metagenomics$new(metaData = metadata_file, treeData = data.frame())
    Condition
      Error:
      ! object 'metadata_file' not found

---

    Code
      metagenomics$new(metaData = metadata_file, treeData = ape::rtree(50))
    Condition
      Error:
      ! object 'metadata_file' not found

---

    Code
      metagenomics$new(metaData = metadata_file, biomData = biom_hdf5, treeData = ape::rtree(
        50))
    Condition
      Error:
      ! object 'metadata_file' not found

