# `colormap()` -- Argument checks

    Code
      colormap(data = as.matrix(taxa$metaData))
    Condition
      Error in `colormap()`:
      ! Data must be a <data.frame> or <data.table>.

---

    Code
      colormap(data = data.frame(), col_name = "nonexisting")
    Condition
      Error in `colormap()`:
      ! The "nonexisting" column does not exist in the provided data.

---

    Code
      colormap(data = data.frame(), col_name = c("col1", "col2"))
    Condition
      Error in `colormap()`:
      ! "col1" and "col2" needs to contain characters with length of 1.

---

    Code
      colormap(data = taxa$metaData, col_name = "CONTRAST_sex", Brewer.palID = 2)
    Condition
      Error in `colormap()`:
      ! The 2 needs to contain characters with length of 1.

---

    Code
      colormap(data = taxa$metaData, col_name = "CONTRAST_sex", Brewer.palID = "colSet")
    Condition
      Error in `colormap()`:
      ! "colSet" is not a valid Brewer pal ID.  Valid options: "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3", "Blues", ..., "YlOrBr", and "YlOrRd".

