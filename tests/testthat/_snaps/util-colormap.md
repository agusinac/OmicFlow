# `colormap()` -- Argument checks

    Code
      colormap(data = as.matrix(taxa$metaData))
    Condition
      Error in `colormap()`:
      ! "data" must be a <data.frame> or <data.table>.

---

    Code
      colormap(data = data.frame(), groups = "nonexisting")
    Condition
      Error in `colormap()`:
      ! The "groups" column does not exist in the provided data.

---

    Code
      colormap(data = data.frame(), groups = c("col1", "col2"))
    Condition
      Error in `colormap()`:
      ! "groups" needs to contain characters with length of 1.

---

    Code
      colormap(data = taxa$metaData, groups = "CONTRAST_sex", Brewer.palID = 2)
    Condition
      Error in `colormap()`:
      ! The "Brewer.palID" needs to contain characters with length of 1.

---

    Code
      colormap(data = taxa$metaData, groups = "CONTRAST_sex", Brewer.palID = "colSet")
    Condition
      Error in `colormap()`:
      ! "Brewer.palID" is not a valid Brewer pal ID.  Valid options: "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3", "Blues", ..., "YlOrBr", and "YlOrRd".

