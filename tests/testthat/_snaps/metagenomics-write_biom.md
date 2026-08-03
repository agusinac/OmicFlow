# `write_biom()` -- Argument check

    Code
      taxa$write_biom(filename = list())
    Condition
      Error in `taxa$write_biom()`:
      ! Needs to contain characters and be length of 1.

---

    Code
      taxa$write_biom(filename = c("file1.biom", "file2.biom"))
    Condition
      Error in `taxa$write_biom()`:
      ! file1.biom and file2.biom Needs to contain characters and be length of 1.

