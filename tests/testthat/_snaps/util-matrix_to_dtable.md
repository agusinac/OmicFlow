# `matrix_to_dtable()` -- Argument checks

    Code
      matrix_to_dtable(x = list())
    Condition
      Error in `matrix_to_dtable()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      matrix_to_dtable(x = data.frame())
    Condition
      Error in `matrix_to_dtable()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

# `matrix_to_dtable()` -- Behavioral checks

    Code
      dt
    Output
            S100  S103  S115  S120
           <num> <num> <num> <num>
        1:     0     0     4     0
        2:     0     0     0     2
        3:     0     2     0     0
        4:     0     0     0     3
        5:     0     0     0    28
       ---                        
      238:     0     0     9     0
      239:     0     0    14     0
      240:     7     0     0     0
      241:     3     0     0     0
      242:     0     0     0     3

