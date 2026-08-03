# `column_exists()` -- Argument checks

    Code
      column_exists(column = 1)
    Condition
      Error in `column_exists()`:
      ! 1 needs to contain characters.

---

    Code
      column_exists(column = "1", table = matrix())
    Condition
      Error in `column_exists()`:
      ! `table` must be a <data.frame> or <data.table>.

