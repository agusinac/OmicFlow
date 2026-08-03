# `combine_conditions()` -- Argument checks

    Code
      combine_conditions(condition1 = matrix())
    Condition
      Error in `combine_conditions()`:
      ! `condition1` must be a <data.frame> or <data.table>.

---

    Code
      combine_conditions(condition1 = list())
    Condition
      Error in `combine_conditions()`:
      ! `condition1` must be a <data.frame> or <data.table>.

---

    Code
      combine_conditions(condition2 = matrix())
    Condition
      Error in `combine_conditions()`:
      ! argument "condition1" is missing, with no default

---

    Code
      combine_conditions(condition2 = list())
    Condition
      Error in `combine_conditions()`:
      ! argument "condition1" is missing, with no default

---

    Code
      combine_conditions(condition1 = matrix(), condition2 = data.frame())
    Condition
      Error in `combine_conditions()`:
      ! `condition1` must be a <data.frame> or <data.table>.

---

    Code
      combine_conditions(condition1 = data.frame(), condition2 = matrix())
    Condition
      Error in `combine_conditions()`:
      ! `condition2` must be a <data.frame> or <data.table>.

