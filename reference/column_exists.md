# Checks if given column name exists in provide table

Mainly used within
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) and
other functions to check if given column name(s) exist in the table and
is not completely empty (containing NAs).

## Usage

``` r
column_exists(column, table)
```

## Arguments

- column:

  A column name in `table`.

- table:

  A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html).

## Value

A boolean value.
