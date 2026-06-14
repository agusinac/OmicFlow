# Checks if column exists in table

Mainly used within
[omics](https://agusinac.github.io/OmicFlow/reference/omics.md) and
other functions to check if given column name does exist in the table
and is not completely empty (containing NAs).

## Usage

``` r
column_exists(column, table)
```

## Arguments

- column:

  A character of length 1.

- table:

  A [data.table](https://rdrr.io/pkg/data.table/man/data.table.html) or
  [data.frame](https://rdrr.io/r/base/data.frame.html).

## Value

A boolean value.
