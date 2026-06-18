# Metadata File Specification

Every class in OmicFlow starts with metadata validation, which checks
that sample IDs match those in the abundance tables.  
The validation is defined via a JSON format in the abstract class
`omics`, it takes as input a **CSV/TSV** file or a `data.table`.  
In both cases, the sample metadata must contain a header (first line if
you supply a file) where **each row = one sample**.  
Additional column names not mentioned here are allowed and will be
ignored during metadata validation.

Below are the exact specifications required for the metadata structure.

## **Minimum requirement**

- **`SAMPLE_ID`** ➡ every row **must** have a unique, non‑empty sample
  identifier.
- No spaces are allowed in IDs — use underscores `_` or dashes `-`
  instead.

Example:

| SAMPLE_ID | SAMPLEPAIR_ID | CONTRAST_Treatment | VARIABLE_Age |
|-----------|---------------|--------------------|--------------|
| S1        | P1            | Drug               | 42           |
| S2        | P1            | Placebo            | 36           |
| S3        | P2            | Drug               | 51           |

------------------------------------------------------------------------

## **Column types and naming rules**

### 🔹 Required column

| Column      | Type   | Rules                                 |
|-------------|--------|---------------------------------------|
| `SAMPLE_ID` | string | Unique, no spaces, one per sample row |

### 🔹 Optional standard columns

| Column | Type | Rules |
|----|----|----|
| `SAMPLEPAIR_ID` | string | Optional — no spaces. Use when samples are paired and belong to an individual source/subject |

### 🔹 Pattern‑based columns

You can define extra variables using special prefixes: -
**`CONTRAST_...`** → grouping/category labels used in differential
comparisons  
Example: `CONTRAST_Treatment` with values `Drug` / `Placebo` -
**`VARIABLE_...`** → numeric or string variables for statistical
analysis  
Example: `VARIABLE_Age` with values `42`, `51`, etc.

The pattern-based columns are only used during the `autoFlow` function.
At the moment only columns with prefix `CONTRAST_` are supported.

## Examples

You can always check metadata up‑front before loading any other data
files by creating an `omics` object with only the `metaData` argument:

``` r

library(OmicFlow)
#> Loading required package: R6
#> Loading required package: data.table
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
#> Loading required package: Matrix

# Check from filepath
metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
test <- omics$new(metaData = metadata_file)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..

# Check from `data.table` object
metadata <- data.table::data.table("SAMPLE_ID" = paste0("Sample_", 1:5))
test <- omics$new(metaData = metadata)
#> ✔ metaData template passed the JSON validation.
#> ℹ Checking for duplicated identifiers ..
```
