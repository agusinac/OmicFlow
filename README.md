[![CRAN](https://www.r-pkg.org/badges/version/OmicFlow)](https://CRAN.R-project.org/package=OmicFlow)
[![Codecov](https://codecov.io/gh/agusinac/OmicFlow/graph/badge.svg)](https://app.codecov.io/gh/agusinac/OmicFlow)
[![R-CMD-check](https://github.com/agusinac/OmicFlow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/agusinac/OmicFlow/actions/workflows/R-CMD-check.yaml)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://hub.docker.com/r/agusinac/autoflow)

OmicFlow
================

## Installation
---

The latest stable version can be installed from CRAN.

``` r
install.packages('OmicFlow', dependencies = TRUE)
```

The development version is available on GitHub.

``` r
install.packages('pak') # if not yet installed
pak::pak('agusinac/OmicFlow')
```

## 📋 Metadata File Specification

OmicFlow expects your sample metadata to follow a **simple, but strict** structure so that all datasets are compatible and validated up‑front. Sample metadata can be supplied as a **CSV/TSV** file or as a `data.table` in R. In both cases the sample metadata should contain a header (this is your first line if you supply a file) where **each row = one sample** Additional column names not mentioned here are allowed and will be ignored during the metadata validation step.

---

### **Minimum requirement**
- **`SAMPLE_ID`** ➡ every row **must** have a unique, non‑empty sample identifier.
- No spaces are allowed in IDs — use underscores `_` or dashes `-` instead.

Example:

| SAMPLE_ID | SAMPLEPAIR_ID | CONTRAST_Treatment | VARIABLE_Age |
|-----------|---------------|--------------------|--------------|
| S1        | P1            | Drug               | 42           |
| S2        | P1            | Placebo            | 36           |
| S3        | P2            | Drug               | 51           |

---

### **Column types and naming rules**

#### 🔹 Required column
| Column       | Type    | Rules                                                |
|--------------|---------|------------------------------------------------------|
| `SAMPLE_ID`  | string  | Unique, no spaces, one per sample row                |

#### 🔹 Optional standard columns
| Column         | Type    | Rules                                                               |
|----------------|---------|---------------------------------------------------------------------|
| `FEATURE_ID`   | string  | Optional — no spaces. Naming of the feature identifiers to include or exclude certain features  |
| `SAMPLEPAIR_ID`| string  | Optional — no spaces. Use when samples are paired and belong to an individual source/subject |

#### 🔹 Pattern‑based columns
You can define extra variables using special prefixes:
- **`CONTRAST_...`** → grouping/category labels used in differential comparisons  
  Example: `CONTRAST_Treatment` with values `Drug` / `Placebo`
- **`VARIABLE_...`** → numeric or string variables for statistical analysis  
  Example: `VARIABLE_Age` with values `42`, `51`, etc.

The pattern-based columns are only used during the `autoFlow` function. At the moment only columns with prefix `CONTRAST_` are supported.
Example: **Outputs a `report.html` file in the current working directory**
```R
taxa$autoFlow(
    normalize = FALSE,
    weighted = TRUE,
    pvalue.threshold = 0.05
)
```
---

## Usage
> [!NOTE]
> Make sure your metadata meets the requirements!
---
```R
library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")

taxa <- metagenomics$new(
    metaData = metadata_file,
    countData = counts_file,
    featureData = features_file
)

taxa$feature_subset(Kingdom == "Bacteria")
taxa$normalize()

```

### Visualisations
---
> [!NOTE]
> All visualizations use by default color-blind palettes!

#### 🔹Alpha diversity
```R
alpha_div <- taxa$alpha_diversity(
    col_name = "treatment",
    metric = "shannon"
)

alpha_div$plot
```
![](docs/figures/alphadiv_readme.png)

#### 🔹Beta diversity
```R
beta_div <- taxa$ordination(
    metric = "unifrac",
    method = "pcoa",
    group_by = "treatment"
)

patchwork::wrap_plots(
    beta_div[c("scree_plot", "anova_plot", "scores_plot")],
    nrow = 1)
```
![](docs/figures/betadiv_readme.png)

#### 🔹Composition
```R
res <- taxa$composition(
    feature_rank = "Genus",
    feature_filter = c("uncultured"),
    feature_top = 15,
    normalize = FALSE,
    col_name = "CONTRAST_sex"
)

composition_plot(
    data = res$data,
    palette = res$palette,
    feature_rank = "Genus",
    group_by = "CONTRAST_sex"
    )
```
![](docs/figures/composition_readme.png)

#### 🔹Volcano plot
```R
res <- taxa$DFE(
    feature_rank = "Genus",
    feature_filter = c("uncultured"),
    paired = FALSE,
    normalize = FALSE,
    condition.group = "CONTRAST_sex",
    condition_A = "male",
    condition_B = "female"
)

res$volcano_plot
```

### Run autoflow via docker 🐳
For additional options please run `autoflow -h`
```bash
docker pull agusinac/autoflow:1.3.0

docker run -it --rm -v \
    "$(pwd)":/data \             # Mount the data in a temporary directory
    -w /data \                   # set working directory
    -u $(id -u):$(id -g) \       # non-root user
    agusinac/autoflow:1.3.0 \
    autoflow \                   # autoflow R script
    -b /data/biom_with_taxonomy_hdf5.biom \
    -m /data/metadata.tsv
```

## Support
If you are having issues, please [create a ticket](https://github.com/agusinac/OmicFlow/issues)
