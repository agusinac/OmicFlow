# OmicFlow

[![CRAN](https://www.r-pkg.org/badges/version/OmicFlow)](https://CRAN.R-project.org/package=OmicFlow)
[![Codecov](https://codecov.io/gh/agusinac/OmicFlow/graph/badge.svg)](https://app.codecov.io/gh/agusinac/OmicFlow)
[![R-CMD-check](https://github.com/agusinac/OmicFlow/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/agusinac/OmicFlow/actions/workflows/R-CMD-check.yaml)
[![run with
conda](https://anaconda.org/agusinac/r-omicflow/badges/version.svg)](https://anaconda.org/agusinac/r-omicflow)
[![run with
docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://hub.docker.com/r/agusinac/autoflow)

## Overview

OmicFlow is a generalised data structure for fast and efficient loading
of various sparse omics data. It can handle metataxonomics/metagenomics
data in text or
[BIOM](https://biom-format.org/documentation/format_versions/biom-2.0.html)
and extends to `proteomics` and other `omics` types. It also supports
non-sparse data, but it’s performance peaks in sparsity.

## Installation

The latest stable version can be installed from CRAN.

``` r

install.packages('OmicFlow', dependencies = TRUE)
```

The development version is available on GitHub.

``` r

install.packages('pak') # if not yet installed
pak::pkg_install('agusinac/OmicFlow@dev')
```

## Usage

Initialize the `metagenomics` or any `omics` object from a file-path or
pre-loaded object.

``` r
library("OmicFlow")

metadata_file <- system.file("extdata", "metadata.tsv", package = "OmicFlow")
counts_file <- system.file("extdata", "counts.tsv", package = "OmicFlow")
features_file <- system.file("extdata", "features.tsv", package = "OmicFlow")
tree_file <- system.file("extdata", "tree.newick", package = "OmicFlow")

taxa <- metagenomics$new(
    metaData = metadata_file,
    countData = counts_file,
    featureData = features_file,
    treeData = tree_file
)

taxa$feature_subset(Kingdom == "Bacteria")
taxa$scale(method = "tss")

# Access variables directly
taxa$metaData
taxa$countData
taxa$featureData
taxa$treeData

# Change variables & enjoy the automatic sync
taxa$featureData <- taxa$featureData[1:100, ]

# Inspect what functions variables are available to the class
taxa$
```

If you are new to `OmicFlow`, the best place to start is the
[Introduction to
OmicFlow](https://agusinac.github.io/OmicFlow/articles/getting-started.html).

## Docker

Example: **Outputs a `report.html` file in current work directory**

``` bash
docker pull agusinac/autoflow:latest

docker run -it --rm -v \
    "$(pwd)":/data \             # Mount the data in a temporary directory
    -w /data \                   # set working directory
    -u $(id -u):$(id -g) \       # non-root user
    agusinac/autoflow:latest \
    autoflow \                   # autoflow R script
    -b /data/biom_with_taxonomy_hdf5.biom \
    -m /data/metadata.tsv
```

## Support

If you are having issues, please [create a
ticket](https://github.com/agusinac/OmicFlow/issues)

## Acknowledgements

- Sticker logo created by [studio
  floorescent](https://www.floorescent.nl/) \| [Floor
  Baas](https://nl.linkedin.com/in/floorbaas)
