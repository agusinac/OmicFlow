# OmicFlow
Inspired by Phyloseq only made to handle large datasets, uses triplet format data structures such as data.table and TsparseMatrix. The OmicFlow started as an automation analysis project to fetch data from a nextflow pipeline and then perform automated analysis. But with phyloseq it was getting very slow and inefficient. The main focus of OmicFlow is to create an efficient way to load and manipulate omics-related data as an object with little memory, therefore the R6 and data.table packages return the object by reference. This approach may not be familiar to most R-users, but it is more efficient when you have large data-sets and you don't want to make many copies of the same.

## Installation [unavailable]
``install.packages("OmicFlow")``

## Usage
Unlike most R packages, the OmicFlow encapsulates a lot of the methods and keeps them in the abstract class; 'tools'. This approach makes it possible to create a more generalized environment for the manipulation and loading of various omics-related data and keeps the naming of functions identical among omics classes. It is easily explandable to add new omics classes, while inheriting from the abstract class.
```r
library("Matrix")
library("ggplot2")
library("OmicFlow")

# Creates a new object
taxa <- metataxonomics$new(metaData = "./tests/testthat/metataxonomics/input/metadata.tsv",
                           biomData = "./tests/testthat/metataxonomics/input/biom_with_taxonomy.biom",
                           treeData = "./tests/testthat/metataxonomics/input/rooted_tree.newick")

# Data files can be easily accessed via '$';
taxa$countData
taxa$treeData

# Perform in place modifications
taxa$feature_subset(Domain == "Bacteria")

# trasnform count data
taxa$transform(function(x) x / sum(x))
taxa$transform(log2)

# Reset all previous changes and return to original loaded data
taxa$reset()

```

## Automatic data analysis
Finally, the automation of data analysis happens via a metadata template, that can be filled in by the user and handled by OmicFlow. Here we need the report folder and autoFlow.R script. This can be easily called from the command-line or wrapped as a autoFlow.nf (nextflow) script.
```bash
Rscript autoFlow.R \
  --metadata "metadata.tsv" \
  --biom "table.biom" \
  --tree "rooted_tree.newick" \
  --cpus 4
```

## Running with Docker
```bash
docker pull agusinac/autoflow:0.0.1

docker run --rm -v $(pwd):/scripts agusinac/autoflow:0.0.1 Rscript scripts/autoFlow.R \
  --metadata "metadata.tsv" \
  --biom "table.biom" \
  --tree "rooted_tree.newick" \
  --cpus 4
```

## metadata template example

