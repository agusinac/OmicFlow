# tmp run command:
# Rscript 00_main.R --metadata ../../Pathology/projects/FFPE_breast_microbiome/data/metadata/metadata.tsv --biom ../../Pathology/projects/FFPE_breast_microbiome/data/BiotaViz/relative-table-with-taxonomy.biom --tree ../../Pathology/projects/FFPE_breast_microbiome/data/phylogeny/phylogenetic_tree/tree.nwk --refseq ../../Pathology/projects/FFPE_breast_microbiome/data/denoise/representative_sequences/unfiltered/sequences.fasta

# Load libraries ---------------------------------------------------------------
library("tidyverse")
library("biomformat")
library("phyloseq")
library("ape")
library("microViz")
library("decontam")
library("Biostrings")
library("patchwork")
library("vegan")
library("ComplexHeatmap")
library("optparse")
library("rmarkdown")

# Load Functions ---------------------------------------------------------------
# Data wrangling functions
source(file = "99_project-functions.R")

# Utils
sourceDir(path = "utils")

# Parse command line
data_00 <- parse_commandline()

# Run all scripts --------------------------------------------------------------
source(file = "01_load.R")
source(file = "02_clean.R")
source(file = "03_data_visualization.R")
source(file = "04_Model-UniFrac-PCoA.R")
source(file = "05_Model-RDA.R")
rmarkdown::render(input = "Report.Rmd",
                  output_file = "trial_report.html",
                  output_dir = ".")
