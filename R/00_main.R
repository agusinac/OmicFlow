# Load libraries ---------------------------------------------------------------
library("tidyverse")
library("biomformat")
library("phyloseq")
library("ape")
library("microViz")
library("Biostrings")
library("patchwork")
library("vegan")
library("optparse")
library("rmarkdown")
library("ComplexHeatmap")

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
rmarkdown::render(input = "../documents/Report.Rmd",
                  output_file = "trial_report.html",
                  output_dir = "../documents")
