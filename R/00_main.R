# Load libraries ---------------------------------------------------------------
library("dplyr")
library("ggplot2")
library("tidyr")
library("biomformat")
library("phyloseq")
library("ape")
library("microViz")
library("Biostrings")
library("patchwork")
library("vegan")
library("optparse")

# Load Functions ---------------------------------------------------------------
# Data wrangling functions
Rscript <- sub("--file=", "", commandArgs()[4])
current_path <- sub(basename(Rscript), "", normalizePath(Rscript))

source(file = paste0(current_path, "99_project-functions.R"))

# Utils
sourceDir(path = paste0(current_path, "utils"))

# Parse command line
data_00 <- parse_commandline()

# Run all scripts --------------------------------------------------------------
source(file = paste0(current_path, "01_load.R"))
source(file = paste0(current_path, "02_clean.R"))
source(file = paste0(current_path, "03_data_visualization.R"))
source(file = paste0(current_path, "04_Model-UniFrac-PCoA.R"))
source(file = paste0(current_path, "05_Model-RDA.R"))
rmarkdown::render(input = paste0(current_path, "../documents/report.Rmd"),
                  output_file = paste0(data_00$outDir, "trial_report.html"))
