# TO DO
library("tidyverse")

getwd()
directory = "../../data-mock"
files <- list.files(path = ".", "\\.txt$", full.names = TRUE)

otu_list <- list()
tax_list <- list()

for (file in files) {
  print(file)
  break
}