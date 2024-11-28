# Load Library -----------------------------------------------------------------
devtools::load_all()
library("Matrix")
library("ggplot2")
library("ggtree")
library("patchwork")

# Fetch current Rscript path
Rscript <- sub("--file=", "", commandArgs()[4])
current_path <- sub(basename(Rscript), "", normalizePath(Rscript))

# Parse command line
options <- parse_commandline()

outfile_path <- normalizePath(options$outdir)

# main -------------------------------------------------------------------------
# switch statement based on omic selected, create object
omics <- metataxonomics$new(metaData = base::file.path(current_path, options$metadata),
                            biomData = base::file.path(current_path, options$biom),
                            treeData = base::file.path(current_path, options$tree))

# Set parameters
feature_ranks = c("Phylum", "Family", "Genus")
distance_metrics = c("unifrac", "bray")
feature_filter = c("uncultured")
rankstat_labels <- sub("RANKSTAT_", "", colnames(omics$metaData)[grepl("RANKSTAT_", colnames(omics$metaData))])

# Initiate OmicFlow
plots <- omics$autoFlow(feature_ranks = feature_ranks,
                        feature_filter = feature_filter,
                        distance_metrics = distance_metrics,
                        cpus = options$cpus)

# Create report
rmarkdown::render(input = paste0(current_path, "report/report.Rmd"),
                  output_file = paste0(outfile_path, "/report.html"))
