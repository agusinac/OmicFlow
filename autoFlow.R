# Load Library -----------------------------------------------------------------
library("OmicFlow")

# Fetch current Rscript path
Rscript <- sub("--file=", "", commandArgs()[4])
current_path <- sub(basename(Rscript), "", normalizePath(Rscript))

# Parse command line
options <- parse_commandline()
outfile_path <- normalizePath(options$outDir)

# main -------------------------------------------------------------------------
# Run autoFlow
# switch statement based on omic selected, create object
omics <- metataxonomics$new()

# Initiate OmixFlow
omic$autoFlow()
