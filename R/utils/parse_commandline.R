parse_commandline <- function() {
  option_list <- list (make_option(c("-m","--metadata"), 
                                   action = "store",
                                   help="tab seperated file"),
                       make_option(c("-b", "--biom"),
                                   action = "store",
                                   help="biom format file"),
                       make_option(c("-t", "--tree"),
                                   action = "store",
                                   help="Phylogenetic tree in newick format"),
                       make_option(c("-r", "--refseq"),
                                   action = "store",
                                   help="Reference sequence in fasta format"),
                       make_option(c("-c", "--colname"),
                                   action = "store",
                                   help="RANKSTAT_name from metadata to be used for analysis")
  )
  
  parser <- OptionParser(option_list = option_list)
  arguments <- optparse::parse_args(parser, positional_arguments=TRUE)
  
  opt <- arguments$options

  if (base::grepl(".tsv$", opt$metadata)) {
    metadata <- phyloseq::import_qiime_sample_data(mapfilename = opt$metadata)
  } else warning("Please provide a tab-separated metadata file!")
  if (base::grepl(".biom$", opt$biom)) {
    biom_data <- phyloseq::import_biom(BIOMfilename = opt$biom)
  } else warning("Please provide a valid biom file")
  if (base::grepl(".nwk$", opt$tree)) {
    tree <- phyloseq::read_tree(treefile = opt$tree)
  } else tree <- NULL
  if (base::grepl(".fasta$", opt$refseq)) {
    refseq <- Biostrings::readDNAStringSet(filepath = opt$refseq)
  } else refseq <- NULL
  
  if (opt$colname) {
    result <- list(
      metaData = metadata,
      biomData = biom_data,
      treeData = tree,
      refseqData = refseq,
      col_name = opt$colname
    )
  } else warning("colname is not provided, this is required for the data analysis!") |> quit(status = 2)
  
  return(result)
}