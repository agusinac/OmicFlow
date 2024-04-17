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
                                   help="Reference sequence in fasta format")
  )
  
  parser <- OptionParser(option_list = option_list)
  arguments <- optparse::parse_args(parser, positional_arguments=TRUE)
  
  opt <- arguments$options
  print(opt$metadata)
  if (base::grepl(".tsv$", opt$metadata)) {
    metadata <- phyloseq::import_qiime_sample_data(mapfilename = opt$metadata)
  }
  if (base::grepl(".biom$", opt$biom)) {
    biom_data <- phyloseq::import_biom(BIOMfilename = opt$biom)
  }
  if (base::grepl(".nwk$", opt$tree)) {
    tree <- phyloseq::read_tree(treefile = opt$tree)
  }
  if (base::grepl(".fasta$", opt$refseq)) {
    refseq <- Biostrings::readDNAStringSet(filepath = opt$refseq)
  }
  
  result <- list(
    metaData = metadata,
    biomData = biom_data,
    treeData = tree,
    refseqData = refseq
  )
  
  return(result)
}