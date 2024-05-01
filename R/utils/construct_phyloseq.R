# Creates a phyloseq S4 class object from different files

construct_phyloseq <- function(biom_data, meta_data, tree_data, refseq_data) {
  ps <- phyloseq::merge_phyloseq(biom_data, meta_data, read_tree(tree_data), refseq_data)
  
  # Rename taxa ranks
  base::colnames(phyloseq::tax_table(ps)) <- c("Domain", "Phylum",
                               "Class", "Order", "Family", 
                               "Genus", "Species")
  # Rename taxa names
  df <- data.frame(base::lapply(data.frame(tax_table(ps)),
                          function(x) gsub("^[dpcofgs]_{2}", "", x)),
                   stringsAsFactors = FALSE)
  row.names(df) <- row.names(tax_table(ps))
  phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(df))
  
  return(ps)
} 