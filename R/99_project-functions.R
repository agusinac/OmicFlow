# Data wrangling functions -----------------------------------------------------#
format_metadata <- function(metadata_file) {
  if (base::grepl(".csv", metadata_file, fixed = TRUE)){
    tsv_metadata_file <- stringr::str_replace(metadata_file, ".csv", ".tsv")
    utils::write.table(x = readr::read_csv(metadata_file),
                file = tsv_metadata_file,
                row.names = FALSE,
                sep = "\t")
    return(tsv_metadata_file)
  } 
  else
    return(metadata_file)
}

# Remove zero's from a phyloseq object
removeZeros <- function(x) {
  if (is(x, "phyloseq")) {
    x <- phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x)
    x <- phyloseq::prune_samples(phyloseq::sample_sums(x) > 0, x)
    return(x)
  }
}

ASVs_filter <- function(ps, ASVs) {
  # Replace phyloseq objects by matching ASVs
  phyloseq::otu_table(ps) <- phyloseq::otu_table(ps)[ASVs, ]
  phyloseq::tax_table(ps) <- phyloseq::tax_table(ps)[ASVs, ]
  ps <- removeZeros(ps)
  return(ps)
}

# Collects all functions from a directory
sourceDir <- function(path, trace = TRUE, ...) {
  op <- base::options(); on.exit(base::options(op)) # to reset after each 
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
    options(op)
  }
}

# Simple log transformation of otu table
logn <- function(otu_tab, scalar=1) {
  otu_tab.log <- ( scalar * otu_tab ) + 1 
  otu_tab.log <- log( otu_tab.log ) 
  otu_tab.sc <- scale(otu_tab.log, center = TRUE, scale = FALSE)
  return(otu_tab.sc)
}

# Transform any taxa rank into a otu_table;
# with rows as samples and cols as taxa (rank to be specified)
get_otu <- function(ps, tax_target="Genus") {
  ps_ref <- ps %>% 
    phyloseq::tax_glom(taxrank = tax_target)
  df <- ps_ref %>% 
    phyloseq::otu_table() %>% 
    t() %>% 
    as.data.frame()
  colnames(df) <- as.data.frame(phyloseq::tax_table(ps_ref))[[ tax_target ]]
  return(df)
}

get_meta <- function(ps) {
  # Fetches meta table in right dataframe format
  meta_tab <- ps %>% 
    phyloseq::sample_data() %>% unclass() %>% as.data.frame()
  
  return(meta_tab)
}

top_taxa <- function(asv.tab, n) {
  # Top Taxa req
  taxaSum_df <- rowSums(asv.tab) %>% 
    as.data.frame() %>% 
    arrange(desc(.)) %>% 
    head(n)
  filt_taxaSums <- asv.tab[rownames(asv.tab) %in% rownames(taxaSum_df),]
  return(filt_taxaSums)
}

phyloseq_from_txt <- function(otu_file, tax_file) {
  otu_tab <- read.delim(otu_file)
  tax_tab <- read.delim(tax_file)
  
  tax_tab <- separate(data = tax_tab, 
                      col = X.OTU.ID, sep = ";", 
                      into = c("Domain", "Phylum", "Class",
                               "Order", "Family", "Genus", "Species"))
  
  OTU <- otu_tab %>% as.matrix %>% phyloseq::otu_table(taxa_are_rows = TRUE)
  TAX <- tax_tab %>% as.matrix %>% phyloseq::tax_table()
  
  ps <- phyloseq::phyloseq(OTU, TAX)
  return(ps)
}

add_metadata <- function(df, meta_tab, meta_col.id, meta_col.add) {
  for (variable in meta_col.add) {
    df <- df %>% 
      rowwise() %>% 
      mutate(
        !!variable := meta_tab[[ {{ variable }} ]][str_detect(as.character(meta_tab[[ {{ meta_col.id }} ]]), as.character(variable))]
      )
  }
  return(df)
}