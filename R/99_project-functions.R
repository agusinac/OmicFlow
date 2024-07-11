# Remove zero's from a phyloseq object
removeZeros <- function(x) {
  if (is(x, "phyloseq")) {
    x <- phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x)
    x <- phyloseq::prune_samples(phyloseq::sample_sums(x) > 0, x)
    return(x)
  }
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

top_taxa <- function(asv.tab, n) {
  # Top Taxa req
  top_taxa_names <- rowSums(asv.tab) %>% 
    sort(decreasing=TRUE) %>% 
    head(n)
  filt_taxaSums <- asv.tab[rownames(asv.tab) %in% names(top_taxa_names),]
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

add_metadata <- function(df_long, meta_tab, meta_col.id, meta_col.add) {
  for (variable in meta_col.add) {
    df_long <- df_long %>% 
      rowwise() %>% 
      mutate(
        !!variable := meta_tab[[ {{ variable }} ]][stringr::str_detect(as.character(meta_tab[[ {{ meta_col.id }} ]]), as.character(`sample-id`))]
      )
  }
  return(df_long)
}

removeDepends <- function(pkg, recursive = FALSE){
  d <- package_dependencies(,installed.packages(), recursive = recursive)
  depends <- if(!is.null(d[[pkg]])) d[[pkg]] else character()
  needed <- unique(unlist(d[!names(d) %in% c(pkg,depends)]))
  toRemove <- depends[!depends %in% needed]
  if(length(toRemove)){
    toRemove <- select.list(c(pkg,sort(toRemove)), multiple = TRUE,
                            title = "Select packages to remove")
    remove.packages(toRemove)
    return(toRemove)
  } else {
    invisible(character())
  }
}
