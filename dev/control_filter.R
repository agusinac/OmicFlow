ps_filter <- function(ps) {
  ps.control <- ps %>% 
    subset_samples(RANKSTAT_treatment == "control") %>% 
    removeZeros()
  
  ps.positive <- ps %>% 
    subset_samples(RANKSTAT_treatment != "control") %>% 
    removeZeros()
  
  ASVs_overlap <- base::intersect(taxa_names(ps.control), taxa_names(ps.positive))
  
  # Create subset phyloseq
  original_ASVs <- taxa_names(ps)
  Keep_ASVs <- base::setdiff(original_ASVs, ASVs_overlap)
  
  ps.keep <- phyloseq::prune_taxa(Keep_ASVs, ps) %>% 
    removeZeros()
  
  return(ps.keep)
} 
