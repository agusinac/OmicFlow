# Data load -------------------------------------------------------------------#
radboud_ps_abs <- readRDS("data/RDS/01_absolute_phyloseq.rds")

# Data wrangling --------------------------------------------------------------#
# In paired samples, should contain an unique sample ID column, then str_extract should be removed!
radboud_ps_abs_bac <- radboud_ps_abs %>% 
  subset_taxa(Domain == "Bacteria") %>% 
  removeZeros()

# normalize by bacterial domain
radboud_ps_rel_bac_norm <- radboud_ps_abs_bac %>% 
  transform_sample_counts(function(x) x / sum(x))


# Save wrangled file ----------------------------------------------------------#
saveRDS(radboud_ps_rel_bac_norm, "data/RDS/02_ps_rel_bac_norm.rds")
saveRDS(radboud_ps_abs_bac, "data/RDS/01_ps_abs_bac.rds")
