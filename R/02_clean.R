# Data wrangling --------------------------------------------------------------#
# In paired samples, should contain an unique sample ID column, then str_extract should be removed!
ps_abs <- ps %>% 
  subset_taxa(Domain == "Bacteria") %>% 
  removeZeros()

# normalize by bacterial domain
ps_rel <- ps_abs %>% 
  transform_sample_counts(function(x) x / sum(x))


# Save wrangled file ----------------------------------------------------------#
saveRDS(ps_rel, "02_ps_rel_bac_norm.rds")
saveRDS(ps_abs, "02_ps_abs_bac.rds")
