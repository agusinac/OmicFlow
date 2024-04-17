# Data load -------------------------------------------------------------------#
radboud_ps_abs <- readRDS("../../Pathology/projects/FFPE_breast_microbiome/data/RDS/01_absolute_phyloseq.rds")

# Data wrangling --------------------------------------------------------------#
sample_names(radboud_ps_abs) <- sample_names(radboud_ps_abs) %>%
  str_extract("(^[^_]*_[^_]*_?\\d)|(C\\d?)") 
radboud_ps_abs_bac <- radboud_ps_abs %>% 
  subset_taxa(Domain == "Bacteria") 

# normalize by bacterial domain
radboud_ps_rel_bac_norm <- radboud_ps_abs_bac %>% 
  transform_sample_counts(function(x) x / sum(x))


# Save wrangled file ----------------------------------------------------------#
saveRDS(radboud_ps_rel_bac_norm, "../../Pathology/projects/FFPE_breast_microbiome/data/RDS/02_ps_rel_bac_norm.rds")
saveRDS(radboud_ps_abs_bac, "../../Pathology/projects/FFPE_breast_microbiome/data/RDS/01_ps_abs_bac.rds")
