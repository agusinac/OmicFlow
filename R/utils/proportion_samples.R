# Creates a barplot of relative abundances for each sample. 
# Displays the first 15 top relative abundant taxa at Genus level.

proportion_samples <- function(ps.norm, taxa_rank="Genus", sample_format=FALSE) {
  if (sample_format == TRUE) {
    # Renaming sample names
    phyloseq::sample_names(ps.norm) <- phyloseq::sample_names(ps.norm) %>% 
      stringr::str_extract("(^[^_]*_[^_]*_?\\d)|(C\\d?)")
  }
  per_group_composition <- ps.norm %>% microViz::tax_fix(anon_unique = FALSE,
                                               verbose = FALSE,
                                               unknowns = c("uncultured")) %>% 
    microViz::comp_barplot(tax_level = {{taxa_rank}}, 
                 n_taxa = 15, 
                 merge_other = FALSE, 
                 sample_order = "asis") +
    coord_flip() +
    guides(fill = guide_legend(ncol = 2, title.position = "top")) +
    theme(
      text=element_text(size=20),
      legend.position = "bottom", legend.key.size = unit(5, "mm")
    )
  return(per_group_composition)
}