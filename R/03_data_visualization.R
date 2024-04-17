# Load phyloseq objects -------------------------------------------------------#
ps_abs <- readRDS("../../Pathology/projects/FFPE_breast_microbiome/data/RDS/01_ps_abs_bac.rds")


## Rankstat
### Classified ASVs per taxonomic rank
#------------------------------------------------------------------------------#
proportion_classified <- plot_proportion_classified(ps_abs)
rankstat_comb <- proportion_classified +
  plot_annotation(title = "Percentage of the number of ASVs classified at each rank") +
  plot_layout(guides = "collect",
              axis_titles = "collect")

ggsave(
  filename = "03_rankstat.png",
  plot = rankstat_comb,
  limitsize = FALSE
)
#------------------------------------------------------------------------------#
# Collect unique treatment conditions
unique_groups <- sample_data(ps_abs)$RANKSTAT_treatment %>% 
  unique()
  
group_by_rank <- ps_abs %>% 
  tax_fix(unknowns = c("uncultured")) %>% 
  merge_samples(group = "RANKSTAT_treatment") %>% 
  comp_barplot(
    tax_level = "Genus", n_taxa = 15,
    sample_order = unique_groups,
    bar_width = 0.8
  ) +
  coord_flip() + 
  labs(x = NULL, y = NULL)

ggsave(
  filename = "03_rankstat_by_rank.png",
  plot = group_by_rank,
  limitsize = FALSE
)

# Spearman correlation heatmap ------------------------------------------------#
spearman_heatmap <- cor_heatmap_plot(ps_abs, tax_level = "Genus")

# Saving heatmap --------------------------------------------------------------#
png(filename = "03_spearman_heatmap.png", width = 7, height = 7, units = 'in', res = 600)
draw(spearman_heatmap)
dev.off()

#------------------------------------------------------------------------------#
# TO DO: heatmap fold log2 change


