# Load phyloseq objects -------------------------------------------------------#
ps_abs <- readRDS("../data/RDS/01_ps_abs_bac.rds")
ps_rel <- readRDS("../data/RDS/02_ps_rel_bac_norm.rds")


## Rankstat
### Classified ASVs per taxonomic rank
#------------------------------------------------------------------------------#
proportion_classified <- ps_rankstat(ps_abs)
rankstat_comb <- proportion_classified +
  plot_annotation(title = "Percentage of the number of ASVs classified at each rank") +
  plot_layout(guides = "collect",
              axis_titles = "collect")

ggsave(
  filename = "../results/03_rankstat.png",
  plot = rankstat_comb,
  limitsize = FALSE
)
#------------------------------------------------------------------------------#
# Collect unique treatment conditions
unique_groups <- sample_data(ps_abs)[[ {{ data_00$col_name }} ]] %>% 
  unique()
  
group_by_rank <- ps_abs %>% 
  tax_fix(unknowns = c("uncultured")) %>% 
  merge_samples(group = {{ data_00$col_name }} ) %>% 
  comp_barplot(
    tax_level = "Genus", n_taxa = 15,
    sample_order = unique_groups,
    bar_width = 0.8
  ) +
  coord_flip() + 
  labs(x = NULL, y = NULL, caption = paste0("metadata column selected: ", data_00$col_name))

ggsave(
  filename = "../results/03_rankstat_by_rank.png",
  plot = group_by_rank,
  limitsize = FALSE
)

# Spearman correlation heatmap ------------------------------------------------#
spearman_heatmap <- cor_heatmap_plot(ps_abs, col_name=data_00$col_name, tax_level = "Genus")

# Saving heatmap --------------------------------------------------------------#
png(filename = "../results/03_spearman_heatmap.png", width = 7, height = 7, units = 'in', res = 600)
draw(spearman_heatmap)
dev.off()

#------------------------------------------------------------------------------#
# For now hardcoded, should be later replaced

ps <- ps_rel %>% 
  subset_samples({{ data_00$col_name }} != "control") %>% 
  subset_taxa(Genus != "Pseudomonas") %>% 
  transform_sample_counts(function(x) x / sum(x)) %>% 
  removeZeros()

plt_heatmap <- ps_heatmap(ps, col_name = data_00$col_name)

ggsave(
  filename = "../results/03_heatmap_2fold.png",
  plot = plt_heatmap,
  limitsize = FALSE,
  width = 10,
  height = 10
)
