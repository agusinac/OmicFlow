# Load phyloseq objects -------------------------------------------------------#
ps_abs <- readRDS("../data/RDS/01_ps_abs_bac.rds")
ps_rel <- readRDS(paste0(dir, "../data/RDS/02_ps_rel_bac_norm.rds"))


## Rankstat
### Classified ASVs per taxonomic rank
#------------------------------------------------------------------------------#
proportion_classified <- ps_rankstat(ps_abs) +
  plot_annotation(title = "Percentage of the number of ASVs classified at each rank") +
  plot_layout(guides = "collect",
              axis_titles = "collect")

ggsave(
  filename = "../results/03_rankstat.png",
  plot = proportion_classified,
  limitsize = FALSE
)
#------------------------------------------------------------------------------#
# Microbiome composition by all samples
comp_rel <- ps_composition(ps = ps_rel,
                           tax_level = "Genus",
                           metadata.columns = c(paste0(data_00$col_name)),
                           taxa_n = 10,
                           excel_path = FALSE
) +
  theme(
    panel.spacing.x = unit(1, "lines"),
    strip.placement = "outside",
    strip.background = element_rect(color = "black",
                                    fill = "white",
                                    linetype = "solid"),
    strip.text.x = element_text(size = 12,
                                color = "black",
                                face = "bold"))

ggsave(
  filename = paste0(dir, "../results/03_sample-compositon.png"),
  plot = comp_rel,
  width = 10,
  height = 10,
  limitsize = FALSE
)

# Spearman correlation heatmap ------------------------------------------------#
spearman_heatmap <- cor_heatmap_plot(ps_abs, col_name=data_00$col_name, tax_level = "Genus")

# Saving heatmap --------------------------------------------------------------#
png(filename = "../results/03_spearman_heatmap.png", width = 7, height = 7, units = 'in', res = 600)
draw(spearman_heatmap)
dev.off()

#------------------------------------------------------------------------------#

# ps <- ps_rel %>% 
#   subset_samples({{ data_00$col_name }} != "control") %>% 
#   subset_taxa(Genus != "Pseudomonas") %>% 
#   transform_sample_counts(function(x) x / sum(x)) %>% 
#   removeZeros()
# 
# plt_heatmap <- ps_heatmap(ps, col_name = data_00$col_name)
# 
# ggsave(
#   filename = "../results/03_heatmap_2fold.png",
#   plot = plt_heatmap,
#   limitsize = FALSE,
#   width = 10,
#   height = 10
# )
