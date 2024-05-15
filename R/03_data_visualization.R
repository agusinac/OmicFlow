# Load phyloseq objects -------------------------------------------------------#
ps_abs <- readRDS("data/RDS/01_ps_abs_bac.rds")
ps_rel <- readRDS("data/RDS/02_ps_rel_bac_norm.rds")


## Rankstat
### Classified ASVs per taxonomic rank
#------------------------------------------------------------------------------#
proportion_classified <- ps_rankstat(ps_abs) +
  plot_annotation(title = "Percentage of the number of ASVs classified at each rank") +
  plot_layout(guides = "collect",
              axis_titles = "collect")

ggsave(
  filename = "results/03_rankstat.png",
  plot = proportion_classified,
  limitsize = FALSE
)
#------------------------------------------------------------------------------#
comp_rel <- ps_composition(ps = ps_rel,
                           tax_level = "Genus",
                           metadata.columns = c("RANKSTAT_treatment", "treatment", "tumorrest"),
                           taxa_n = 10,
                           excel_path = "../results"
) +
  facet_grid(. ~ treatment + tumorrest,
             scales = "free", space = "free") +
  theme(
    axis.title.y = element_text(hjust=-0.5),
    panel.spacing.x = unit(1, "lines"),
    strip.placement = "outside",
    strip.background = element_rect(color = "black",
                                    fill = "white",
                                    linetype = "solid"),
    strip.text.x = element_text(size = 12,
                                color = "black",
                                face = "bold"))

# composition of absolute counts per sample
comp_abs <- sampleSums_plot(ps_abs, tax_level = "Genus") + theme(axis.title.y = element_text(hjust=-0.5))

comp_plot <- (comp_rel + plot_spacer() + comp_abs) +
  plot_layout(guides = "collect",
              widths = c(5, -0.12, 1.8))

ggsave(
  filename = "results/03_sample-compositon.png",
  plot = comp_rel,
  width = 15,
  height = 10,
  limitsize = FALSE
)

# Spearman correlation heatmap ------------------------------------------------#
#spearman_heatmap <- cor_heatmap_plot(ps_abs, col_name=paste0(data_00$col_name), tax_level = "Genus")

# Saving heatmap --------------------------------------------------------------#
# png(filename = "../results/03_spearman_heatmap.png", width = 7, height = 7, units = 'in', res = 600)
# draw(spearman_heatmap)
# dev.off()

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
