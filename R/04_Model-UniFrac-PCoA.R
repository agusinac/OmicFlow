# Load data -------------------------------------------------------------------#
radboud_ps_rel_bac_norm <- readRDS("../../Pathology/projects/FFPE_breast_microbiome/data/RDS/02_ps_rel_bac_norm.rds")

# PCoA with weighted UniFrac Analysis -----------------------------------------#
comb_plot <- ps_pcoa(radboud_ps_rel_bac_norm)

ggsave(
  filename = "04_pcoa_permanova.png",
  plot = comb_plot,
  width = 10,
  height = 10,
  limitsize = FALSE
)