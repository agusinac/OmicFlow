# Load data -------------------------------------------------------------------#
radboud_ps_rel_bac_norm <- readRDS("../data/RDS/02_ps_rel_bac_norm.rds")

# PCoA with weighted UniFrac Analysis -----------------------------------------#
comb_plot <- ps_pcoa(ps = radboud_ps_rel_bac_norm, dist.metric = "wunifrac", col_name = data_00$col_name)

ggsave(
  filename = "../results/04_pcoa_permanova.png",
  plot = comb_plot,
  width = 10,
  height = 10,
  limitsize = FALSE
)