# PCoA with weighted UniFrac Analysis -----------------------------------------#
comb_plot <- ps_pcoa(ps = ps_rel, dist.metric = "wunifrac", col_name = data_00$col_name)

ggsave(
  filename = "results/04_pcoa_permanova.png",
  plot = comb_plot,
  width = 15,
  height = 5,
  limitsize = FALSE
)