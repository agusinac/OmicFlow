# PCoA with weighted UniFrac Analysis -----------------------------------------#

#####################
### RANKSTAT PCoA ###
#####################

if (RANKSTAT_ncol > 0) {
  metrics <- c("wunifrac", "uunifrac")
  nrow <- length(metrics)
  pcoa_plots <- matrix(list(), RANKSTAT_ncol, nrow)
  nmds_plots <- matrix(list(), RANKSTAT_ncol, nrow)
  
  for (i in 1:RANKSTAT_ncol) {
    
    col_name <- colnames(RANKSTAT_data)[i]
    
    for (j in 1:nrow) {
      
      pcoa_plots[[i, j]] <- patchwork::wrap_plots(
        ps_ordination(ps = ps_rel,
                      ordination_method = "pcoa",
                      dist.metric = metrics[j], 
                      col_name = col_name)[c("scree_plot", "scores_plot", "anova_plot")],
        nrow = 1) +
        plot_layout(widths = c(5, 5, 5),
                    guides = "collect")
      
      nmds_plots[[i, j]] <- patchwork::wrap_plots(
        ps_ordination(ps = ps_rel,
                      ordination_method = "nmds",
                      dist.metric = metrics[j], 
                      col_name = col_name)[c("scores_plot", "anova_plot")],
        nrow = 1) + 
        plot_layout(widths = c(5, 5),
                    guides = "collect")
    }
    
  }
}
    
# ggsave(
#   filename = "automated-omics-analysis/results/04_pcoa_permanova.png",
#   plot = patchwork::wrap_plots(pcoa_plots,
#                                ncol = RANKSTAT_ncol,
#                                nrow = nrow),
#   width = 30,
#   height = 10,
#   limitsize = FALSE,
#   dpi = 300,
#   scaling = 1
# )
