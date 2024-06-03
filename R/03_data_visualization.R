## Rankstat
### Classified ASVs per taxonomic rank
#------------------------------------------------------------------------------#
proportion_classified <- ps_rankstat(ps_abs) +
  plot_annotation(title = "Percentage of the number of ASVs classified at each rank") +
  plot_layout(guides = "collect",
              axis_titles = "collect")

ggsave(
  filename = "03_rankstat.png",
  plot = proportion_classified,
  limitsize = FALSE
)
#------------------------------------------------------------------------------#
# Microbiome composition by all samples
taxa_name <- "Genus"

res <- ps_composition(ps = ps_rel,
                      tax_level = taxa_name,
                      metadata.columns = data_00$col_name,
                      taxa_n = 10
)

comp_plot <- composition_plot(df = res$df,
                              palette = res$palette,
                              tax_level = res$tax_rank,
                              title_name = ""
)

ggsave(filename = "03_composition.png",
       plot = comp_plot,
       width = 30,
       height = 18,
       limitsize = FALSE)