# Creates a spearman correlation heatmap for the specified group.
# Automatically transforms the phyloseq for k number of unique groups.

cor_heatmap_plot <- function(ps, tax_level, col_name, tax_cor="spearman") {
  # Define unique groups
  unique_groups <- phyloseq::sample_data(ps)[[ {{ col_name }} ]] %>% 
    base::unique()
  
  # Tax fixes the phyloseq object
  psq <- ps %>% 
    microViz::tax_fix(
      unknowns = c("uncultured"),
      anon_unique = FALSE,
      verbose = FALSE)
  
  # Creates groups with binary values
  for (variable in unique_groups) {
    psq <- psq %>% 
      microViz::ps_mutate(
        !!variable := if_else({{ col_name }} == variable, true = 1, false = 0)
      )
  }
  # performs a correlation heatmap plot
  plt <- psq %>% 
    microViz::tax_agg({{tax_level}}) %>%
    microViz::cor_heatmap(
      taxa = tax_top(psq, rank = {{tax_level}}, n = 15, by = max),
      cor = {{tax_cor}},
      vars = unique_groups,
      colors = heat_palette("Green-Orange", rev = TRUE, sym = TRUE)
    )
  
  return(plt)
}