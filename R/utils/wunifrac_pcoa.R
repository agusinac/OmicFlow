# Takes in a phyloseq object and computes the weighted UniFrac.

# - INPUT: Phyloseq object, must include phylo tree.
# - OUTPUT: n x n plot, depending on number of plots generated and combined.

# Plots a pairwise PCoA of different dimensions that sum to at least 80%.
# Produces pairwise PERMANOVA statistics on RANKSTAT_treatment (to be modified).

wunifrac_pcoa <- function(phyloseq_obj) {
  radboud_unifrac <- phyloseq_obj %>% 
    phyloseq::distance(
      method = "wunifrac"
    )
  
  group_names <- phyloseq_obj %>% 
    phyloseq::sample_data() %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(RANKSTAT_treatment) %>% 
    as.data.frame
  
  # PCoA via SVD
  pcs <- vegan::wcmdscale(
    d = radboud_unifrac,
    k = 15,
    eig = TRUE
  )
  
 #scree_plot <- stats::screeplot(pcs)
  
  norm_pcs <- pcs$eig %>% 
    purrr::map(function(x) x / sum(pcs$eig) * 100) %>% 
    as.matrix
  
  
  df_pcs_points <- data.frame(
    pcs$points
  )
  df_pcs_points["treatment"] <- group_names
  df_pcs_points["samples"] <- row.names(df_pcs_points)
  
  # Permanova test
  permanova_results <- pairwise.adonis(as.matrix(radboud_unifrac),
                                       phyloseq::sample_data(phyloseq_obj)$RANKSTAT_treatment)
  
  # Finds dimensions required for at least 80% distance explained
  n_dimensions = 0
  sum_eig = 0
  for (eig in norm_pcs) {
    if (sum_eig < 80) {
      sum_eig <- sum_eig + eig
      n_dimensions <- n_dimensions + 1
    } else break
  }
  
  # Creates combinations and a list of plots
  n_dim_pairs <- utils::combn(seq(n_dimensions), 2)
  plot_list <- list()
  for (i in seq(ncol(n_dim_pairs))) {
    pair <- n_dim_pairs[, i]
    
    plot_title = paste0("PCoA - Dim",pair[1] ," vs Dim", pair[2])
    x_label = paste0("PC", pair[1], " ", round(as.numeric(norm_pcs[pair[1]]), 2), "%")
    y_label = paste0("PC", pair[2], " ", round(as.numeric(norm_pcs[pair[2]]), 2), "%")
    
    plt <- df_pcs_points %>% 
      ggplot(mapping = aes(x = base::get(paste0("Dim", pair[1])),
                           y = base::get(paste0("Dim", pair[2])),
                           color = treatment,
                           linetype = treatment)) +
      geom_point(alpha = 5) +
      stat_ellipse(type = "norm") +
      theme_bw() +
      labs(title = plot_title,
           subtitle = "Distances are computed via weighted UniFrac",
           x = x_label,
           y = y_label
      )
    
    plot_list[[i]] <- plt
  }
  
  # Creates permanova plot
  permanova_plot <- permanova_results %>% 
    ggplot(mapping=aes(x = pairs, 
                       y = F.Model,
                       label = p.adjusted)) +
    geom_bar(stat = "identity", 
             fill = "blue") +
    geom_label(nudge_y = 7) +
    labs(title = "PERMANOVA Results", 
         subtitle = "P adjusted significant scores are shown above each bar",
         x = "pair", 
         y = "F test statistic") +
    theme_bw()
  
  plot_list[[length(plot_list)+1]] <- permanova_plot
  #plot_list[[length(plot_list)+1]] <- scree_plot
  
  comb_plot <- patchwork::wrap_plots(plotlist = plot_list) +
    plot_layout(guides = "collect")
  
  return(comb_plot)
}
