
# Takes in a phyloseq object and computes the weighted UniFrac.

# - INPUT: Phyloseq object, must include phylo tree.
# - OUTPUT: n x n plot, depending on number of plots generated and combined.

# Plots a pairwise PCoA of different dimensions that sum to at least 80%.
# Produces pairwise PERMANOVA statistics on RANKSTAT_treatment (to be modified).

ps_pcoa <- function(ps, dist.metric, col_name, pairwise = FALSE, dist.explained=80, choice_dim=c(1,2), outdir="results") {
  # General PCoA plot via ggplot
  pcoa_plot <- function(df_pcs_points, pcs, pair, dist.metric) {
    plot_title = paste0("PCoA - PC",pair[1] ," vs PC", pair[2])
    x_label = paste0("PC", pair[1], " ", round(as.numeric(pcs$eig_norm[pair[1]]), 2), "%")
    y_label = paste0("PC", pair[2], " ", round(as.numeric(pcs$eig_norm[pair[2]]), 2), "%")
    
    plt <- df_pcs_points %>% 
      ggplot(mapping = aes(x = base::get(paste0("Dim", pair[1])),
                           y = base::get(paste0("Dim", pair[2])),
                           color = groups,
                           linetype = groups)) +
      geom_point(alpha = 5) +
      stat_ellipse(type = "norm") +
      theme_bw() +
      labs(title = plot_title,
           subtitle = paste0("Distances are computed via ", dist.metric, " metric"),
           x = x_label,
           y = y_label
      )
    return(plt)
  }
  
  # Switch case to compute different dissimilarity matrices
  dist_mat <- switch(
    dist.metric,
    "wunifrac" = if (!is.null(phy_tree(ps))) {
      phyloseq::distance(ps, method = "wunifrac")},
    "uunifrac" = if (!is.null(phy_tree(ps))) {
      phyloseq::distance(ps, method = "uunifrac")},
    "bray" = phyloseq::distance(ps, method = "bray"),
    "jaccard" = phyloseq::distance(ps, method = "jaccard"),
    "jsd" = phyloseq::distance(ps, method = "jsd")
  )
  
  # Fetch sample groups
  group_names <- sample_data(ps)[[ {{col_name}} ]]
  
  # PCoA via Singular Value Decomposition (SVD)
  pcs <- vegan::wcmdscale(
    d = dist_mat,
    k = 15,
    eig = TRUE
  )
  
  # Normalized eigenvalues
  pcs$eig_norm <- pcs$eig %>% 
    purrr::map(function(x) x / sum(pcs$eig) * 100) %>% 
    unlist()
  
  # Collects loading scores into dataframe
  df_pcs_points <- data.frame(
    pcs$points
  )
  # Adds relevant data
  df_pcs_points["groups"] <- group_names
  df_pcs_points["samples"] <- row.names(df_pcs_points)
  
  plot_list <- list()
  # Pairwise dimensions
  if (pairwise) {
    # Finds number of dimensions that explain 80% of distances
    n_dimensions = 0
    sum_eig = 0
    for (eig in pcs$eig_norm) {
      if (sum_eig < dist.explained) {
        sum_eig <- sum_eig + eig
        n_dimensions <- n_dimensions + 1
      } else break
    }
    
    # Creates paired combinations of dimensions into a list of plots
    n_dim_pairs <- utils::combn(seq(n_dimensions), 2)
    pdf(paste0(outdir, "/pairwise_PCoA.pdf"))
    for (i in seq(ncol(n_dim_pairs))) {
      pair <- n_dim_pairs[, i]
      print(pcoa_plot(df_pcs_points, pcs, pair, dist.metric))
    }
    dev.off()
  }
  
  # Scree plot of first 10 dimensions
  plot_list[[1]] <- data.frame(
    dims = seq(length(pcs$eig_norm[1:10])),
    dims.explained = pcs$eig_norm[1:10]
  ) %>% 
    ggplot(mapping = aes(x = dims,
                         y = dims.explained)) +
    geom_col() +
    theme_bw() +
    scale_x_continuous(breaks=seq(1, 10, 1)) +
    scale_y_continuous(breaks=seq(0, 100, 10)) +
    labs(title = "Screeplot of first 10 PCs",
         x = "Principal Components (PCs)",
         y = "dissimilarity explained [%]")
  
  
  # PCoA plot according to user specified dimensions
  plot_list[[2]] <- pcoa_plot(df_pcs_points, pcs, pair=choice_dim, dist.metric)
  
  # Permanova test
  permanova_results <- pairwise.adonis(as.matrix(dist_mat),
                                       phyloseq::sample_data(ps)[[ {{ col_name }} ]])
  
  # Creates permanova plot
  plot_list[[3]] <- permanova_results %>% 
    ggplot(mapping=aes(x = pairs, 
                       y = F.Model,
                       label = p.adjusted)) +
    geom_bar(stat = "identity", 
             fill = "blue") +
    geom_label(nudge_y = 0) +
    labs(title = "PERMANOVA Results", 
         subtitle = "P adjusted significant scores are shown above each bar",
         x = "groups", 
         y = "F test statistic") +
    theme_bw()
  
  # Combines plots
  comb_plot <- patchwork::wrap_plots(plotlist = plot_list,
                                     ncol = 3,
                                     nrow = 1) +
    plot_layout(guides = "collect")
  
  return(comb_plot)
}