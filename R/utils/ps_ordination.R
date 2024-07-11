# Takes in a phyloseq object and computes the weighted UniFrac.

# - INPUT: Phyloseq object, must include phylo tree.
# - OUTPUT: n x n plot, depending on number of plots generated and combined.

# Plots a pairwise PCoA of different dimensions that sum to at least 80%.
# Produces pairwise PERMANOVA statistics on RANKSTAT_treatment (to be modified).

ps_ordination <- function(ps, dist.metric, ordination_method, col_name, pairwise = FALSE, dist.explained=80, choice_dim=c(1,2), outdir="results") {
  
  # General PCoA plot via ggplot
  ordination_plot <- function(df_pcs, pcs, pair, dist.metric, col_name) {
    # Creates labels
    plot_title = paste0("Distance metric used: ", dist.metric)
    if (!is.null(pcs$eig_norm)) {
      x_label = paste0(pair[1], " (", round(as.numeric(pcs$eig_norm[1]), 2), "%)")
      y_label = paste0(pair[2], " (", round(as.numeric(pcs$eig_norm[2]), 2), "%)")
    } else {
      x_label = paste0(pair[1])
      y_label = paste0(pair[2])
    }
    
    return(
      df_pcs %>% 
        ggplot(mapping = aes(x = base::get(pair[1], df_pcs),
                             y = base::get(pair[2], df_pcs),
                             color = groups,
                             linetype = groups)) +
        geom_point(alpha = 5) +
        stat_ellipse(type = "norm") +
        theme_bw() +
        theme(text=element_text(size=14),
              legend.text = element_text(size=12),
              legend.title = element_text(size=14),
              axis.text = element_text(size=12),
              axis.text.y = element_text(size=12),
              axis.text.x = element_text(size=12)
        ) +
        labs(title = NULL,
             subtitle = NULL,
             x = x_label,
             y = y_label
        )
    )
  }
  
  # Universal plot for PERMANOVA/ANOSIM/ANOVA etc.
  stats_plot <- function(stats_df, X, Y, Label, Y_title, plot.title) {
    return(
      stats_df %>%
        ggplot(mapping=aes(x = base::get(X, stats_df),
                           y = base::get(Y, stats_df),
                           label = base::get(Label, stats_df))) +
        geom_bar(stat = "identity",
                 fill = "blue") +
        geom_label(nudge_y = 0) +
        labs(title = plot.title,
             subtitle = "P adjusted significant scores are shown above each bar",
             x = "groups",
             y = Y_title) +
        theme_bw()
    )
  }
  
  # Switch case to compute (dis)similarity matrices
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
  
  # Switch case to compute loading scores
  pcs <- switch(
    ordination_method,
    "pcoa" = vegan::wcmdscale(d = dist_mat, 
                              k = 15, 
                              eig = TRUE),
    "nmds" = vegan::metaMDS(dist_mat, 
                            trace = FALSE, 
                            autotransform = FALSE)
  )
  
  # Switch case to compute relevant statistics
  stats_results <- switch(
    ordination_method,
    "pcoa" = pairwise.adonis(as.matrix(dist_mat),
                             phyloseq::sample_data(ps)[[ {{ col_name }} ]]),
    "nmds" = pairwise.anosim(as.matrix(dist_mat),
                             phyloseq::sample_data(ps)[[ {{ col_name }} ]])
  )
  
  # Normalization of eigenvalues
  if (ordination_method == "pcoa") {
    pcs$eig_norm <- pcs$eig %>% 
      purrr::map(function(x) x / sum(pcs$eig) * 100) %>% 
      unlist()
    
    # Collects loading scores into dataframe
    df_pcs_points <- data.frame(pcs$points)
    colnames(df_pcs_points) <- base::sub("Dim", "PC", colnames(df_pcs_points))
  } else if (ordination_method == "nmds") {
    df_pcs_points <- data.frame(pcs$points)
    df_pcs_points$stress <- pcs$stress
  }
  
  # Adds relevant data
  df_pcs_points["groups"] <- sample_data(ps)[[ {{col_name}} ]]
  df_pcs_points["samples"] <- row.names(df_pcs_points)
  
  # Pairwise dimensions
  if (pairwise & ordination_method == "pcoa") {
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
  
  # Creates a list of plots
  plot_list <- list(scree_plot = NULL,
                    anova_plot = NULL,
                    scores_plot = NULL)
  if (ordination_method == "pcoa") {
    # Scree plot of first 10 dimensions
    plot_list$scree_plot <- data.frame(
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
    
    # PERMANOVA
    plot_list$anova_plot <- stats_plot(stats_results, 
                                       X = "pairs",
                                       Y = "F.Model",
                                       Label = "p.adjusted",
                                       Y_title = "Pseudo F test statistic",
                                       plot.title = "PERMANOVA")
    # Loading score plot
    plot_list$scores_plot <- ordination_plot(df_pcs_points, 
                                             pcs, 
                                             pair=c("PC1", "PC2"), 
                                             dist.metric, 
                                             col_name)
    
  } else if (ordination_method == "nmds") {
    plot_list$anova_plot <- stats_plot(stats_results, 
                                       X = "pairs",
                                       Y = "anosimR",
                                       Label = "p.adj",
                                       Y_title = "ANOSIM R statistic",
                                       plot.title = "ANOSIM")
    
    plot_list$scores_plot <- ordination_plot(df_pcs_points, 
                                             pcs, 
                                             pair=c("MDS1", "MDS2"), 
                                             dist.metric, 
                                             col_name)
  }
  
  return(plot_list)
}