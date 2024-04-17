pairwise_triplot <- function(model, target_col, metadata, pairwise = FALSE, outdir=".", choice_dim=c("RDA1", "PC1")) {
  # Nested function #1: Find dimensions for 80% variance explained
  eigen_80 <- function(eig_explained) {
    sum_variance = 0
    counter = 0
    for (i in 1:length(eig_explained)) {
      sum_variance <- sum_variance + eig_explained[i]
      counter <- counter + 1
      if (sum_variance >= 80) break
    }
    
    return(counter)
  }
  
  # Nested function #2: Create a triplot by specific dimensions
  triplot <- function(model, target_col, metadata, perc_explained, 
                      scores_species, scores_species_explained, scores_sites, pc1, pc2) {
    #Fetch groups
    mygroups <- get(target_col, metadata)
    # Fetch number of unique groups and assign a color
    mygroups.unique <- base::unique( get(target_col, metadata) )
    mygroups.color <- rainbow(length(mygroups.unique))  
    
    # Create empty plot with labels
    graphics::plot(
      x = model,
      type = "none",
      frame = FALSE,
      main = paste0("Triplot RDA without scaling, ", pc1 , " vs ", pc2 ),
      sub = paste0("Selected column: ", target_col),
      xlab = paste0(pc1, " (", perc_explained[ pc1 ], "%)"),
      ylab = paste0(pc2, " (", perc_explained[ pc2 ], "%)")
    ) 
    
    # Add points for sites scores
    graphics::points(
      x = scores_sites[, c(pc1, pc2)],
      pch = 21,
      col = "black",
      bg = "steelblue",
      cex = 1.2
    )
    
    # Add points for species scores
    graphics::points(
      x = scores_species[, c(pc1, pc2)],
      pch = 22,
      col = "black",
      bg = "#f2bd33",
      cex = 1.2
    )
    # Add arrows based on species scores explained
    arrows(
      x0 = 0,
      y0 = 0,
      x1 = scores_species_explained[, pc1],
      y1 = scores_species_explained[, pc2],
      col = "red",
      lwd = 2
    )
    
    # Add text for top explained species scores
    text(
      x = scores_species_explained[, pc1] + 0.2,
      y = scores_species_explained[, pc2] - 0.05,
      labels = base::rownames(scores_species_explained),
      col = "red",
      cex = 0.5,
      font = 1
    )
    
    
    if (grepl("RDA", pc1, fixed = TRUE) | grepl("RDA", pc1, fixed = TRUE)) {
      subset_model <- model
    } else {
      # Subsets left (u) and right (v) ordination components
      subset_model <- model
      subset_model$CA$u <- model$CA$u[, c(pc1, pc2)]
      subset_model$CA$v <- model$CA$v[, c(pc1, pc2)]
    }

    # Add ordihull
    vegan::ordihull(subset_model, 
                    groups = mygroups,
                    show.group = mygroups,
                    col = mygroups.color,
                    draw = 'polygon',
                    alpha = 50)
    
    # Add labels of groups
    vegan::ordispider(subset_model,
                      groups = mygroups,
                      show.group = mygroups.unique,
                      col = 'black', lty = 'blank',
                      label = TRUE)
  }
  
  subset_by_dimensions <- function(model, dimensions) {
    perc_explained <- round(100*(summary(model)$cont$importance[2, dimensions]),2)
    n_dim_pairs <- utils::combn(dimensions[1:eigen_80(perc_explained)], 2)
    
    results <- list(
      var_explained = perc_explained,
      n_dim_pairs = n_dim_pairs
    )
    return(results)
  }
  
  subset_by_species <- function(model, pc) {
    species_explained <- utils::head(base::sort(round(100*scores_species[, pc]^2, 3), decreasing = TRUE))
    scores_species_explained <- scores_species[rownames(scores_species) %in% names(species_explained),]
    return(scores_species_explained) 
  }
  
  # Subset species and sites scores
  scores_species <- vegan::scores(x = model, display = "species", choices = c(1:15))
  scores_sites <- vegan::scores(x = model, display = "sites", choices = c(1:15))
  
  # Main pairwise code
  if (pairwise == TRUE) {
    # Creates a vector of 10 dimensions (PC1 - PC10)
    pairwise_dims <- sprintf("PC%d", seq(1:11))
    subset.result <- subset_by_dimensions(model, pairwise_dims)
    
    # save pdf
    pdf(paste0(outdir, "/test_pairwise.pdf"))
    for (i in seq(ncol(subset.result$n_dim_pairs))) {
      pair <- subset.result$n_dim_pairs[, i]
      scores_species_explained <- subset_by_species(model, pc = pair[1])
      
      triplot(model, target_col, metadata, subset.result$var_explained, scores_species, 
              scores_species_explained, scores_sites,
              pc1 = pair[1],
              pc2 = pair[2])
    }
    dev.off()
    
  } else {
    # Subset species most fitted/captured by user defined dimensions
    choice_dim.scores_species_explained <- subset_by_species(model, choice_dim[1])
    choice_dim.explained <- subset_by_dimensions(model, choice_dim)$var_explained
    
    # Create custom plot by user specification
    triplot(model, target_col, metadata, choice_dim.explained, scores_species, 
            choice_dim.scores_species_explained, scores_sites, choice_dim[1], choice_dim[2])
  }
  
}
