# To be modified and integrated into class-tools or stays as a standalone script. To be decided!

triplot <- function(counts, metadata, metadata.col = NA, choice_dim = c("RDA1", "PC1"), pairwise = FALSE, Brewer.palID = "Set2", counts.scalar = 1) {
  # Nested functions, later to be combined within tools-class
  logn <- function(otu_tab, scalar=1) {
    # log-transform, center
    # Y' = log ( A * Y + 1 ) ; where A is the 'strength' of the log transformation : 1, 10, 100, 1000, etc., default = 1
    otu_tab.log <- ( scalar * otu_tab ) + 1 
    otu_tab.log <- log( otu_tab.log ) 
    otu_tab.sc <- scale(otu_tab.log, center = TRUE, scale = FALSE)
    return(otu_tab.sc)
  }
  
  eigen_80 <- function(eig_explained) {
    sum_variance = 0
    counter = 1
    for (i in 1:length(eig_explained)) {
      sum_variance <- sum_variance + eig_explained[i]
      counter <- counter + 1
      if (sum_variance >= 80) break
    }
    
    return(counter)
  }
  
  subset_by_dimensions <- function(model, dimensions) {
    perc_explained <- round(100*(summary(model)$cont$importance[2, dimensions]),2)
    n_dim_pairs <- dimensions[1:eigen_80(perc_explained)]
    return(perc_explained)
  }
  
  subset_by_species <- function(model, scores_species, pc) {
    species_explained <- utils::head(base::sort(round(100*scores_species[, pc]^2, 3), decreasing = TRUE))
    scores_species_explained <- scores_species[rownames(scores_species) %in% names(species_explained),]
    
    result <- list(
      scores = scores_species_explained,
      explained_PC1 = species_explained
    )
    
    return(result) 
  }
  # Main pairwise code
  if (pairwise == TRUE) {
    # Creates a vector of 10 dimensions (PC1 - PC10)
    pairwise_dims <- sprintf("PC%d", seq(1:11))
    subset.result <- subset_by_dimensions(model, pairwise_dims)
    
    # save pdf
    pdf(paste0(outdir, "/pairwise_PCA.pdf"))
    for (i in seq(ncol(subset.result$n_dim_pairs))) {
      pair <- subset.result$n_dim_pairs[, i]
      scores_species_explained <- subset_by_species(model, scores_species, pc = pair[1])
      
      triplot(model, target_col, metadata, subset.result$var_explained, scores_species, 
              scores_species_explained, scores_sites,
              pc1 = pair[1],
              pc2 = pair[2])
    }
    dev.off()
  }
  # Subsets user specified dimensions
  pc1 <- choice_dim[1]
  pc2 <- choice_dim[2]
  
  # Transformation of counts and modelling to RDA
  counts.log <- logn(counts, scalar = counts.scalar)
  model <- vegan::rda(counts.log ~ get(metadata.col, metadata) + Condition(NULL), 
                      data = metadata, 
                      scale = FALSE, 
                      na.action = na.fail, 
                      subset = NULL)
  
  # MAIN code
  # Subset species and sites scores
  model.dims <- dim(model$CCA$u)[2] + dim(model$CA$u)[2]
  scores_species <- vegan::scores(x = model, display = "species", choices = c(1:model.dims), scaling=0)
  scores_sites <- vegan::scores(x = model, display = "sites", choices = c(1:model.dims))
  
  # Subset species most fitted/captured by user defined dimensions
  choice_dim.scores_species_explained <- subset_by_species(model, scores_species, pc = choice_dim[1])
  choice_dim.explained <- subset_by_dimensions(model, choice_dim)
  
  # Include relative abundance and significant groups in scores_sites
  rel_abun <- colSums(counts)
  Explained_species <- rownames(choice_dim.scores_species_explained$scores)
  scores_species_merged <- data.table::data.table(cbind(scores_species, rel_abun))
  scores_species_merged$taxa <- rownames(scores_species)
  
  # Creating color palette
  chosen_palette <- RColorBrewer::brewer.pal(length(Explained_species), Brewer.palID)
  colors <- stats::setNames(chosen_palette, Explained_species)
  
  # include groups for labelling and size
  scores_species_merged[, explained_species_label := ifelse(taxa %in% Explained_species, taxa, "")]
  scores_species_merged[, explained_species_size := ifelse(taxa %in% Explained_species, rel_abun, 0)]
  
  #Fetch groups
  mygroups <- get(metadata.col, metadata)
  
  # to be named: scores_sites
  dt <- data.table::data.table(data.frame(pc1 = scores_sites[, pc1],
                                          pc2 = scores_sites[, pc2],
                                          group = mygroups))
  
  # Get centroid centers for annotation
  df_mean.ord <- stats::aggregate(dt, by=list(dt$group),mean)
  colnames(df_mean.ord) <- c("Group", "x", "y")
  df_mean.ord <- df_mean.ord[df_mean.ord$Group %in% mygroups, ]
  
  rslt.hull <- vegan::ordihull(scores_sites[, c(pc1, pc2)], 
                               groups = mygroups,
                               show.group = mygroups,
                               draw = 'none')
  
  # Initialize an empty data.table
  df_hull <- data.table::data.table(Group = character(), x = numeric(), y = numeric())
  
  # Loop through groups and bind data
  for (g in mygroups) {
    x <- rslt.hull[[g]][ , 1]
    y <- rslt.hull[[g]][ , 2]
    Group <- rep(g, length(x))
    df_temp <- data.table::data.table(Group = Group, x = x, y = y)
    df_hull <- rbind(df_hull, df_temp, use.names = TRUE, fill = TRUE)
  }
  
  # Convert to data.table
  data.table::setDT(df_hull)
  
  # Returns plot
  return(
    ggplot() +
      geom_point(data = scores_sites, aes(x = .data[[pc1]], 
                                          y = .data[[pc2]]),
                 shape = 21,
                 fill = "steelblue",
                 col = "black") +
      geom_point(data = scores_species_merged, aes(x = .data[[pc1]], 
                                                   y = .data[[pc2]],
                                                   size = .data[["explained_species_size"]],
                                                   col = .data[["explained_species_label"]]),
                 show.legend = TRUE,
                 stroke = ifelse(scores_species_merged$explained_species_label != "", 1.5, 0.5),
                 shape = 22) +
      geom_polygon(data=df_hull, aes(x=x, y=y, fill=Group),
                   alpha = 0.2, color = "gray40", show.legend = FALSE) +
      geom_text(data=df_mean.ord, aes(x=x, y=y, label=Group), 
                color = "black", show.legend = FALSE) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
            axis.text.y = element_text(size=12),
            axis.text = element_text(size=12),
            text = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
      scale_size_continuous(name = "rel. Abun.") +
      scale_color_manual(name = paste0(pc1, " explained species"),
                         values = colors) +
      labs(x = paste0(pc1, " (", choice_dim.explained[ pc1 ], "%)"),
           y = paste0(pc2, " (", choice_dim.explained[ pc2 ], "%)")) +
      guides(fill = "none", colour = guide_legend(override.aes = list(stroke = 1.5)))
  )
}