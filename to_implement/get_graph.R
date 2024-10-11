get_graph <- function(ps, taxa_level="Genus", top_taxa = 50, cutoff=0.5, cor_method="spearman", method_threshold=0.8) {
  # Create correlations from otu table
  otu_tab <- get_otu(ps_rel_bac, top_n = top_taxa, tax_target = taxa_level, filter_taxa = TRUE)
  taxa.cor <- abs(cor(otu_tab, method = cor_method))
  
  # Constructing graph
  net <- igraph::graph_from_adjacency_matrix(taxa.cor, mode="lower", weighted = TRUE, diag = FALSE)
  net.filt <- igraph::delete_edges(net, which(igraph::E(net)$weight <= cutoff))
  net.final <- igraph::delete_vertices(net,igraph::degree(net.filt)==0) #remove nodes without edges
  
  # Collects Taxa that meet method_threshold
  subset_correlated_vertices <- igraph::as_ids(igraph::E(net.final)[igraph::E(net.final)$weight >= method_threshold])
  subset_correlated_labels <- unique(unlist(strsplit(subset_correlated_vertices, "\\|")))
  
  # Adding color palette
  correlated_colors <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                         "#920000","#924900","#db6d00","#24ff24","#ffff6d")[1:length(subset_correlated_labels)]
  color_mapping <- setNames(correlated_colors, subset_correlated_labels)
  
  # Create a background color
  igraph::V(net.final)$color <- "lightgrey"
  
  # Replace the colors if the condition is met
  igraph::V(net.final)$color[igraph::V(net.final)$name %in% names(color_mapping)] <- color_mapping
  
  # Returns a list of data structures
  result <- list(
    network = net.final,
    legend.name = ifelse(igraph::V(net.final)$name %in% subset_correlated_labels, igraph::V(net.final)$name, NA),
    legend.color = ifelse(igraph::V(net.final)$name %in% subset_correlated_labels, igraph::V(net.final)$color, NA),
    plot.edge_color = ifelse(igraph::E(net.final)$weight >= method_threshold, "blue", "lightgrey")
  )
  return(result)
}