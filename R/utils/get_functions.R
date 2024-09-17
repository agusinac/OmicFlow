# Transform any taxa rank into a otu_table;
# with rows as samples and cols as taxa (rank to be specified)
get_otu <- function(ps, tax_target="Genus", top_n=FALSE, filter_taxa = FALSE) {
  ps_ref <- ps %>% 
    phyloseq::tax_glom(taxrank = tax_target)
  df <- ps_ref %>% 
    phyloseq::otu_table() %>% 
    t() %>% 
    as.data.frame()
  
  # Checks presence of hash code and replaces by value of tax_target
  ref_df <- data.frame(as(phyloseq::tax_table(ps_ref), "matrix"))
  tax_labels <- colnames(df)
  
  for (i in 1:length(tax_labels)) {
    if (tax_labels[i] %in% rownames(ref_df)) {
      tax_labels[i] <- ref_df[tax_labels[i], tax_target]
    }
  }
  
  colnames(df) <- tax_labels
  
  if (filter_taxa) {
    # Filtering undesired taxa by 'database_filter' file
    toMatch <- read.table(file = "C:/Users/Z289224/Documents/general_scripts/R/automated-omics-analysis/documents/database_filter")$V1
    df <- df[, !grepl(paste(toMatch, collapse = "|"), colnames(df))] 
  }
  # Fetch top taxa
  if (top_n == FALSE) {
    return(df)
  } else if (length(colnames(df)) <= top_n) {
    return(df)
  } else {
    top_taxa_names <- colSums(df) %>% 
      sort(decreasing = TRUE) %>% 
      head(top_n)
    filt_taxaSums <- df[, colnames(df) %in% names(top_taxa_names)]
    return(filt_taxaSums)
  }
}

get_tree <- function(ps, tax_target) {
  # Prune tree by otu genera
  ps_ref <- ps %>% 
    phyloseq::tax_glom(taxrank = tax_target)
  
  # Retrieve Tree
  tree <- phy_tree(ps_ref)
  
  # Create data tables to match hash and tax_target labels
  tip.labels <- tree$tip.label
  ref_df <- as.data.frame(as(phyloseq::tax_table(ps_ref), "matrix"))
  
  # Replace hash code by genus name
  for (i in 1:length(tip.labels)) {
    if (tip.labels[i] %in% rownames(ref_df)) {
      tip.labels[i] <- ref_df[tip.labels[i], tax_target]
    }
  }
  # Rewrite tree tip with tax_target labels
  tree$tip.label <- tip.labels
  # Cleans up tip labels
  tree.filt <- drop.tip(tree, tree$tip.label[grep("uncultured|unidentified", tree$tip.label)])
  
  return(tree.filt)
}

get_meta <- function(ps) {
  # Fetches meta table in right dataframe format
  meta_tab <- ps %>% 
    phyloseq::sample_data() %>% unclass() %>% as.data.frame()
  
  return(meta_tab)
}

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
