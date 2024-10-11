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