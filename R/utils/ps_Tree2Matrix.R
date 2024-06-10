ps_Tree2Matrix <- function(ps, tax_level = "Genus") {
  tree_mat <- ape::cophenetic.phylo(phy_tree(ps))
  taxa <- as.data.frame(tax_table(ps))[tax_level]
  tax_tab <- phyloseq::tax_table(ps) %>% as.data.frame()
  
  df_rows <- data.frame(rows = rownames(tree_mat))
  df_cols <- data.frame(cols = colnames(tree_mat))
  
  new_rows <- list()
  new_cols <- list()
  
  for (i in 1:length(df_rows$rows)) {
    new_rows[i] <- tax_tab[[ {{tax_level}} ]][stringr::str_detect(as.character(rownames(tax_tab)),
                                                                  as.character(df_rows$rows[i]))]
    new_cols[i] <- tax_tab[[ {{tax_level}} ]][stringr::str_detect(as.character(rownames(tax_tab)),
                                                                  as.character(df_cols$cols[i]))]
  }
  rownames(tree_mat) <- new_rows
  colnames(tree_mat) <- new_cols
  
  tree_mat.clean <- tree_mat[ !is.na(rownames(tree_mat)), !is.na(colnames(tree_mat))]
  tree_mat.final <- tree_mat.clean[ !grepl("uncultured", rownames(tree_mat.clean)), !grepl("uncultured", colnames(tree_mat.clean))]
  
  return(tree_mat.final)
}