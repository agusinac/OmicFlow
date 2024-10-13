fetch_palette <- function(metadata, col_name, Brewer.palID="Set2") {
  # Creating color palette
  unique_groups <- unique(metadata[[col_name]])
  chosen_palette <- RColorBrewer::brewer.pal(length(unique_groups), Brewer.palID)
  colors <- stats::setNames(chosen_palette, unique_groups)
  return(colors)
}
