#' Composition plot with ggplot2
#'
#' @description Creates a composition of features works on the output of \link[OmicFlow]{metataxonomics} with method \code{composition} from abstract class \link[OmicFlow]{tools}
#'
#' @param data A \code{data.frame} or \code{data.table}.
#' @param palette An object with names and hexcode or color names, see \link{fetch_palette} or \link[stats]{setNames}.
#' @param feature_rank A character variable of the feature column, by default set to "Genus".
#' @param title_name A character variable to set the \code{ggtitle} of the ggplot, by default set to NULL.
#' @param group_by A character variable to aggregate the stacked bars by group.
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#'
#' @examples
#' # Initialize a new object 'taxa_class'
#' taxa_class <- metataxonomics$new(metaData = "metadata.tsv",
#'                                  biomData = "biom_with_taxonomy.biom",
#'                                  treeData = "rooted_tree.newick")
#'
#' # Compute the composition for the top 10 features
#' result <- taxa_class$composition(feature_rank = "Genus",
#'                                  feature_filter = c("uncultured"),
#'                                  feature_top = 10)
#'
#' # Create a ggplot() graph with composition_plot
#' composition_plot(data = result$data,
#'                  palette = result$palette,
#'                  feature_rank = "Genus")
#'
#' @export

composition_plot <- function(data, palette, feature_rank = "Genus", title_name = NULL, group_by = FALSE) {
  # Generates a stacked barplot as base with custome palette
  if (group_by != FALSE) {
    plt <- data %>%
      ggplot(mapping = aes(y = value,
                           x = base::get(group_by, data),
                           fill = base::get(feature_rank, data)))
  } else {
    plt <- data %>%
      ggplot(mapping = aes(y = value,
                           x = `SAMPLE-ID`,
                           fill = base::get(feature_rank, data)))
  }
  # Required for stacked barplot
  plt <- plt +
    geom_bar(position = "fill",
             stat = "identity")

  if (group_by == FALSE) {
    plt <- plt +
      coord_flip()
  }
  plt <- plt +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 13),
      axis.text.x = element_text(angle = 90, size = 12,
                                 vjust = 0.5, hjust=1,
                                 colour = "black"),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12, vjust=0.5),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12, colour = "black"),
      axis.text.y = element_text(colour = "black", size = 12)
    )

  if (group_by == FALSE) {
    plt <- plt +
      scale_x_discrete(limits = rev(levels(as.factor(data$`SAMPLE-ID`))))
  }
  plt <- plt +
    scale_fill_manual(values = palette, name = feature_rank) +
    labs(y = "Rel. Abun.",
         x = NULL) +
    ggtitle(title_name)

  return(plt)
}
