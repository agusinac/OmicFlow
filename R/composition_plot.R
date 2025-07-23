#' Composition plot with ggplot2
#'
#' @description Creates a composition of features works on the output of \link[OmicFlow]{metagenomics} with method \code{composition} from abstract class \link[OmicFlow]{omics}
#'
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table}.
#' @param palette An object with names and hexcode or color names, see \link[OmicFlow]{fetch_palette}.
#' @param feature_rank A character variable of the feature column (Default: `"Genus"`).
#' @param title_name A character variable to set the \code{ggtitle} of the ggplot (Default: NULL).
#' @param group_by A character variable to aggregate the stacked bars by group.
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#'
#' @examples
#' # Initialize a new object 'taxa_class'
#' taxa_class <- metagenomics$new(metaData = "metadata.tsv",
#'                                biomData = "biom_with_taxonomy.biom",
#'                                treeData = "rooted_tree.newick")
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

composition_plot <- function(data,
                             palette,
                             feature_rank = "Genus",
                             title_name = NULL,
                             group_by = FALSE) {
  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") || !inherits(data, "data.table"))
    cli::cli_abort("Data must be a data.frame or data.table.")

  if (!is.character(palette))
    cli::cli_abort("{palette} needs to contain characters.")

  if (!is.character(feature_rank) && length(feature_rank) != 1)
    cli::cli_abort("{feature_rank} needs to contain characters with length of 1.")

  if (!is.null(title_name) && !is.character(title_name))
    cli::cli_abort("{title_name} needs to be of type character.")

  if (!column_exists("SAMPLE_ID", data))
    cli::cli_abort("SAMPLE_ID needs to exist within the provided data.frame/data.table.")

  ## MAIN
  #--------------------------------------------------------------------#

  # Generates a stacked barplot as base with custome palette
  if (group_by != FALSE) {
    plt <- data %>%
      ggplot(mapping = aes(y = value,
                           x = base::get(group_by, data),
                           fill = base::get(feature_rank, data)))
  } else {
    plt <- data %>%
      ggplot(mapping = aes(y = value,
                           x = SAMPLE_ID,
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
