#' Compositional plot
#'
#' @description Creates a stacked barchart of features. It is possible to both show barcharts for each sample or group them by a categorical variable.
#' The function is compatible with the class \link{omics} method \code{composition()}.
#'
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table}.
#' @param palette An object with names and hexcode or color names, see \link{colormap}.
#' @param feature_rank A character variable of the feature column.
#' @param title_name A character to set the \code{ggtitle} of the \link[ggplot2]{ggplot}, (Default: NULL).
#' @param group_by A character variable to aggregate the stacked bars by group (Default: NULL).
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#'
#' @importFrom ggplot2 ggplot aes .data geom_bar coord_flip theme_bw theme element_text scale_x_discrete scale_fill_manual labs ggtitle
#' 
#' @examples
#' library("ggplot2")
#' 
#' # Create mock_data as data.frame (data.table is also supported)
#' mock_data <- data.frame(
#'   SAMPLE_ID = rep(paste0("Sample", 1:10), each = 5),
#'   Genus = rep(c("GenusA","GenusB","GenusC","GenusD","GenusE"), times = 10),
#'   value = c(
#'     0.1119, 0.1303, 0.0680, 0.5833, 0.1065,      # Sample1
#'     0.2080, 0.1179, 0.0211, 0.4578, 0.1951,      # Sample2
#'     0.4219, 0.1189, 0.2320, 0.1037, 0.1235,      # Sample3
#'     0.4026, 0.0898, 0.1703, 0.1063, 0.2309,      # Sample4
#'     0.1211, 0.0478, 0.5721, 0.1973, 0.0618,      # Sample5
#'     0.2355, 0.0293, 0.2304, 0.1520, 0.3528,      # Sample6
#'     0.2904, 0.0347, 0.3651, 0.0555, 0.2544,      # Sample7
#'     0.4138, 0.0299, 0.0223, 0.4996, 0.0345,      # Sample8
#'     0.4088, 0.0573, 0.0155, 0.2888, 0.2296,      # Sample9
#'     0.4941, 0.0722, 0.2331, 0.1023, 0.0983       # Sample10
#'   ),
#'   Group = rep(c("Group1","Group2","Group1",
#'                 "Group1","Group2","Group2",
#'                  "Group1","Group1","Group1","Group2"),
#'                each = 5)
#' )
#' 
#' # Create a colormap
#' mock_palette <- c(
#'   GenusA = "#1f77b4",  # blue
#'   GenusB = "#ff7f0e",  # orange
#'   GenusC = "#2ca02c",  # green
#'   GenusD = "#d62728",  # red
#'   GenusE = "#9467bd"   # purple
#' )
#' 
#' # Optionally: Use OmicFlow::colormap()
#' mock_palette <- colormap(
#'   data = mock_data,
#'   col_name = "Genus",
#'   Brewer.palID = "RdYlBu"
#' )
#' 
#' composition_plot(
#'   data = mock_data,
#'   palette = mock_palette,
#'   feature_rank = "Genus",
#'   title_name = "Mock Genus Composition"
#' )
#' 
#' composition_plot(
#'   data = mock_data,
#'   palette = mock_palette,
#'   feature_rank = "Genus",
#'   title_name = "Mock Genus Composition by Group",
#'   group_by = "Group"
#' )
#' @export

composition_plot <- function(data,
                             palette,
                             feature_rank,
                             title_name = NULL,
                             group_by = NULL) {
  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("Data must be a data.frame or data.table.")

  if (!is.character(palette))
    cli::cli_abort("{palette} needs to contain characters.")

  if (!is.character(feature_rank) && length(feature_rank) != 1)
    cli::cli_abort("{feature_rank} needs to contain characters with length of 1.")

  if (!is.null(title_name) && !is.character(title_name))
    cli::cli_abort("{title_name} needs to be of type character.")

  if (!is.null(group_by)) {
    if (!is.character(group_by) && length(group_by) != 1) {
      cli::cli_abort("{group_by} must be a character and of length 1")
    } else if (!column_exists(group_by, data)) {
      cli::cli_abort("The specified {group_by} does not exist in the metaData.")
    }
  }

  if (!column_exists("SAMPLE_ID", data))
    cli::cli_abort("SAMPLE_ID needs to exist within the provided data.frame/data.table.")

  ## MAIN
  #--------------------------------------------------------------------#

  # Generates a stacked barplot as base with custome palette
  if (!is.null(group_by)) {
    plt <- data %>%
      ggplot(mapping = aes(y = .data[["value"]],
                           x = as.factor(.data[[ group_by ]]),
                           fill = .data[[ feature_rank ]]))
  } else {
    plt <- data %>%
      ggplot(mapping = aes(y = .data[["value"]],
                           x = .data[["SAMPLE_ID"]],
                           fill = base::get(feature_rank, data)))
  }
  # Required for stacked barplot
  plt <- plt +
    geom_bar(
      position = "fill",
      stat = "identity"
    )

  if (is.null(group_by)) {
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

  if (is.null(group_by)) {
    plt <- plt +
      scale_x_discrete(limits = rev(levels(as.factor(data[["SAMPLE_ID"]]))))
  }
  plt <- plt +
    scale_fill_manual(values = palette, name = feature_rank) +
    labs(y = "Rel. Abun.",
         x = NULL) +
    ggtitle(title_name)

  return(plt)
}
