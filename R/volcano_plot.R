#' Volcano plot
#'
#' @description Creates a Volcano plot from the output of \link{foldchange}, it plots the foldchanges on the x-axis, 
#' log10 trasnformed p-values on the y-axis and adjusts the scatter size based on the percentage abundance of the features. 
#' This function is built into the class \link{omics} with method \code{DFE()} and inherited by other omics classes, such as;
#' \link{metagenomics} and \link{proteomics}.
#'
#' @param data A \link[data.table]{data.table}.
#' @param logfold_col A column name of a continuous variable.
#' @param pvalue_col A column name of a continuous variable.
#' @param abundance_col A column name of a continuous variable.
#' @param feature_rank A character variable of the feature column.
#' @param logfold.threshold A Log2(A/B) Fold Change threshold (default: 0.6).
#' @param pvalue.threshold A P-value threshold (default: 0.05).
#' @param abundance.threshold An abundance threshold (default: 0.01).
#' @param label_A A character to describe condition A.
#' @param label_B A character to describe condition B.
#' 
#' @return A \link[ggplot2]{ggplot2} object to be further modified.
#' @importFrom ggplot2 ggplot aes .data theme_bw theme element_text labs geom_point geom_vline geom_hline element_rect scale_color_gradient2 scale_size_continuous
#' @examples 
#' library(data.table)
#' library(ggplot2)
#' 
#' # Create mock data frame
#' mock_volcano_data <- data.table(
#' 
#'   # Feature names (feature_rank)
#'   Feature = paste0("Gene", 1:20),   
#' 
#'   # Log2 fold changes (X)        
#'   log2FC = c(1.2, -1.5, 0.3, -0.7, 2.3,    
#'              -2.0, 0.1, 0.5, -1.0, 1.8,
#'              -0.4, 0.7, -1.4, 1.5, 0.9,
#'              -2.1, 0.2, 1.0, -0.3, -1.8),
#'   
#'   # P-values (Y)
#'   pvalue = c(0.001, 0.02, 0.3, 0.04, 0.0005, 
#'              0.01, 0.7, 0.5, 0.02, 0.0008,
#'              0.15, 0.06, 0.01, 0.005, 0.3,
#'              0.02, 0.8, 0.04, 0.12, 0.03),
#'   
#'   # Mean (relative) abundance for point sizing
#'   rel_abun = runif(20, 0.01, 0.1)            
#' )
#' 
#' volcano_plot(
#'   data = mock_volcano_data,
#'   logfold_col = "log2FC",
#'   pvalue_col = "pvalue",
#'   abundance_col = "rel_abun",
#'   feature_rank = "Feature",
#' )
#' @export

volcano_plot <- function(data,
                         logfold_col,
                         pvalue_col,
                         feature_rank,
                         abundance_col,
                         pvalue.threshold = 0.05,
                         logfold.threshold = 0.6,
                         abundance.threshold = 0.01,
                         label_A = "A",
                         label_B = "B") {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.table"))
    cli::cli_abort("Data must be a data.table.")

  if (!is.character(logfold_col) && length(logfold_col) != 1) {
    cli::cli_abort("{logfold_col} needs to contain characters with length of 1.")
  } else if (!column_exists(logfold_col, data)) {
    cli::cli_abort("The {logfold_col} column does not exist in the provided data.")
  }

  if (!is.character(pvalue_col) && length(pvalue_col) != 1) {
    cli::cli_abort("{pvalue_col} needs to contain characters with length of 1.")
  } else if (!column_exists(pvalue_col, data)) {
    cli::cli_abort("The {pvalue_col} column does not exist in the provided data.")
  }

  if (!is.character(abundance_col) && length(abundance_col) != 1) {
    cli::cli_abort("{abundance_col} needs to contain characters with length of 1.")
  } else if (!column_exists(abundance_col, data)) {
    cli::cli_abort("The {abundance_col} column does not exist in the provided data.")
  }

  if (!is.character(feature_rank) && length(feature_rank) != 1) {
    cli::cli_abort("{feature_rank} needs to contain characters with length of 1.")
  } else if (!column_exists(feature_rank, data)) {
    cli::cli_abort("The {feature_rank} column does not exist in the provided data.")
  }

  if (!is.numeric(pvalue.threshold))
    cli::cli_abort("{pvalue.threshold} need to be numeric.")

  if (!is.numeric(logfold.threshold))
    cli::cli_abort("{logfold.threshold} need to be numeric.")

  if (!is.numeric(abundance.threshold))
    cli::cli_abort("{abundance.threshold} need to be numeric.")

  if (!is.character(label_A))
    cli::cli_abort("{label_A} needs to contain characters.")

  if (!is.character(label_B))
    cli::cli_abort("{label_B} needs to contain characters.")

  ## MAIN
  #--------------------------------------------------------------------#

  # copies data.table
  tmpdt <- data.table::copy(data)

  # Creates labels for significant and non-significant differential expression
  tmpdt[, (pvalue_col) := -log10(base::get(pvalue_col))]
  tmpdt[, "diffexpressed" := ifelse(
    base::get(logfold_col) > logfold.threshold &
    base::get(pvalue_col) > -log10(pvalue.threshold) &
    base::get(abundance_col) >= abundance.threshold, 
    "Upregulated",
    ifelse(
      base::get(logfold_col) < -logfold.threshold &
      base::get(pvalue_col) > -log10(pvalue.threshold) &
      base::get(abundance_col) >= abundance.threshold, 
      "Downregulated", 
      "non-significant"
      )
    )]
  tmpdt[, "diffexpressed_labels" := ifelse(base::get("diffexpressed") != "non-significant", base::get(feature_rank), "")]

  return(
    tmpdt %>%
      ggplot(mapping = aes(x = .data [[ logfold_col ]],
                           y = .data [[ pvalue_col ]],
                           label = .data[["diffexpressed_labels"]],
                           color = .data [[ logfold_col ]])) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
            axis.text.y = element_text(size=12),
            axis.text = element_text(size=12),
            text = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF")) +
      geom_vline(xintercept = c(-logfold.threshold, logfold.threshold),
                 col = "black", linetype = 'dashed') +
      geom_hline(yintercept = -log10(pvalue.threshold),
                 col = "black", linetype = 'dashed') +
      scale_color_gradient2(name = "foldchange",
                            low = "blue",
                            mid = "black",
                            high = "red",
                            na.value = "grey80") +
      ggrepel::geom_label_repel(show.legend = FALSE,
                                max.overlaps = getOption("ggrepel.max.overlaps", default = Inf),
                                color = "black") +
      geom_point(aes(size = as.numeric(ifelse(.data[["diffexpressed"]] != "non-significant", .data[[ abundance_col ]]*100, 0))),
                 shape = 16, alpha = 0.5) +
      scale_size_continuous(name = "Mean Abundance (%)") +
      labs(x = paste0("Fold Change log2( ", label_A," / ", label_B," )"),
           y = paste0("-log10( ", pvalue_col ," )"))
  )
}
