#' Diversity plot
#' 
#' @description Creates an Alpha diversity plot. This function is built into the class \link{omics} with method \code{alpha_diversity()}.
#' It computes the pairwise wilcox test, paired or non-paired, given a data frame and adds useful labelling.
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table} computed from \link{diversity}.
#' @param values A column name of a continuous variable.
#' @param col_name A column name of a categorical variable.
#' @param palette An object with names and hexcode or color names, see \link{colormap}.
#' @param method A character variable indicating what method is used to compute the diversity.
#' @param paired A boolean value to perform paired analysis in \link[stats]{wilcox.test}.
#' @param p.adjust.method A character variable to specify the p.adjust.method to be used (Default: fdr).
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#'
#' @importFrom ggplot2 ggplot aes .data theme_bw theme element_text scale_colour_manual labs
#' @importFrom stats p.adjust.methods
#' 
#' @examples
#' library("ggplot2")
#'  
#' n_row <- 1000
#' n_col <- 100
#' density <- 0.2
#' num_entries <- n_row * n_col
#' num_nonzero <- round(num_entries * density)
#' 
#' set.seed(123)
#' positions <- sample(num_entries, num_nonzero, replace=FALSE)
#' row_idx <- ((positions - 1) %% n_row) + 1
#' col_idx <- ((positions - 1) %/% n_row) + 1
#' 
#' values <- runif(num_nonzero, min = 0, max = 1)
#' sparse_mat <- Matrix::sparseMatrix(
#'    i = row_idx,
#'    j = col_idx,
#'    x = values,
#'    dims = c(n_row, n_col)
#'  )
#' 
#' div <- OmicFlow::diversity(
#'   x = sparse_mat,
#'   metric = "shannon"
#' )
#' 
#' dt <- data.table::data.table(
#'   "values" = div,
#'   "treatment" = c(rep("healthy", n_col / 2), rep("tumor", n_col / 2))
#' )
#' 
#' colors <- OmicFlow::colormap(dt, "treatment")
#' 
#' diversity_plot(
#'  data = dt,
#'  values = "values",
#'  col_name = "treatment",
#'  palette = colors,
#'  method = "shannon",
#'  paired = FALSE,
#'  p.adjust.method = "fdr"
#' )
#' @export

diversity_plot <- function(data,
                           values,
                           col_name,
                           palette,
                           method,
                           paired = FALSE,
                           p.adjust.method = "fdr") {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("data must be a data.frame or data.table.")

  if (!is.character(palette))
    cli::cli_abort("{palette} needs to contain characters.")

  if (!is.character(method) && length(method) != 1)
    cli::cli_abort("{method} needs to contain characters with length of 1.")

  if (!is.character(values) && length(values) != 1) {
    cli::cli_abort("Column name: {values} needs to contain characters with length of 1.")
  } else if (!column_exists(values, data)) {
    cli::cli_abort("The {values} column does not exist in the provided data.")
  }

  if (!is.character(col_name) && length(col_name) != 1) {
    cli::cli_abort("values needs to contain characters with length of 1.")
  } else if (!column_exists(col_name, data)) {
    cli::cli_abort("The {col_name} column does not exist in the provided data.")
  }

  if (!c(p.adjust.method %in% p.adjust.methods))
    cli::cli_abort("Specified {p.adjust.method} is not valid. \nValid options: {p.adjust.methods}")

  ## MAIN
  #--------------------------------------------------------------------#

  result <- list()

  pvalues_adjusted <- data %>%
    rstatix::pairwise_wilcox_test(formula = stats::reformulate(col_name, response = values),
                                  p.adjust.method = p.adjust.method,
                                  paired = paired) %>%
    rstatix::add_significance() %>%
    rstatix::add_xy_position(x = col_name)

  pvalues_adjusted.filtered <- pvalues_adjusted[grepl("\\*", pvalues_adjusted$p.adj.signif) ,]

  plt <- data %>%
    ggplot(mapping = aes(x = as.factor(.data[[ col_name ]]),
                         y = .data[[ values ]])) +
    gghalves::geom_half_boxplot() +
    gghalves::geom_half_point_panel(aes(color = as.factor(.data[[ col_name ]]))) +
    theme_bw() +
    theme(legend.position = "none",
          text=element_text(size=14),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14),
          axis.text = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12)) +
    scale_colour_manual(name = "groups",
                        values = palette) +
    ggpubr::stat_pvalue_manual(pvalues_adjusted.filtered,
                               label = "p.adj",
                               step.increase = 0.05) +
    labs(title = NULL,
         subtitle = paste0(
          "Attribute: ", col_name,
          ", test: ", ifelse(paired, "Wilcox signed rank test", "Mann-Whitney U test"),
          ", p.adjusted by ", p.adjust.method),
         x = "sample groups",
         y = paste0("Alpha diversity metric: ", method))

  result <- list(
    plot = plt,
    stats = pvalues_adjusted
  )

  return(result)
}
