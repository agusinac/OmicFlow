#' Diversity plot
#' 
#' @description Creates an Alpha diversity plot. This function is built into the class \link{omics} with method \code{alpha_diversity()}.
#' It computes the pairwise wilcox test, paired or non-paired, given a data frame and adds useful labelling.
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table} computed from \link{diversity}.
#' @param values A column name of a continuous variable.
#' @param col_name A column name of a categorical variable.
#' @param group_by A column name to perform grouped statistical test (default: NULL).
#' @param palette An object with names and hexcode or color names, see \link{colormap}.
#' @param method A character variable indicating what method is used to compute the diversity.
#' @param paired A boolean value to perform paired analysis in \link[stats]{wilcox.test}.
#' @param p.adjust.method A character variable to specify the p.adjust.method to be used (Default: fdr).
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#'
#' @importFrom ggplot2 ggplot aes .data theme_bw theme element_text scale_colour_manual scale_x_continuous labs geom_boxplot geom_segment geom_point facet_wrap position_nudge 
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
#' plt <- OmicFlow::diversity_plot(
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
                           group_by = NULL,
                           palette,
                           method,
                           paired = FALSE,
                           p.adjust.method = "fdr") {

  ## Error handling
  #--------------------------------------------------------------------#

  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("Data must be a {.cls data.frame} or {.cls data.table}.")

  if (!is.character(palette))
    cli::cli_abort("{.val {palette}} needs to contain characters.")

  if (!is.character(method) && length(method) != 1)
    cli::cli_abort("{.val {method}} needs to contain characters with length of 1.")

  if (!is.character(values) && length(values) != 1) {
    cli::cli_abort("{.val {values}} needs to contain characters with length of 1.")
  } else if (!column_exists(values, data)) {
    cli::cli_abort("The {.val values}} column does not exist in the provided {.arg data}.")
  }

  if (!is.character(col_name) && length(col_name) != 1) {
    cli::cli_abort("{.val {col_name}} needs to contain characters with length of 1.")
  } else if (!column_exists(col_name, data)) {
    cli::cli_abort("The {.val {col_name}} column does not exist in the provided {.arg data}.")
  }

  if (!is.null(group_by)) {
    if (!is.character(group_by) && length(group_by) != 1) {
      cli::cli_abort("{.val {group_by}} needs to contain characters with length of 1.")
    } else if (!column_exists(group_by, data)) {
      cli::cli_abort("The {.val {group_by}} column does not exist in the provided {.arg data}.")
    }
  }

  if (!c(p.adjust.method %in% p.adjust.methods))
    cli::cli_abort("Specified {.val {p.adjust.method}} is not valid. \nValid options: <{.val {p.adjust.methods}}>")

  ## MAIN
  #--------------------------------------------------------------------#

  result <- list()

  if (!is.null(group_by)) {
      pvalues_adjusted <- data[, {
        tmp <- rstatix::pairwise_wilcox_test(
          data   = .SD,
          formula = stats::reformulate(col_name, response = values),
          p.adjust.method = "fdr",
          paired = TRUE
        )
        tmp <- rstatix::add_significance(tmp)
        tmp <- rstatix::add_xy_position(tmp, x = group_by)
        tmp
      }, by = group_by]

      # Creates box_stats for half geom_box/geom_point
      setnames(data, old = group_by, new = "group_col")
      group_by <- "group_col"
      box_stats <- data[, .(
        ymin = min(get(values)),
        ymax = max(get(values)),
        lower = quantile(get(values), 0.25),
        middle = median(get(values)),
        upper = quantile(get(values), 0.75)
      ), by = .(group_numeric = as.numeric(as.factor(get(col_name))), group_col)]
    } else {
      pvalues_adjusted <- data %>%
        rstatix::pairwise_wilcox_test(formula = stats::reformulate(col_name, response = values),
                                      p.adjust.method = p.adjust.method,
                                      paired = paired) %>%
        rstatix::add_significance() %>%
        rstatix::add_xy_position(x = col_name)

      # Creates box_stats for half geom_box/geom_point
      box_stats <- data[, .(
        ymin = min(get(values)),
        ymax = max(get(values)),
        lower = quantile(get(values), 0.25),
        middle = median(get(values)),
        upper = quantile(get(values), 0.75)
      ), by = .(group_numeric = as.numeric(as.factor(get(col_name))))]
    }
    pvalues_adjusted.filtered <- pvalues_adjusted[grepl("\\*", pvalues_adjusted$p.adj.signif) ,]

    plt <- data %>%
      ggplot(mapping = aes(x = as.numeric(as.factor(.data[[col_name]])), y = .data[[values]]))
    # Custom half-boxplot using pre-computed stats
    if (!is.null(group_by)) {
      plt <- plt + geom_boxplot(
        data = box_stats,
        aes(x = group_numeric - 0.2,
            ymin = lower, ymax = upper,
            lower = lower, middle = middle, upper = upper,
            width = 0.4,
            group = interaction(group_numeric, group_col)),
        stat = "identity",
        fill = "white", color = "black",
        alpha = 0.8,
        inherit.aes = FALSE 
      )
    } else {
      plt <- plt + geom_boxplot(
        data = box_stats,
        aes(x = group_numeric - 0.2,
            ymin = lower, ymax = upper,
            lower = lower, middle = middle, upper = upper,
            width = 0.4,
            group = interaction(group_numeric)),
        stat = "identity",
        fill = "white", color = "black",
        alpha = 0.8,
        inherit.aes = FALSE 
      )
    }
    plt <- plt +
      # Points on right side
      geom_point(aes(color = as.factor(.data[[col_name]])), 
                position = position_nudge(x = 0.2), 
                shape = 20, size = 2) +
      geom_segment(
        data = box_stats,
        aes(x = group_numeric, y = ymin,
            xend = group_numeric, yend = ymax),
        color = "black", size = 0.3
      ) +
      # Top horizontal segment
      geom_segment(
        data = box_stats,
        aes(x = group_numeric - 0.1, y = ymax,
            xend = group_numeric, yend = ymax),
        color = "black", size = 0.3
      ) +
      # Bottom horizontal segment  
      geom_segment(
        data = box_stats,
        aes(x = group_numeric - 0.1, y = ymin,
            xend = group_numeric, yend = ymin),
        color = "black", size = 0.3
      )
    
    if (!is.null(group_by)) {
      plt <- plt +
        facet_wrap(~.data[[ group_by ]])
    }
      
    plt <- plt +
      theme_bw() +
      theme(legend.position = "none",
            text=element_text(size=14),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            axis.text = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.text.x = element_text(size=12)) +
      # Restore proper x-axis labels
      scale_x_continuous(breaks = seq_along(unique(data[[col_name]])),
                        labels = levels(as.factor(data[[col_name]]))) +
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
