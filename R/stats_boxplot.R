#' Create a boxplot with jittered points and statistical tests
#' 
#' @description
#' Visualizes group differences using a half-boxplot, half-jitter plot with automatic statistical testing and p-value annotation. 
#' Thus far only the wilcoxon rank test (non-paired) and wilcoxon signed rank test (paired) are supported.
#' This function is built into the class \link{omics} with method \code{alpha_diversity()}.
#' @param data A \link[base]{data.frame} or \link[data.table]{data.table} table.
#' @param values A column name of a continuous variable.
#' @param col_name A column name of a categorical variable.
#' @param group_by A column name to perform grouped statistical test (default: \code{NULL}).
#' @param palette An object with names and hexcode or color names, see \link{colormap}.
#' @param method A character variable indicating what method is used to compute the diversity (default: \code{"custom"}).
#' @param paired A boolean value to perform paired analysis with \link[matrixTests]{row_wilcoxon_paired} (default: \code{FALSE}).
#' @param p.adjust.method A character variable to specify the p.adjust.method to be used (default: \code{"fdr"}).
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#' 
#' @examples  
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
#'   "shannon" = div,
#'   "treatment" = c(rep("healthy", n_col / 2), rep("tumor", n_col / 2)),
#'   "sex" = c(rep("male", n_col / 4), rep("female", n_col / 4))
#' )
#' 
#' colors <- OmicFlow::colormap(dt, "treatment")
#' 
#' # Comparing two groups
#' plt <- OmicFlow::boxjitter_test(
#'  data = dt,
#'  values = "shannon",
#'  col_name = "treatment",
#'  palette = colors,
#'  method = "shannon",
#'  paired = FALSE,
#'  p.adjust.method = "fdr"
#' )
#' 
#' # Performing a test while stratifying the plot in two groups
#' plt <- OmicFlow::boxjitter_test(
#'  data = dt,
#'  values = "shannon",
#'  col_name = "treatment",
#'  group_by = "sex",
#'  palette = colors,
#'  method = "shannon",
#'  paired = FALSE,
#'  p.adjust.method = "fdr"
#' )
#' @export

boxjitter_test <- function(
  data,
  values,
  col_name,
  group_by = NULL,
  palette,
  method = "custom",
  paired = FALSE,
  p.adjust.method = "fdr"
  ) {
  
  ## Error handling
  #--------------------------------------------------------------------#
  
  if (!inherits(data, "data.frame") && !inherits(data, "data.table"))
    cli::cli_abort("{.val data} must be a {.cls data.frame} or {.cls data.table}.")

  if (!is.character(values) || length(values) != 1) {
    cli::cli_abort("{.val values} needs to contain characters with length of 1.")
  } else if (!column_exists(values, data)) {
    cli::cli_abort("The {.val values} column does not exist in the provided {.arg data}.")
  }
  
  if (!is.character(col_name) || length(col_name) != 1) {
    cli::cli_abort("{.val col_name} needs to contain characters with length of 1.")
  } else if (!column_exists(col_name, data)) {
    cli::cli_abort("The {.val col_name} column does not exist in the provided {.arg data}.")
  }
  
  if (!is.null(group_by)) {
    if (!is.character(group_by) || length(group_by) != 1) {
      cli::cli_abort("{.val group_by} needs to contain characters with length of 1.")
    } else if (!column_exists(group_by, data)) {
      cli::cli_abort("The {.val group_by} column does not exist in the provided {.arg data}.")
    }
  }

  if (!is.character(palette)) {
    cli::cli_abort("{.val palette} needs to contain characters.")
  } else if (!is.color(palette)) {
    cli::cli_abort("{.val palette} contains invalid colors.")
  }
    
  if (!is.character(method)) {
    cli::cli_abort("{.val method} needs to be a character {.cls vector}.")
  }

  if (!is.logical(paired))
    cli::cli_abort("{.val paired} needs to be either `TRUE` or `FALSE`.")

  if (!c(p.adjust.method %in% stats::p.adjust.methods))
    cli::cli_abort("{.val {p.adjust.method}} is not a valid option. \nValid options: {.val {p.adjust.methods}}")
  
  ## MAIN
  #--------------------------------------------------------------------#

  data_tmp <- data.table::copy(data)
  result <- list()

  if (!is.null(group_by)) {
    data.table::setnames(data_tmp, old = group_by, new = "group_by")
    group_by <- "group_by"

    pvalues_adjusted <- data_tmp[, {
      tmp <- pairwise_wilcox_test(
          data = .SD,
          x_col = values,
          g_col = col_name,
          p.adjust.method = p.adjust.method,
          paired = paired
      )
      tmp
    }, by = group_by]

    # Creates box_stats for half geom_box
    box_stats <- data_tmp[, .(
      ymin = base::min(base::get(values)),
      ymax = base::max(base::get(values)),
      lower = stats::quantile(base::get(values), 0.25),
      middle = stats::median(base::get(values)),
      upper = stats::quantile(base::get(values), 0.75)
    ), by = .(group_numeric = as.numeric(as.factor(base::get(col_name))), group_by)]
  } else {
    pvalues_adjusted <- data_tmp[, {
      tmp <- pairwise_wilcox_test(
        data = .SD,
        x_col = values,
        g_col = col_name,
        p.adjust.method = p.adjust.method,
        paired = paired
      )
      tmp
    }]

    # Creates box_stats for half geom_box
    box_stats <- data_tmp[, .(
      ymin = base::min(base::get(values)),
      ymax = base::max(base::get(values)),
      lower = stats::quantile(base::get(values), 0.25),
      middle = stats::median(base::get(values)),
      upper = stats::quantile(base::get(values), 0.75)
    ), by = .(group_numeric = as.numeric(as.factor(get(col_name))))]
  }

  ## Check if any significant pairs
  signif_rows <- pvalues_adjusted$p.adj < 0.05
  xmin <- xmax <- y.position <- p.adj <- NULL
  if (any(signif_rows)) {
    pvalues_adjusted.filtered <- pvalues_adjusted[signif_rows, ]
    pvalues_dt <- data.table::as.data.table(pvalues_adjusted.filtered)
    pvalues_dt[, `:=` (
      xmid = (xmin + xmax) / 2,
      bracket_height = y.position - 0.02 * diff(range(y.position, na.rm = TRUE)),
      label = round(p.adj, 2)
      )
    ]
  }

  plt <- ggplot2::ggplot(
    data = data_tmp,
    mapping = ggplot2::aes(
      x = as.numeric(as.factor(.data[[col_name]])),
      y = .data[[values]]
    )
  )
  # Custom half-boxplot using pre-computed stats
  if (!is.null(group_by)) {
    suppressWarnings(
      plt <- plt + ggplot2::geom_boxplot(
        data = box_stats,
        mapping = ggplot2::aes(x = .data$group_numeric - 0.2,
            ymin = .data$lower, ymax = .data$upper,
            lower = .data$lower, middle = .data$middle, upper = .data$upper,
            width = 0.4,
            group = base::interaction(.data$group_numeric, .data$group_by)),
        stat = "identity",
        fill = "white", color = "black",
        alpha = 0.8,
        inherit.aes = FALSE 
      )
    )
  } else {
    suppressWarnings(
      plt <- plt + ggplot2::geom_boxplot(
        data = box_stats,
        mapping = ggplot2::aes(x = .data$group_numeric - 0.2,
            ymin = .data$lower, ymax = .data$upper,
            lower = .data$lower, middle = .data$middle, upper = .data$upper,
            width = 0.4,
            group = base::interaction(.data$group_numeric)),
        stat = "identity",
        fill = "white", color = "black",
        alpha = 0.8,
        inherit.aes = FALSE 
      )
    )
  }

  plt <- plt +
    # Points on right side
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        x = as.numeric(as.factor(.data[[col_name]])) + 0.2, 
        color = as.factor(.data[[col_name]])
      ), 
      position = ggplot2::position_jitter(width = 0.1, height = 0, seed = 1970), 
      shape = 20, size = 2
    ) +
    ggplot2::geom_segment(
      data = box_stats,
      mapping = ggplot2::aes(
        x = .data$group_numeric, 
        y = .data$ymin,
        xend = .data$group_numeric, 
        yend = .data$ymax
      ),
      color = "black", 
      linewidth = 0.3
    ) +
    # Top horizontal segment
    ggplot2::geom_segment(
      data = box_stats,
      mapping = ggplot2::aes(
        x = .data$group_numeric - 0.1, 
        y = .data$ymax,
        xend = .data$group_numeric, 
        yend = .data$ymax
      ),
      color = "black", 
      linewidth = 0.3
    ) +
    # Bottom horizontal segment  
    ggplot2::geom_segment(
      data = box_stats,
      mapping = ggplot2::aes(
        x = .data$group_numeric - 0.1, 
        y = .data$ymin,
        xend = .data$group_numeric, 
        yend = .data$ymin
      ),
      color = "black", 
      linewidth = 0.3
    )

  if (!is.null(group_by)) {
    plt <- plt +
      ggplot2::facet_wrap(~.data[[ group_by ]])
  }

  plt <- plt + 
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      text = ggplot2::element_text(size=14),
      legend.text = ggplot2::element_text(size=12),
      legend.title = ggplot2::element_text(size=14),
      axis.text = ggplot2::element_text(size=12),
      axis.text.y = ggplot2::element_text(size=12),
      axis.text.x = ggplot2::element_text(size=12)
    ) +
    # Restore proper x-axis labels
    ggplot2::scale_x_continuous(
      breaks = seq_along(unique(data_tmp[[col_name]])),
      labels = levels(as.factor(data_tmp[[col_name]]))
    ) +
    ggplot2::scale_colour_manual(
      name = col_name,
      values = palette
    )

  ## Adding significant bars
  if (any(signif_rows)) {
    plt <- plt + 
    # horizontal line
    ggplot2::geom_segment(
      data = pvalues_dt,
      mapping = ggplot2::aes(
        x = .data$xmin, xend = .data$xmax,
        y = .data$y.position, yend = .data$y.position
      ),
      inherit.aes = FALSE
    ) +
    # Left vertical line
    ggplot2::geom_segment(
      data = pvalues_dt,
      mapping = ggplot2::aes(
        x = .data$xmin, xend = .data$xmin,
        y = .data$y.position, yend = .data$y.position - 0.003
      ),
      inherit.aes = FALSE
    ) +
    # Right vertical line
    ggplot2::geom_segment(
      data = pvalues_dt,
      mapping = ggplot2::aes(
        x = .data$xmax, xend = .data$xmax,
        y = .data$y.position, yend = .data$y.position - 0.003
      ),
      inherit.aes = FALSE
    ) +
    # Label
    ggplot2::geom_text(
      data = pvalues_dt,
      mapping = ggplot2::aes(
        x = .data$xmid,
        y = .data$y.position,
        label = .data$label
      ),
      vjust = -0.4,
      inherit.aes = FALSE
    )
  }

  plt <- plt +
    ggplot2::labs(
      title = NULL,
      subtitle = paste0(
      "Attribute: ", col_name,
      ", test: ", ifelse(paired, "Wilcox signed rank test", "Mann-Whitney U test"),
      ", p.adjusted by ", p.adjust.method),
      x = "sample groups",
      y = paste0("Alpha diversity metric: ", method)
    )

  result <- list(
    plot = plt,
    stats = pvalues_adjusted
  )
  
  return(result)
}

#' Diversity plot
#'
#' @description Creates an Alpha diversity plot. This function is built into the class \link{omics} with method \code{alpha_diversity()}.
#' It computes the pairwise wilcox test, paired or non-paired, given a data frame and adds useful labelling.
#' 
#' This function has been changed throughout iterations and become more generic than it's original use. Therefore the function name is changed to `boxjitter_test`, the current `diversity_plot` will be deprecated within the near future.
#'
#' @param ... arguments passed to \link{boxjitter_test}
#' @export
diversity_plot <- function(...) {
  lifecycle::deprecate_warn(
    when = "1.6.0",
    what = "diversity_plot()",
    with = "boxjitter_plot()",
    details = "This function has been generalized for all group comparisons, not just alpha diversity."
  )
  boxjitter_test(...)
}
