#' Diversity plot with ggplot2
#'
#' @description Creates an Alpha diversity plot. This function is built into the \code{alpha_diversity} method from the abstract class \link[OmicFlow]{tools} and inherited by other omics classes, such as;
#' \link[OmicFlow]{metataxonomics}, transcriptomics, metabolomics and proteomics.
#'
#' @param dt A \code{data.frame} or \code{data.table} computed from \link[vegan]{diversity}.
#' @param values A column name of a continuous variable.
#' @param col_name A column name of a categorical variable.
#' @param palette An object with names and hexcode or color names, see \link{fetch_palette} or \link[stats]{setNames}.
#' @param method A character variable indicating what method is used to compute the diversity.
#' @param p.adjust.method A character variable to specify the p.adjust.method to be used, default is 'fdr'.
#' @return A \link[ggplot2]{ggplot2} object to be further modified
#'
#' @export

diversity_plot <- function(dt, values, col_name, palette, method, paired = FALSE, p.adjust.method = "fdr") {
  pvalues_adjusted <- dt %>%
    rstatix::pairwise_wilcox_test(formula = stats::reformulate(col_name, response = values),
                                  p.adjust.method = p.adjust.method,
                                  paired = paired) %>%
    rstatix::add_significance() %>%
    rstatix::add_xy_position(x = col_name)

  pvalues_adjusted.filtered <- pvalues_adjusted[grepl("\\*", pvalues_adjusted$p.adj.signif) ,]

  return(dt %>%
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
                subtitle = paste0("Attribute: ", col_name, ", test: Mann-whitney, p.adjusted by ", p.adjust.method),
                x = "sample groups",
                y = paste0("Alpha diversity: ", method))
  )
}
