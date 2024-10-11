# fold_plot <- function(dt, X, Y, method, title = NULL, taxa_labels = FALSE, pvalues = NULL, pvalue.col = "pvalue", pvalue.threshold = 0.01) {
#     plt <- dt %>% 
#       ggplot(mapping = aes(x = .data[[ X ]],
#                            y = .data[[ Y ]])) +
#       geom_boxplot()
#   }
#   
#   plt <- plt +
#     theme_bw() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
#           axis.text.y = element_text(size=12),
#           axis.text = element_text(size=12),
#           text = element_text(size=12),
#           legend.text = element_text(size=12),
#           legend.title = element_text(size=14),
#           legend.position = "none",
#           axis.title.y = element_blank(),
#           strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"))
#   plt <- plt +
#     scale_fill_gradient2(name = paste0("log2( A / B )"),
#                          low = "blue",
#                          mid = "white",
#                          high = "red",
#                          na.value = "grey80") +
#     scale_y_discrete(limits = rev(levels(as.factor(dt[[ Y ]])))) +
#     labs(x = paste(X),
#          y = paste(Y)) +
#     ggtitle(title)
#   
#   return(plt)
# }
