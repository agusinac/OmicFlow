# test_that("Testing autoFlow while ignoring dynamic date", {
#   taxa <- metagenomics$new(
#     biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
#     metaData = "input/metagenomics/metadata.tsv",
#     treeData = "input/metagenomics/rooted_tree.newick"
#   )
  
#   taxa$autoFlow(filename = paste0(getwd(), "/report.html"))
  
#   # Define a transform function to replace dates (e.g., YYYY-MM-DD) with a placeholder
#   transform_fn <- function(x) {
#     gsub("\\d{4}-\\d{2}-\\d{2}", "<DATE>", x)
#   }
  
#   expect_snapshot_file("report.html", transform = transform_fn)
#   file.remove("report.html")
# })