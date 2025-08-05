test_that("Tests autoFlow report creation", {
    output_file <- "report.html"
    output_file <- normalizePath(output_file)

    taxa <- metagenomics$new(
        biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
        metaData = "input/metagenomics/metadata.tsv",
        treeData = "input/metagenomics/rooted_tree.newick"
    )
    suppressWarnings(taxa$autoFlow(filename = output_file))

    expect_true(file.exists(output_file))
    
    file.remove(output_file)
})