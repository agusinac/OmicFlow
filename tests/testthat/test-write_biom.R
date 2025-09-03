test_that("Writes biom in expected format", {
    output_file <- paste0(tempdir(),"/test.biom")

    taxa <- metagenomics$new(
        biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
        metaData = "input/metagenomics/metadata.tsv",
        treeData = "input/metagenomics/rooted_tree.newick"
    )
    taxa$write_biom(filename = output_file)

    expect_true(file.exists(output_file))
    
    taxa1 <- metagenomics$new(
        biomData = output_file,
        metaData = "input/metagenomics/metadata.tsv",
        treeData = "input/metagenomics/rooted_tree.newick"
    )
    
    expect_equal(taxa$print(), taxa1$print())
    file.remove(output_file)
})
