taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`write_biom()` -- Argument check", {

})

test_that("`write_biom()` -- Behavior check", {
    output_file <- paste0(tempdir(),"/test.biom")

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
