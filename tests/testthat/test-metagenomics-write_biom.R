## Load example data
taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`write_biom()` -- Argument check", {
    expect_error(taxa$write_biom(), 'argument "filename" is missing, with no default')
    expect_snapshot(taxa$write_biom(filename = list()), error = TRUE)
    expect_snapshot(taxa$write_biom(filename = c("file1.biom", "file2.biom")), error = TRUE)

    output_file <- paste0(tempdir(), "/test.biom")
    data.table::fwrite(x = data.table::data.table(A = c(2, 3, 1)), file = output_file)
    expect_snapshot(taxa$write_biom(filename = output_file), error = TRUE)
    file.remove(output_file)
})

test_that("`write_biom()` -- Behavior check", {
    output_file <- paste0(tempdir(),"/test_new.biom")

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
