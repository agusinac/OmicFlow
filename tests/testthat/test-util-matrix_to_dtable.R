## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`matrix_to_dtable()` -- Argument checks", {

})

test_that("`matrix_to_dtable()` -- Behavioral checks", { 
 
})

test_that("Testing exported utils functions", {
    # Load test data
    taxa <- metagenomics$new(
        biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
        metaData = "input/metagenomics/metadata.tsv",
        treeData = "input/metagenomics/rooted_tree.newick"
    )

    # Testing `matrix_to_dtable`
    expect_no_error(matrix_to_dtable(taxa$countData))
    expect_error(matrix_to_dtable(dt_sparse))
})