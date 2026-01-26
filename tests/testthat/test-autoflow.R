test_that("Testing autoFlow on metagenomics data", {
    taxa <- metagenomics$new(
        biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
        metaData = "input/metagenomics/metadata.tsv",
        treeData = "input/metagenomics/rooted_tree.newick"
    )
  
    suppressWarnings(
        taxa_autoflow <- taxa$autoFlow(
            feature_contrast = c("Phylum", "Family", "Genus"),
            report = FALSE
        )
    )
    expect_true(length(taxa_autoflow$plots) == length(taxa_autoflow$data))
})

test_that("Testing autoFlow on proteomics data", {
    prot <- proteomics$new(
        metaData = "input/proteomics/metadata.csv",
        countData = "input/proteomics/counts.csv",
        treeData = "input/proteomics/tree.newick"
    )

    suppressWarnings(
        prot_autoflow <- prot$autoFlow(
            report = FALSE
        )
    )
    expect_true(length(prot_autoflow$plots) == length(prot_autoflow$data))
})