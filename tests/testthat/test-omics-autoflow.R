## Load example data
taxa <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`omics$autoFlow()` -- Argument checks", {

})

test_that("`omics$autoFlow()` -- Behavioral checks", { 
    taxa <- metagenomics$new(
        biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
        metaData = "input/metagenomics/metadata.tsv",
        treeData = "input/metagenomics/rooted_tree.newick"
    )
  
    suppressWarnings(
        taxa_autoflow <- taxa$autoFlow(
            feature_contrast = c("Phylum", "Family", "Genus"),
            pvalue.threshold = 1,
            report = FALSE
        )
    )
    expect_true(length(taxa_autoflow$plots) == length(taxa_autoflow$data))

    prot <- proteomics$new(
        metaData = "input/proteomics/metadata.csv",
        countData = "input/proteomics/counts.csv",
        treeData = "input/proteomics/tree.newick"
    )
    prot$scale(method = "clr")

    suppressWarnings(
        prot_autoflow <- prot$autoFlow(
            pvalue.threshold = 1,
            distance_metrics = "euclidean",
            report = FALSE,
        )
    )
    expect_true(length(prot_autoflow$plots) == length(prot_autoflow$data))
})