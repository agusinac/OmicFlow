test_that("Testing dissimilarity metrics", {
    taxa <- metagenomics$new(
    biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
    metaData = "input/metagenomics/metadata.tsv",
    treeData = "input/metagenomics/rooted_tree.newick"
    )
    taxa$normalize()

    ## Testing Weighted Normalized UniFrac
    wunifrac_n <- unifrac(
        x = taxa$countData,
        tree = taxa$treeData,
        weighted = TRUE,
        normalized = TRUE
    )
    expect_snapshot(cat(wunifrac_n))

    ## Testing Weighted UniFrac
    wunifrac <- unifrac(
        x = taxa$countData,
        tree = taxa$treeData,
        weighted = TRUE,
        normalized = FALSE
    )
    expect_snapshot(cat(wunifrac))

    ## Testing Unweighted UniFrac
    uunifrac <- unifrac(
        x = taxa$countData,
        tree = taxa$treeData,
        weighted = FALSE,
        normalized = FALSE
    )
    expect_snapshot(cat(uunifrac))

    ## Testing Cosine
    cos <- cosine(
        x = taxa$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(cos))

    ## Testing Jaccard
    jac <- jaccard(
        x = taxa$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(jac))

    ## Testing Bray
    br <- bray(
        x = taxa$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(br))

    ## Testing Canberra
    can <- canberra(
        x = taxa$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(can))

    ## Testing Jensen-Shannon Divergence
    jsd_res <- jsd(
        x = taxa$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(jsd_res))

    ## Testing Manhattan
    man <- manhattan(
        x = taxa$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(man))
    }
)