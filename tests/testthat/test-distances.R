test_that("Testing dissimilarity metrics on sparse data", {
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

test_that("Testing dissimilarity metrics on dense data", {
    prot_file <- proteomics$new(
        metaData = "input/proteomics/metadata.csv",
        countData = "input/proteomics/counts.csv",
        treeData = "input/proteomics/tree.newick"
    )
    prot$normalize()

    ## Testing Weighted Normalized UniFrac
    wunifrac_n <- unifrac(
        x = prot$countData,
        tree = prot$treeData,
        weighted = TRUE,
        normalized = TRUE
    )
    expect_snapshot(cat(wunifrac_n))

    ## Testing Weighted UniFrac
    wunifrac <- unifrac(
        x = prot$countData,
        tree = prot$treeData,
        weighted = TRUE,
        normalized = FALSE
    )
    expect_snapshot(cat(wunifrac))

    ## Testing Unweighted UniFrac
    uunifrac <- unifrac(
        x = prot$countData,
        tree = prot$treeData,
        weighted = FALSE,
        normalized = FALSE
    )
    expect_snapshot(cat(uunifrac))

    ## Testing Cosine
    cos <- cosine(
        x = prot$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(cos))

    ## Testing Jaccard
    jac <- jaccard(
        x = prot$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(jac))

    ## Testing Bray
    br <- bray(
        x = prot$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(br))

    ## Testing Canberra
    can <- canberra(
        x = prot$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(can))

    ## Testing Jensen-Shannon Divergence
    jsd_res <- jsd(
        x = prot$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(jsd_res))

    ## Testing Manhattan
    man <- manhattan(
        x = prot$countData,
        weighted = TRUE
    )
    expect_snapshot(cat(man))
    }
)