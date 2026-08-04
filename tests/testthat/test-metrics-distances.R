## Load example data
test <- metagenomics$new(
  biomData = "input/metagenomics/biom_with_taxonomy_hdf5.biom",
  metaData = "input/metagenomics/metadata.tsv",
  treeData = "input/metagenomics/rooted_tree.newick"
)

test_that("`bray()` -- Argument checks", {
    expect_snapshot(bray(x = c(1, 2, 3)), error = TRUE)
    expect_snapshot(bray(x = data.frame()), error = TRUE)
    expect_snapshot(bray(x = test$countData, weighted = "FALSE"), error = TRUE)
    expect_snapshot(bray(x = test$countData, threads = "1"), error = TRUE)
    expect_snapshot(bray(x = test$countData, threads = 1.9), error = TRUE)
    expect_snapshot(bray(x = test$countData, threads = c(1, 2)), error = TRUE)
})

test_that("`bray()` -- Behavioral checks", { 
    expect_snapshot(bray(test$countData))
    expect_snapshot(bray(test$countData, weighted = FALSE))
    expect_snapshot(bray(as.matrix(test$countData)))
})


test_that("`jaccard()` -- Argument checks", {
    expect_snapshot(jaccard(x = c(1, 2, 3)), error = TRUE)
    expect_snapshot(jaccard(x = data.frame()), error = TRUE)
    expect_snapshot(jaccard(x = test$countData, weighted = "FALSE"), error = TRUE)
    expect_snapshot(jaccard(x = test$countData, threads = "1"), error = TRUE)
    expect_snapshot(jaccard(x = test$countData, threads = 1.9), error = TRUE)
    expect_snapshot(jaccard(x = test$countData, threads = c(1, 2)), error = TRUE)
})

test_that("`jaccard()` -- Behavioral checks", {
    expect_snapshot(jaccard(test$countData))
    expect_snapshot(jaccard(test$countData, weighted = FALSE))
    expect_snapshot(jaccard(as.matrix(test$countData)))
})

test_that("`cosine()` -- Argument checks", {
    expect_snapshot(cosine(x = c(1, 2, 3)), error = TRUE)
    expect_snapshot(cosine(x = data.frame()), error = TRUE)
    expect_snapshot(cosine(x = test$countData, weighted = "FALSE"), error = TRUE)
    expect_snapshot(cosine(x = test$countData, threads = "1"), error = TRUE)
    expect_snapshot(cosine(x = test$countData, threads = 1.9), error = TRUE)
    expect_snapshot(cosine(x = test$countData, threads = c(1, 2)), error = TRUE)
})

test_that("`cosine()` -- Behavioral checks", { 
    expect_snapshot(cosine(test$countData))
    expect_snapshot(cosine(test$countData, weighted = FALSE))
    expect_snapshot(cosine(as.matrix(test$countData)))
})

test_that("`manhattan()` -- Argument checks", {
    expect_snapshot(manhattan(x = c(1, 2, 3)), error = TRUE)
    expect_snapshot(manhattan(x = data.frame()), error = TRUE)
    expect_snapshot(manhattan(x = test$countData, weighted = "FALSE"), error = TRUE)
    expect_snapshot(manhattan(x = test$countData, threads = "1"), error = TRUE)
    expect_snapshot(manhattan(x = test$countData, threads = 1.9), error = TRUE)
    expect_snapshot(manhattan(x = test$countData, threads = c(1, 2)), error = TRUE)
})

test_that("`manhattan()` -- Behavioral checks", { 
    expect_snapshot(manhattan(test$countData))
    expect_snapshot(manhattan(test$countData, weighted = FALSE))
    expect_snapshot(manhattan(as.matrix(test$countData)))
})

test_that("`jsd()` -- Argument checks", {
    expect_snapshot(jsd(x = c(1, 2, 3)), error = TRUE)
    expect_snapshot(jsd(x = data.frame()), error = TRUE)
    expect_snapshot(jsd(x = test$countData, weighted = "FALSE"), error = TRUE)
    expect_snapshot(jsd(x = test$countData, threads = "1"), error = TRUE)
    expect_snapshot(jsd(x = test$countData, threads = 1.9), error = TRUE)
    expect_snapshot(jsd(x = test$countData, threads = c(1, 2)), error = TRUE)

    test$countData[1,1] <- -1
    expect_snapshot(jsd(x = test$countData), error = TRUE)
    test$reset()
})

test_that("`jsd()` -- Behavioral checks", { 
    expect_snapshot(jsd(test$countData))
    expect_snapshot(jsd(test$countData, weighted = FALSE))
    expect_snapshot(jsd(as.matrix(test$countData)))
})

test_that("`canberra()` -- Argument checks", {
    expect_snapshot(canberra(x = c(1, 2, 3)), error = TRUE)
    expect_snapshot(canberra(x = data.frame()), error = TRUE)
    expect_snapshot(canberra(x = test$countData, weighted = "FALSE"), error = TRUE)
    expect_snapshot(canberra(x = test$countData, threads = "1"), error = TRUE)
    expect_snapshot(canberra(x = test$countData, threads = 1.9), error = TRUE)
    expect_snapshot(canberra(x = test$countData, threads = c(1, 2)), error = TRUE)
})

test_that("`canberra()` -- Behavioral checks", { 
    expect_snapshot(canberra(test$countData))
    expect_snapshot(canberra(test$countData, weighted = FALSE))
    expect_snapshot(canberra(as.matrix(test$countData)))
})

test_that("`unifrac()` -- Argument checks", {
    expect_snapshot(unifrac(tree = data.frame()), error = TRUE)
    expect_snapshot(unifrac(tree = test$treeData, x = c(1, 2, 3)), error = TRUE)
    expect_snapshot(unifrac(tree = test$treeData, x = data.frame()), error = TRUE)
    expect_snapshot(unifrac(tree = test$treeData, x = test$countData, weighted = "FALSE"), error = TRUE)
    expect_snapshot(unifrac(tree = test$treeData, x = test$countData, normalize = "FALSE"), error = TRUE)
    expect_snapshot(unifrac(tree = test$treeData, x = test$countData, threads = "1"), error = TRUE)
    expect_snapshot(unifrac(tree = test$treeData, x = test$countData, threads = 1.9), error = TRUE)
    expect_snapshot(unifrac(tree = test$treeData, x = test$countData, threads = c(1, 2)), error = TRUE)

    test$countData[1,1] <- -1
    expect_snapshot(unifrac(tree = test$treeData, x = test$countData), error = TRUE)
    test$reset()
})

test_that("`unifrac()` -- Behavioral checks", {
    ## Weighted Normalised UniFrac
    expect_snapshot(unifrac(x = test$countData, tree = test$treeData, weighted = TRUE, normalize = TRUE))
    expect_snapshot(unifrac(x = as.matrix(test$countData), tree = test$treeData, weighted = TRUE, normalize = TRUE))

    ## Weighted UniFrac
    expect_snapshot(unifrac(x = test$countData, tree = test$treeData, weighted = TRUE, normalize = FALSE))
    expect_snapshot(unifrac(x = as.matrix(test$countData), tree = test$treeData, weighted = TRUE, normalize = FALSE))

    ## Unweighted UniFrac
    expect_snapshot(unifrac(x = test$countData, tree = test$treeData, weighted = FALSE, normalize = FALSE))
    expect_snapshot(unifrac(x = as.matrix(test$countData), tree = test$treeData, weighted = FALSE, normalize = FALSE))
})

test_that("`euclidean()` -- Argument checks", {
    expect_snapshot(euclidean(x = c(1, 2, 3)), error = TRUE)
    expect_snapshot(euclidean(x = data.frame()), error = TRUE)
    expect_snapshot(euclidean(x = test$countData, weighted = "FALSE"), error = TRUE)
    expect_snapshot(euclidean(x = test$countData, threads = "1"), error = TRUE)
    expect_snapshot(euclidean(x = test$countData, threads = 1.9), error = TRUE)
    expect_snapshot(euclidean(x = test$countData, threads = c(1, 2)), error = TRUE)
})

test_that("`euclidean()` -- Behavioral checks", { 
    expect_snapshot(euclidean(test$countData))
    expect_snapshot(euclidean(test$countData, weighted = FALSE))
    expect_snapshot(euclidean(as.matrix(test$countData)))
})