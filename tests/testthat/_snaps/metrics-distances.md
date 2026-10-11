# `bray()` -- Argument checks

    Code
      bray(x = c(1, 2, 3))
    Condition
      Error in `bray()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      bray(x = data.frame())
    Condition
      Error in `bray()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      bray(x = test$countData, weighted = "FALSE")
    Condition
      Error in `bray()`:
      ! "weighted" needs to be either `TRUE` or `FALSE`.

---

    Code
      bray(x = test$countData, threads = "1")
    Condition
      Error in `bray()`:
      ! "threads" must be a whole number.

---

    Code
      bray(x = test$countData, threads = 1.9)
    Condition
      Error in `bray()`:
      ! "threads" must be a whole number.

---

    Code
      bray(x = test$countData, threads = c(1, 2))
    Condition
      Error in `bray()`:
      ! "threads" must be a single whole number.

# `bray()` -- Behavioral checks

    Code
      bray(test$countData)
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.8845188 1.0000000          
      S120 1.0000000 0.9470058 1.0000000

---

    Code
      bray(test$countData, weighted = FALSE)
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.9344262 1.0000000          
      S120 1.0000000 0.9682540 1.0000000

---

    Code
      bray(as.matrix(test$countData))
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.8845188 1.0000000          
      S120 1.0000000 0.9470058 1.0000000

# `jaccard()` -- Argument checks

    Code
      jaccard(x = c(1, 2, 3))
    Condition
      Error in `jaccard()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      jaccard(x = data.frame())
    Condition
      Error in `jaccard()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      jaccard(x = test$countData, weighted = "FALSE")
    Condition
      Error in `jaccard()`:
      ! "weighted" needs to be either `TRUE` or `FALSE`.

---

    Code
      jaccard(x = test$countData, threads = "1")
    Condition
      Error in `jaccard()`:
      ! "threads" must be a whole number.

---

    Code
      jaccard(x = test$countData, threads = 1.9)
    Condition
      Error in `jaccard()`:
      ! "threads" must be a whole number.

---

    Code
      jaccard(x = test$countData, threads = c(1, 2))
    Condition
      Error in `jaccard()`:
      ! "threads" must be a single whole number.

# `jaccard()` -- Behavioral checks

    Code
      jaccard(test$countData)
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.9387211 1.0000000          
      S120 1.0000000 0.9727817 1.0000000

---

    Code
      jaccard(test$countData, weighted = FALSE)
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.9661017 1.0000000          
      S120 1.0000000 0.9838710 1.0000000

---

    Code
      jaccard(as.matrix(test$countData))
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.9387211 1.0000000          
      S120 1.0000000 0.9727817 1.0000000

# `cosine()` -- Argument checks

    Code
      cosine(x = c(1, 2, 3))
    Condition
      Error in `cosine()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      cosine(x = data.frame())
    Condition
      Error in `cosine()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      cosine(x = test$countData, weighted = "FALSE")
    Condition
      Error in `cosine()`:
      ! "weighted" needs to be either `TRUE` or `FALSE`.

---

    Code
      cosine(x = test$countData, threads = "1")
    Condition
      Error in `cosine()`:
      ! "threads" must be a whole number.

---

    Code
      cosine(x = test$countData, threads = 1.9)
    Condition
      Error in `cosine()`:
      ! "threads" must be a whole number.

---

    Code
      cosine(x = test$countData, threads = c(1, 2))
    Condition
      Error in `cosine()`:
      ! "threads" must be a single whole number.

# `cosine()` -- Behavioral checks

    Code
      cosine(test$countData)
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.7913906 1.0000000          
      S120 1.0000000 0.9282048 1.0000000

---

    Code
      cosine(test$countData, weighted = FALSE)
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.9333333 1.0000000          
      S120 1.0000000 0.9682380 1.0000000

---

    Code
      cosine(as.matrix(test$countData))
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.7913906 1.0000000          
      S120 1.0000000 0.9282048 1.0000000

# `manhattan()` -- Argument checks

    Code
      manhattan(x = c(1, 2, 3))
    Condition
      Error in `manhattan()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      manhattan(x = data.frame())
    Condition
      Error in `manhattan()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      manhattan(x = test$countData, weighted = "FALSE")
    Condition
      Error in `manhattan()`:
      ! "weighted" needs to be either `TRUE` or `FALSE`.

---

    Code
      manhattan(x = test$countData, threads = "1")
    Condition
      Error in `manhattan()`:
      ! "threads" must be a whole number.

---

    Code
      manhattan(x = test$countData, threads = 1.9)
    Condition
      Error in `manhattan()`:
      ! "threads" must be a whole number.

---

    Code
      manhattan(x = test$countData, threads = c(1, 2))
    Condition
      Error in `manhattan()`:
      ! "threads" must be a single whole number.

# `manhattan()` -- Behavioral checks

    Code
      manhattan(test$countData)
    Output
           S100 S103 S115
      S103 2058          
      S115 2114 2272     
      S120 2005 1787 2219

---

    Code
      manhattan(test$countData, weighted = FALSE)
    Output
           S100 S103 S115
      S103  115          
      S115  114  137     
      S120  111  122  133

---

    Code
      manhattan(as.matrix(test$countData))
    Output
           S100 S103 S115
      S103 2058          
      S115 2114 2272     
      S120 2005 1787 2219

# `jsd()` -- Argument checks

    Code
      jsd(x = c(1, 2, 3))
    Condition
      Error in `jsd()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      jsd(x = data.frame())
    Condition
      Error in `jsd()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      jsd(x = test$countData, weighted = "FALSE")
    Condition
      Error in `jsd()`:
      ! "weighted" needs to be either `TRUE` or `FALSE`.

---

    Code
      jsd(x = test$countData, threads = "1")
    Condition
      Error in `jsd()`:
      ! "threads" must be a whole number.

---

    Code
      jsd(x = test$countData, threads = 1.9)
    Condition
      Error in `jsd()`:
      ! "threads" must be a whole number.

---

    Code
      jsd(x = test$countData, threads = c(1, 2))
    Condition
      Error in `jsd()`:
      ! "threads" must be a single whole number.

---

    Code
      jsd(x = test$countData)
    Condition
      Error in `jsd()`:
      ! "x" must be non-negative.

# `jsd()` -- Behavioral checks

    Code
      jsd(test$countData)
    Output
               S100     S103     S115
      S103 1029.000                  
      S115 1002.005 1136.000         
      S120 1002.500  889.269 1109.500

---

    Code
      jsd(test$countData, weighted = FALSE)
    Output
           S100 S103 S115
      S103 57.5          
      S115 57.0 68.5     
      S120 55.5 61.0 66.5

---

    Code
      jsd(as.matrix(test$countData))
    Output
               S100     S103     S115
      S103 1029.000                  
      S115 1002.005 1136.000         
      S120 1002.500  889.269 1109.500

# `canberra()` -- Argument checks

    Code
      canberra(x = c(1, 2, 3))
    Condition
      Error in `canberra()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      canberra(x = data.frame())
    Condition
      Error in `canberra()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      canberra(x = test$countData, weighted = "FALSE")
    Condition
      Error in `canberra()`:
      ! "weighted" needs to be either `TRUE` or `FALSE`.

---

    Code
      canberra(x = test$countData, threads = "1")
    Condition
      Error in `canberra()`:
      ! "threads" must be a whole number.

---

    Code
      canberra(x = test$countData, threads = 1.9)
    Condition
      Error in `canberra()`:
      ! "threads" must be a whole number.

---

    Code
      canberra(x = test$countData, threads = c(1, 2))
    Condition
      Error in `canberra()`:
      ! "threads" must be a single whole number.

# `canberra()` -- Behavioral checks

    Code
      canberra(test$countData)
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.9760365 1.0000000          
      S120 1.0000000 0.9852530 1.0000000

---

    Code
      canberra(test$countData, weighted = FALSE)
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.9661017 1.0000000          
      S120 1.0000000 0.9838710 1.0000000

---

    Code
      canberra(as.matrix(test$countData))
    Output
                S100      S103      S115
      S103 1.0000000                    
      S115 0.9760365 1.0000000          
      S120 1.0000000 0.9852530 1.0000000

# `unifrac()` -- Argument checks

    Code
      unifrac(tree = data.frame())
    Condition
      Error in `unifrac()`:
      ! "tree" must be a <phylo>.

---

    Code
      unifrac(tree = test$treeData, x = c(1, 2, 3))
    Condition
      Error in `unifrac()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      unifrac(tree = test$treeData, x = data.frame())
    Condition
      Error in `unifrac()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      unifrac(tree = test$treeData, x = test$countData, weighted = "FALSE")
    Condition
      Error in `unifrac()`:
      ! "weighted" needs to be either `TRUE` or `FALSE`.

---

    Code
      unifrac(tree = test$treeData, x = test$countData, normalize = "FALSE")
    Condition
      Error in `unifrac()`:
      ! "normalize" needs to be either `TRUE` or `FALSE`.

---

    Code
      unifrac(tree = test$treeData, x = test$countData, threads = "1")
    Condition
      Error in `unifrac()`:
      ! "threads" must be a whole number.

---

    Code
      unifrac(tree = test$treeData, x = test$countData, threads = 1.9)
    Condition
      Error in `unifrac()`:
      ! "threads" must be a whole number.

---

    Code
      unifrac(tree = test$treeData, x = test$countData, threads = c(1, 2))
    Condition
      Error in `unifrac()`:
      ! "threads" must be a single whole number.

---

    Code
      unifrac(tree = test$treeData, x = test$countData)
    Condition
      Error in `unifrac()`:
      ! "x" must be non-negative.

# `unifrac()` -- Behavioral checks

    Code
      unifrac(x = test$countData, tree = test$treeData, weighted = TRUE, normalize = TRUE)
    Output
                S100      S103      S115
      S103 0.3722812                    
      S115 0.1069508 0.3299658          
      S120 0.3381820 0.1645516 0.3038010

---

    Code
      unifrac(x = as.matrix(test$countData), tree = test$treeData, weighted = TRUE,
      normalize = TRUE)
    Output
                S100      S103      S115
      S103 0.3722812                    
      S115 0.1069508 0.3299658          
      S120 0.3381820 0.1645516 0.3038010

---

    Code
      unifrac(x = test$countData, tree = test$treeData, weighted = TRUE, normalize = FALSE)
    Output
                S100      S103      S115
      S103 0.5576522                    
      S115 0.1570318 0.5058064          
      S120 0.4956745 0.2518210 0.4559069

---

    Code
      unifrac(x = as.matrix(test$countData), tree = test$treeData, weighted = TRUE,
      normalize = FALSE)
    Output
                S100      S103      S115
      S103 0.5576522                    
      S115 0.1570318 0.5058064          
      S120 0.4956745 0.2518210 0.4559069

---

    Code
      unifrac(x = test$countData, tree = test$treeData, weighted = FALSE, normalize = FALSE)
    Output
                S100      S103      S115
      S103 0.8498807                    
      S115 0.7214772 0.8362140          
      S120 0.8465006 0.6551804 0.8027923

---

    Code
      unifrac(x = as.matrix(test$countData), tree = test$treeData, weighted = FALSE,
      normalize = FALSE)
    Output
                S100      S103      S115
      S103 0.8498807                    
      S115 0.7214772 0.8362140          
      S120 0.8465006 0.6551804 0.8027923

# `euclidean()` -- Argument checks

    Code
      euclidean(x = c(1, 2, 3))
    Condition
      Error in `euclidean()`:
      ! "x" must be a <matrix>, <denseMatrix> or <sparseMatrix>, not a <vector>.

---

    Code
      euclidean(x = data.frame())
    Condition
      Error in `euclidean()`:
      ! "x" isn't a <matrix>, <denseMatrix> or <sparseMatrix>.

---

    Code
      euclidean(x = test$countData, weighted = "FALSE")
    Condition
      Error in `euclidean()`:
      ! "weighted" needs to be either `TRUE` or `FALSE`.

---

    Code
      euclidean(x = test$countData, threads = "1")
    Condition
      Error in `euclidean()`:
      ! "threads" must be a whole number.

---

    Code
      euclidean(x = test$countData, threads = 1.9)
    Condition
      Error in `euclidean()`:
      ! "threads" must be a whole number.

---

    Code
      euclidean(x = test$countData, threads = c(1, 2))
    Condition
      Error in `euclidean()`:
      ! "threads" must be a single whole number.

# `euclidean()` -- Behavioral checks

    Code
      euclidean(test$countData)
    Output
               S100     S103     S115
      S103 280.9199                  
      S115 300.8322 299.7099         
      S120 283.0530 227.4357 301.7101

---

    Code
      euclidean(test$countData, weighted = FALSE)
    Output
               S100     S103     S115
      S103 10.72381                  
      S115 10.67708 11.70470         
      S120 10.53565 11.04536 11.53256

---

    Code
      euclidean(as.matrix(test$countData))
    Output
               S100     S103     S115
      S103 280.9199                  
      S115 300.8322 299.7099         
      S120 283.0530 227.4357 301.7101

