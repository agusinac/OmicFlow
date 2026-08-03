# `omics$distance()` -- Argument checks

    Code
      test$distance()
    Condition
      Error in `test$distance()`:
      ! "metric" must be specified!

---

    Code
      test$distance(metric = "nothing")
    Condition
      Error in `test$distance()`:
      ! "nothing" is not a valid metric.  Valid options: "bray", "jaccard", "cosine", "manhattan", "jsd", "canberra", "unifrac", "euclidean", and "aitchison"

---

    Code
      test$distance(metric = c("n1", "n2"))
    Condition
      Error in `test$distance()`:
      ! "n1" and "n2" needs to be a character with a length of 1

---

    Code
      test$distance(metric = 1)
    Condition
      Error in `test$distance()`:
      ! 1 needs to be a character with a length of 1

---

    Code
      test$distance(metric = "unifrac")
    Condition
      Error in `test$distance()`:
      ! The specified "unifrac" is invalid since no treeData is supplied.

---

    Code
      test$distance(metric = "bray", threads = "1")
    Condition
      Error in `test$distance()`:
      ! "1" need to be a whole number!

---

    Code
      test$distance(metric = "bray", threads = 50.2)
    Condition
      Error in `test$distance()`:
      ! 50.2 need to be a whole number!

# `omics$ordination()` -- Argument checks

    Code
      test$ordination()
    Condition
      Error in `test$ordination()`:
      ! "group_by" must be specified!

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", metric = "nothing")
    Condition
      Error in `self$distance()`:
      ! "nothing" is not a valid metric.  Valid options: "bray", "jaccard", "cosine", "manhattan", "jsd", "canberra", "unifrac", "euclidean", and "aitchison"

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", metric = c("n1", "n2"))
    Condition
      Error in `self$distance()`:
      ! "n1" and "n2" needs to be a character with a length of 1

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", metric = 1)
    Condition
      Error in `self$distance()`:
      ! 1 needs to be a character with a length of 1

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", method = "nothing")
    Condition
      Error in `test$ordination()`:
      ! "nothing" is not a valid method.  Valid options: "pcoa" and "nmds"

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", method = c("n1", "n2"))
    Condition
      Error in `test$ordination()`:
      ! "n1" and "n2" needs to be a character with a length of 1

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", method = 1)
    Condition
      Error in `test$ordination()`:
      ! 1 needs to be a character with a length of 1

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", perm = "999")
    Condition
      Error in `test$ordination()`:
      ! Permutations "999" need to be a whole number.

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", perm = 50.2)
    Condition
      Error in `test$ordination()`:
      ! Permutations 50.2 need to be a whole number.

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", threads = "999")
    Condition
      Error in `self$distance()`:
      ! "999" need to be a whole number!

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", threads = 50.2)
    Condition
      Error in `self$distance()`:
      ! 50.2 need to be a whole number!

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", perm_design = list())
    Condition
      Error in `test$ordination()`:
      ! perm_design must be a function.

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", perm_design = function(x)
        print(x))
    Output
          SAMPLE_ID SAMPLEPAIR_ID CONTRAST_treatment
             <char>        <char>             <char>
       1: Sample_01          S001              tumor
       2: Sample_02          S002              tumor
       3: Sample_03          S003              tumor
       4: Sample_04          S004              tumor
       5: Sample_05          S005              tumor
       6: Sample_06          S001            control
       7: Sample_07          S002            control
       8: Sample_08          S010            control
       9: Sample_09          S004            control
      10: Sample_10          S012            control
    Condition
      Error in `Math.data.frame()`:
      ! non-numeric-alike variable(s) in data frame: SAMPLE_ID, SAMPLEPAIR_ID, CONTRAST_treatment

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", distmat = list())
    Condition
      Error in `test$ordination()`:
      ! `distmat` need to be <Matrix> or <dist>

---

    Code
      test$ordination(group_by = "CONTRAST_treatment", distmat = dist(c(2, 1, 2)))
    Condition
      Error in `test$ordination()`:
      ! None "SAMPLE_ID" from metaData match the "distmat" colnames!

# `omics$ordination()` -- Behavioral checks

    Code
      res$anova_data
    Output
                   pairs Df   SumsOfSqs  F.Model        R2 p.value p.adj
      1 tumor vs control  1 0.006240019 1.747221 0.1792532   0.048 0.048

---

    Code
      res$dist
    Output
                 Sample_01  Sample_02  Sample_03  Sample_04  Sample_05  Sample_06
      Sample_01 0.00000000 0.07933219 0.08058596 0.08757372 0.07874846 0.08086030
      Sample_02 0.07933219 0.00000000 0.08314432 0.07240245 0.09406856 0.09353530
      Sample_03 0.08058596 0.08314432 0.00000000 0.06990762 0.08091202 0.10921049
      Sample_04 0.08757372 0.07240245 0.06990762 0.00000000 0.08728288 0.10084309
      Sample_05 0.07874846 0.09406856 0.08091202 0.08728288 0.00000000 0.09045771
      Sample_06 0.08086030 0.09353530 0.10921049 0.10084309 0.09045771 0.00000000
      Sample_07 0.08037143 0.07754409 0.09187468 0.07365468 0.07144072 0.07715973
      Sample_08 0.07297827 0.09569631 0.11082168 0.10370835 0.08862186 0.09298599
      Sample_09 0.08594647 0.09028899 0.10424058 0.09718785 0.08314188 0.08158129
      Sample_10 0.09088034 0.08573821 0.09778241 0.08786596 0.08495964 0.09235266
                 Sample_07  Sample_08  Sample_09  Sample_10
      Sample_01 0.08037143 0.07297827 0.08594647 0.09088034
      Sample_02 0.07754409 0.09569631 0.09028899 0.08573821
      Sample_03 0.09187468 0.11082168 0.10424058 0.09778241
      Sample_04 0.07365468 0.10370835 0.09718785 0.08786596
      Sample_05 0.07144072 0.08862186 0.08314188 0.08495964
      Sample_06 0.07715973 0.09298599 0.08158129 0.09235266
      Sample_07 0.00000000 0.08935864 0.08680636 0.08128468
      Sample_08 0.08935864 0.00000000 0.08246781 0.09139688
      Sample_09 0.08680636 0.08246781 0.00000000 0.09526405
      Sample_10 0.08128468 0.09139688 0.09526405 0.00000000

---

    Code
      res$pcs
    Output
                    PC1          PC2           PC3          PC4          PC5
                  <num>        <num>         <num>        <num>        <num>
       1:  0.0072101776  0.029714822  0.0018471905 -0.014536154  0.027838251
       2: -0.0214624306 -0.013484928  0.0154181397 -0.036493483 -0.006591635
       3: -0.0522260108  0.026312754 -0.0112301580  0.003636017  0.003308725
       4: -0.0426991088 -0.010705457  0.0008322959 -0.008072932 -0.010296509
       5: -0.0004576682  0.015508271 -0.0182265575  0.041803074 -0.003298443
       6:  0.0354108168 -0.025917453 -0.0304518736 -0.012168930  0.023749743
       7: -0.0007701352 -0.023972168 -0.0112475934  0.008664586  0.008737484
       8:  0.0414329816  0.027906898  0.0330893060 -0.002498912  0.001390973
       9:  0.0329538797  0.005532739 -0.0171634014 -0.008808913 -0.045606729
      10:  0.0006074980 -0.030895477  0.0371326517  0.028475647  0.000768141
                    PC6          PC7           PC8          PC9  groups samples
                  <num>        <num>         <num>        <num>  <char>  <char>
       1: -0.0064199160  0.003073528  0.0230138376  0.008760282   tumor       1
       2: -0.0004950177  0.028766704 -0.0092615818  0.004685848   tumor       2
       3: -0.0181172652 -0.004637543 -0.0071902419 -0.018412224   tumor       3
       4:  0.0169078288 -0.028136707  0.0008658901  0.014970275   tumor       4
       5:  0.0062675321  0.017090263 -0.0085259781  0.015236377   tumor       5
       6: -0.0146750367 -0.010227517 -0.0128787822  0.001741224 control       6
       7:  0.0339371353  0.008005643  0.0111758238 -0.016653213 control       7
       8:  0.0173604831 -0.010061199 -0.0142200671 -0.005611360 control       8
       9: -0.0104523469 -0.001638505  0.0103980004 -0.003620876 control       9
      10: -0.0243133968 -0.002234667  0.0066230993 -0.001096332 control      10

---

    Code
      res$anova_data
    Output
                   pairs Df   SumsOfSqs  F.Model        R2 p.value p.adj
      1 tumor vs control  1 0.006298002 1.738349 0.1785056   0.059 0.059

---

    Code
      res$dist
    Output
                 Sample_01  Sample_02  Sample_03  Sample_04  Sample_05  Sample_06
      Sample_01 0.00000000 0.07944700 0.08274841 0.08826924 0.07964405 0.08158715
      Sample_02 0.07944700 0.00000000 0.08259256 0.07386369 0.09422101 0.09355350
      Sample_03 0.08274841 0.08259256 0.00000000 0.07116931 0.08104109 0.11042979
      Sample_04 0.08826924 0.07386369 0.07116931 0.00000000 0.08887421 0.10275058
      Sample_05 0.07964405 0.09422101 0.08104109 0.08887421 0.00000000 0.09122483
      Sample_06 0.08158715 0.09355350 0.11042979 0.10275058 0.09122483 0.00000000
      Sample_07 0.08135964 0.07809416 0.09338074 0.07473912 0.07332313 0.07786200
      Sample_08 0.07280218 0.09450806 0.11102348 0.10357251 0.08871931 0.09329472
      Sample_09 0.08657947 0.09081649 0.10582457 0.09784647 0.08366459 0.08298722
      Sample_10 0.09075501 0.08573181 0.09807939 0.08854323 0.08542345 0.09351215
                 Sample_07  Sample_08  Sample_09  Sample_10
      Sample_01 0.08135964 0.07280218 0.08657947 0.09075501
      Sample_02 0.07809416 0.09450806 0.09081649 0.08573181
      Sample_03 0.09338074 0.11102348 0.10582457 0.09807939
      Sample_04 0.07473912 0.10357251 0.09784647 0.08854323
      Sample_05 0.07332313 0.08871931 0.08366459 0.08542345
      Sample_06 0.07786200 0.09329472 0.08298722 0.09351215
      Sample_07 0.00000000 0.08955225 0.08723115 0.08213706
      Sample_08 0.08955225 0.00000000 0.08210798 0.09106892
      Sample_09 0.08723115 0.08210798 0.00000000 0.09555412
      Sample_10 0.08213706 0.09106892 0.09555412 0.00000000

---

    Code
      res$pcs
    Output
                    PC1          PC2          PC3          PC4           PC5
                  <num>        <num>        <num>        <num>         <num>
       1:  0.0076978094  0.027992736  0.004740304  0.016720513  0.0280936316
       2: -0.0206953966 -0.013054214  0.024573418  0.029512417 -0.0066984615
       3: -0.0535756290  0.024538327 -0.013587871  0.001324876  0.0051634716
       4: -0.0428225504 -0.011508857  0.004143655  0.007441152 -0.0136845873
       5: -0.0008548591  0.017115708 -0.030671411 -0.034483555 -0.0004535603
       6:  0.0365697077 -0.028449141 -0.024916498  0.019309100  0.0227011914
       7:  0.0001346456 -0.026431913 -0.012299065 -0.004687886  0.0054253771
       8:  0.0399135733  0.029795010  0.031099414 -0.004163203  0.0013869169
       9:  0.0335146390  0.006224459 -0.013029208  0.009178771 -0.0462909592
      10:  0.0001180601 -0.026222114  0.029947263 -0.040152186  0.0043569797
                   PC6          PC7           PC8          PC9  groups samples
                 <num>        <num>         <num>        <num>  <char>  <char>
       1:  0.003224642  0.001583015  0.0257856528  0.004482240   tumor       1
       2:  0.005652085 -0.030264302 -0.0062459844  0.007157178   tumor       2
       3:  0.017132410  0.004060149 -0.0102393400 -0.017536602   tumor       3
       4: -0.018302051  0.027474825  0.0009535253  0.014304402   tumor       4
       5: -0.004958133 -0.016814909 -0.0043372577  0.016651787   tumor       5
       6:  0.014295543  0.010753442 -0.0131264996  0.003448149 control       6
       7: -0.034375840 -0.010272532  0.0082198403 -0.017306069 control       7
       8: -0.018465642  0.006629405 -0.0172652852 -0.003614373 control       8
       9:  0.013076141  0.002512118  0.0097888770 -0.005322126 control       9
      10:  0.022720846  0.004338788  0.0064664715 -0.002264587 control      10

