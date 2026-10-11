# `omics$sample_subset()` -- Argument checks

    Code
      test$sample_subset(CONTRAST_treatment == "nothing")
    Condition
      Error in `test$sample_subset()`:
      ! The expression resulted in all rows in metaData to be `FALSE`.

# `omics$sample_subset()` -- Behavioral checks

    Code
      test
    Message
      
      -- <proteomics> object 
      metaData: 3 variables x 5 samples
      countData: 5 samples x 50 features
      featureData: 0 attributes x 50 features

---

    Code
      test
    Message
      
      -- <proteomics> object 
      metaData: 3 variables x 10 samples
      countData: 10 samples x 50 features
      featureData: 0 attributes x 50 features

---

    Code
      test$metaData
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

---

    Code
      test$countData
    Output
      50 x 10 sparse Matrix of class "dgCMatrix"
    Message
        [[ suppressing 10 column names 'Sample_01', 'Sample_02', 'Sample_03' ... ]]
    Output
                                                                        
      P00001 24.62 24.96 24.91 27.58 25.54 24.53 32.74 20.31 27.67 23.98
      P00002 23.87 22.98 25.84 33.90 34.34 34.35 29.80 22.56 24.17 27.77
      P00003 28.28 23.54 20.62 22.08 33.70 25.85 34.26 29.60 25.41 31.07
      P00004 20.85 24.12 25.42 22.54 32.35 25.58 29.26 22.47 26.56 31.69
      P00005 27.03 28.87 28.56 29.10 24.79 32.65 27.39 25.32 32.05 21.72
      P00006 27.26 23.80 30.27 32.24 33.17 32.71 34.64 22.80 27.81 29.19
      P00007 32.19 21.85 34.57 32.65 32.01 24.79 27.36 33.46 30.44 33.66
      P00008 25.55 23.45 30.53 31.82 29.17 21.99 29.83 23.57 32.72 29.51
      P00009 28.20 28.96 20.17 20.29 21.09 29.26 28.98 34.77 32.69 24.12
      P00010 22.55 23.17 28.03 30.47 26.32 31.87 34.21 20.32 25.88 25.29
      P00011 29.37 26.96 32.55 32.23 25.17 25.07 25.52 21.61 22.30 30.06
      P00012 33.23 29.71 32.10 28.53 31.28 33.58 33.17 23.67 29.59 33.60
      P00013 24.21 34.41 21.21 27.20 23.28 22.96 26.82 30.84 24.31 30.95
      P00014 25.98 30.15 23.58 22.42 24.38 31.91 27.45 20.49 33.93 33.80
      P00015 31.44 26.68 34.49 21.34 25.34 31.32 26.91 28.24 22.33 25.74
      P00016 30.04 25.37 20.57 22.44 29.48 33.67 29.27 30.24 34.46 32.81
      P00017 23.07 26.84 33.75 20.40 33.39 24.84 29.06 24.50 20.02 20.82
      P00018 25.36 26.68 30.89 30.64 31.16 21.29 31.79 25.83 30.61 23.81
      P00019 25.39 23.68 23.01 31.42 26.86 33.67 28.33 30.99 29.45 25.51
      P00020 30.35 30.42 32.60 32.86 20.54 34.33 31.53 34.44 31.60 25.41
      P00021 28.04 26.18 25.95 26.56 28.53 30.07 26.06 31.41 33.39 24.21
      P00022 30.66 24.92 25.89 26.26 26.59 31.17 27.66 28.73 27.67 25.83
      P00023 28.08 28.59 27.09 28.78 29.01 26.60 27.85 26.94 31.24 25.22
      P00024 31.23 34.50 28.75 32.37 34.27 21.73 34.89 25.34 33.90 24.18
      P00025 26.30 29.93 25.29 31.90 24.04 30.13 26.44 25.76 21.38 32.02
      P00026 22.57 29.37 20.43 24.89 29.88 30.97 34.94 23.09 27.45 34.93
      P00027 31.55 32.85 34.93 34.37 21.15 27.26 31.80 22.08 22.96 31.30
      P00028 33.23 31.62 34.36 29.80 21.07 22.57 27.74 25.84 34.91 29.05
      P00029 28.24 32.51 28.27 26.88 25.55 30.14 27.54 23.99 20.55 28.91
      P00030 24.17 21.37 21.53 29.09 24.46 23.94 33.67 30.55 23.26 28.73
      P00031 27.32 26.89 23.57 24.33 28.27 25.12 23.97 26.12 33.36 23.58
      P00032 33.93 28.99 32.90 30.77 25.55 23.15 22.61 23.98 21.78 26.81
      P00033 25.23 33.80 31.07 33.85 32.69 20.24 26.00 26.01 25.54 34.93
      P00034 34.31 34.74 27.46 30.12 29.31 25.65 28.08 22.96 22.44 21.24
      P00035 30.43 20.57 28.70 22.79 25.99 28.43 23.67 32.45 22.50 20.78
      P00036 33.34 28.67 20.24 25.22 24.49 30.20 25.64 27.92 34.59 27.41
      P00037 22.71 31.00 27.08 21.85 25.72 31.18 28.70 25.93 32.16 20.81
      P00038 29.44 23.73 20.64 21.62 30.52 34.25 23.13 28.61 34.30 34.26
      P00039 34.84 24.51 26.95 24.47 34.20 22.45 32.02 34.56 26.00 29.50
      P00040 21.95 31.00 29.45 32.57 31.59 24.87 29.58 29.75 32.78 24.05
      P00041 24.96 33.60 30.10 34.89 23.29 21.99 31.04 25.00 20.97 34.95
      P00042 32.98 23.15 21.31 26.56 30.74 29.56 26.58 24.97 24.54 25.38
      P00043 31.66 25.37 22.15 23.04 29.96 24.96 28.70 34.14 27.36 25.58
      P00044 32.41 26.72 33.63 34.46 30.59 29.74 23.84 25.72 29.69 24.95
      P00045 29.05 33.60 21.84 29.91 24.19 24.54 26.96 28.46 27.66 21.57
      P00046 27.37 25.84 30.93 24.47 30.68 21.07 22.53 27.67 31.96 28.70
      P00047 31.71 27.76 34.26 21.79 29.91 29.94 29.21 22.08 28.48 25.00
      P00048 33.26 21.88 20.65 29.00 20.62 31.39 34.37 23.60 25.30 20.06
      P00049 23.12 20.45 20.29 21.80 20.92 28.30 27.17 30.77 29.89 34.92
      P00050 24.61 31.58 22.99 31.83 24.20 28.09 31.27 24.46 23.59 25.54

---

    Code
      test$treeData
    Output
      NULL

---

    Code
      test
    Message
      
      -- <proteomics> object 
      metaData: 3 variables x 10 samples
      countData: 10 samples x 50 features
      featureData: 0 attributes x 50 features

# `omics$feature_subset()` -- Argument checks

    Code
      test$feature_subset(grepl("nothing", FEATURE_ID))
    Condition
      Error in `test$feature_subset()`:
      ! The expression resulted in all rows in featureData to be `FALSE`.

# `omics$feature_subset()` -- Behavioral checks

    Code
      test
    Message
      
      -- <proteomics> object 
      metaData: 3 variables x 10 samples
      countData: 10 samples x 26 features
      featureData: 0 attributes x 26 features

---

    Code
      test
    Message
      
      -- <proteomics> object 
      metaData: 3 variables x 10 samples
      countData: 10 samples x 50 features
      featureData: 0 attributes x 50 features

# `omics$samplepair_subset()` -- Argument checks

    Code
      test$samplepair_subset()
    Condition
      Error in `test$samplepair_subset()`:
      ! "SAMPLEPAIR_ID" doesn't exist in the metaData.

---

    Code
      test$samplepair_subset(num_unique_pairs = 2)
    Condition
      Error in `test$samplepair_subset()`:
      ! "SAMPLEPAIR_ID" doesn't exist in the metaData.

# `omics$samplepair_subset()` -- Behavioral checks

    Code
      test
    Message
      
      -- <proteomics> object 
      metaData: 3 variables x 6 samples
      countData: 6 samples x 50 features
      featureData: 0 attributes x 50 features

