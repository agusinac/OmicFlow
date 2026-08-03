# `omics$composition()` -- Argument checks

    Code
      taxa$composition(feature_rank = "Genus", col_name = "1")
    Condition
      Error in `taxa$composition()`:
      ! The specified "1" does not exist in the metaData.

---

    Code
      taxa$composition(feature_rank = "Genus", col_name = c("treatment", "sex"))
    Condition
      Error in `taxa$composition()`:
      ! "treatment" and "sex" must be a character and of length 1.

---

    Code
      taxa$composition(feature_rank = "Genus", col_name = "nonexisting")
    Condition
      Error in `taxa$composition()`:
      ! The specified "nonexisting" does not exist in the metaData.

---

    Code
      taxa$composition(feature_rank = "Genus", feature_top = "10")
    Condition
      Error in `taxa$composition()`:
      ! "10" must be a whole number.

---

    Code
      taxa$composition(feature_rank = "Genus", feature_top = c(10, 15))
    Condition
      Error in `taxa$composition()`:
      ! "10 and 15 must be a single element."

---

    Code
      taxa$composition(feature_rank = "Genus", feature_top = 16)
    Condition
      Error in `taxa$composition()`:
      ! The "feature_top" cannot be higher than 15. This may lead that colors are difficult to be distinguished for color-blind people, therefore the limit is set to 15.

# `omics$composition()` -- Behavioral checks

    Code
      res$data
    Output
      Key: <SAMPLE_ID>
          SAMPLE_ID                   Genus      value CONTRAST_sex
             <char>                  <fctr>      <num>       <char>
       1:      S100         Woesearchaeales 0.02022059         male
       2:      S100            Sideroxydans 0.19209559         male
       3:      S100              uncultured 0.00000000         male
       4:      S100           Acinetobacter 0.16727941         male
       5:      S100  Candidatus_Omnitrophus 0.05147059         male
       6:      S100             Pseudomonas 0.09099265         male
       7:      S100          Nitrosarchaeum 0.06341912         male
       8:      S100 Candidatus_Peribacteria 0.00000000         male
       9:      S100              Nitrospira 0.03033088         male
      10:      S100           Sulfurifustis 0.00000000         male
      11:      S100                   Other 0.12867647         male
      12:      S103         Woesearchaeales 0.28865979       female
      13:      S103            Sideroxydans 0.00000000       female
      14:      S103              uncultured 0.11752577       female
      15:      S103           Acinetobacter 0.02268041       female
      16:      S103  Candidatus_Omnitrophus 0.06391753       female
      17:      S103             Pseudomonas 0.00000000       female
      18:      S103          Nitrosarchaeum 0.00000000       female
      19:      S103 Candidatus_Peribacteria 0.04123711       female
      20:      S103              Nitrospira 0.00000000       female
      21:      S103           Sulfurifustis 0.03402062       female
      22:      S103                   Other 0.20618557       female
      23:      S115         Woesearchaeales 0.03149002       female
      24:      S115            Sideroxydans 0.20890937       female
      25:      S115              uncultured 0.04531490       female
      26:      S115           Acinetobacter 0.09677419       female
      27:      S115  Candidatus_Omnitrophus 0.08986175       female
      28:      S115             Pseudomonas 0.04685100       female
      29:      S115          Nitrosarchaeum 0.09831029       female
      30:      S115 Candidatus_Peribacteria 0.00000000       female
      31:      S115              Nitrospira 0.03917051       female
      32:      S115           Sulfurifustis 0.01228879       female
      33:      S115                   Other 0.17204301       female
      34:      S120         Woesearchaeales 0.12540894         male
      35:      S120            Sideroxydans 0.00000000         male
      36:      S120              uncultured 0.15158124         male
      37:      S120           Acinetobacter 0.00000000         male
      38:      S120  Candidatus_Omnitrophus 0.06434024         male
      39:      S120             Pseudomonas 0.05779716         male
      40:      S120          Nitrosarchaeum 0.00000000         male
      41:      S120 Candidatus_Peribacteria 0.05125409         male
      42:      S120              Nitrospira 0.01853871         male
      43:      S120           Sulfurifustis 0.03707743         male
      44:      S120                   Other 0.20828790         male
          SAMPLE_ID                   Genus      value CONTRAST_sex
             <char>                  <fctr>      <num>       <char>

---

    Code
      res$palette
    Output
              Woesearchaeales            Sideroxydans              uncultured 
                    "#A50026"               "#D73027"               "#F46D43" 
                Acinetobacter  Candidatus_Omnitrophus             Pseudomonas 
                    "#FDAE61"               "#FEE090"               "#E0F3F8" 
               Nitrosarchaeum Candidatus_Peribacteria              Nitrospira 
                    "#ABD9E9"               "#74ADD1"               "#4575B4" 
                Sulfurifustis                   Other 
                    "#313695"             "lightgrey" 

---

    Code
      res$data
    Output
      Key: <SAMPLE_ID>
          SAMPLE_ID                      Genus       value CONTRAST_sex
             <char>                     <fctr>       <num>       <char>
       1:      S100            Woesearchaeales 0.020220588         male
       2:      S100               Sideroxydans 0.192095588         male
       3:      S100                 uncultured 0.000000000         male
       4:      S100              Acinetobacter 0.167279412         male
       5:      S100     Candidatus_Omnitrophus 0.051470588         male
       6:      S100                Pseudomonas 0.090992647         male
       7:      S100             Nitrosarchaeum 0.063419118         male
       8:      S100    Candidatus_Peribacteria 0.000000000         male
       9:      S100                 Nitrospira 0.030330882         male
      10:      S100              Sulfurifustis 0.000000000         male
      11:      S100              Aquabacterium 0.034007353         male
      12:      S100                      WOR-1 0.001838235         male
      13:      S100  Candidatus_Falkowbacteria 0.000000000         male
      14:      S100 Candidatus_Kerfeldbacteria 0.000000000         male
      15:      S100              Parcubacteria 0.000000000         male
      16:      S100                      Other 0.092830882         male
      17:      S103            Woesearchaeales 0.288659794       female
      18:      S103               Sideroxydans 0.000000000       female
      19:      S103                 uncultured 0.117525773       female
      20:      S103              Acinetobacter 0.022680412       female
      21:      S103     Candidatus_Omnitrophus 0.063917526       female
      22:      S103                Pseudomonas 0.000000000       female
      23:      S103             Nitrosarchaeum 0.000000000       female
      24:      S103    Candidatus_Peribacteria 0.041237113       female
      25:      S103                 Nitrospira 0.000000000       female
      26:      S103              Sulfurifustis 0.034020619       female
      27:      S103              Aquabacterium 0.000000000       female
      28:      S103                      WOR-1 0.026804124       female
      29:      S103  Candidatus_Falkowbacteria 0.036082474       female
      30:      S103 Candidatus_Kerfeldbacteria 0.032989691       female
      31:      S103              Parcubacteria 0.017525773       female
      32:      S103                      Other 0.092783505       female
      33:      S115            Woesearchaeales 0.031490015       female
      34:      S115               Sideroxydans 0.208909370       female
      35:      S115                 uncultured 0.045314900       female
      36:      S115              Acinetobacter 0.096774194       female
      37:      S115     Candidatus_Omnitrophus 0.089861751       female
      38:      S115                Pseudomonas 0.046850998       female
      39:      S115             Nitrosarchaeum 0.098310292       female
      40:      S115    Candidatus_Peribacteria 0.000000000       female
      41:      S115                 Nitrospira 0.039170507       female
      42:      S115              Sulfurifustis 0.012288786       female
      43:      S115              Aquabacterium 0.025345622       female
      44:      S115                      WOR-1 0.000000000       female
      45:      S115  Candidatus_Falkowbacteria 0.000000000       female
      46:      S115 Candidatus_Kerfeldbacteria 0.000000000       female
      47:      S115              Parcubacteria 0.001536098       female
      48:      S115                      Other 0.145161290       female
      49:      S120            Woesearchaeales 0.125408942         male
      50:      S120               Sideroxydans 0.000000000         male
      51:      S120                 uncultured 0.151581243         male
      52:      S120              Acinetobacter 0.000000000         male
      53:      S120     Candidatus_Omnitrophus 0.064340240         male
      54:      S120                Pseudomonas 0.057797165         male
      55:      S120             Nitrosarchaeum 0.000000000         male
      56:      S120    Candidatus_Peribacteria 0.051254089         male
      57:      S120                 Nitrospira 0.018538713         male
      58:      S120              Sulfurifustis 0.037077426         male
      59:      S120              Aquabacterium 0.000000000         male
      60:      S120                      WOR-1 0.028353326         male
      61:      S120  Candidatus_Falkowbacteria 0.020719738         male
      62:      S120 Candidatus_Kerfeldbacteria 0.022900763         male
      63:      S120              Parcubacteria 0.027262814         male
      64:      S120                      Other 0.109051254         male
          SAMPLE_ID                      Genus       value CONTRAST_sex
             <char>                     <fctr>       <num>       <char>

---

    Code
      res$palette
    Output
                 Woesearchaeales               Sideroxydans 
                       "#000000"                  "#004949" 
                      uncultured              Acinetobacter 
                       "#009292"                  "#ff6db6" 
          Candidatus_Omnitrophus                Pseudomonas 
                       "#ffb6db"                  "#490092" 
                  Nitrosarchaeum    Candidatus_Peribacteria 
                       "#006ddb"                  "#b66dff" 
                      Nitrospira              Sulfurifustis 
                       "#6db6ff"                  "#b6dbff" 
                   Aquabacterium                      WOR-1 
                       "#920000"                  "#924900" 
       Candidatus_Falkowbacteria Candidatus_Kerfeldbacteria 
                       "#db6d00"                  "#24ff24" 
                   Parcubacteria                      Other 
                       "#ffff6d"                "lightgrey" 

---

    Code
      res$data
    Output
                            Genus SAMPLE_ID      value
                           <fctr>    <char>      <num>
       1:         Woesearchaeales      S100 0.02022059
       2:            Sideroxydans      S100 0.19209559
       3:              uncultured      S100 0.00000000
       4:           Acinetobacter      S100 0.16727941
       5:  Candidatus_Omnitrophus      S100 0.05147059
       6:             Pseudomonas      S100 0.09099265
       7:          Nitrosarchaeum      S100 0.06341912
       8: Candidatus_Peribacteria      S100 0.00000000
       9:              Nitrospira      S100 0.03033088
      10:           Sulfurifustis      S100 0.00000000
      11:                   Other      S100 0.12867647
      12:         Woesearchaeales      S103 0.28865979
      13:            Sideroxydans      S103 0.00000000
      14:              uncultured      S103 0.11752577
      15:           Acinetobacter      S103 0.02268041
      16:  Candidatus_Omnitrophus      S103 0.06391753
      17:             Pseudomonas      S103 0.00000000
      18:          Nitrosarchaeum      S103 0.00000000
      19: Candidatus_Peribacteria      S103 0.04123711
      20:              Nitrospira      S103 0.00000000
      21:           Sulfurifustis      S103 0.03402062
      22:                   Other      S103 0.20618557
      23:         Woesearchaeales      S115 0.03149002
      24:            Sideroxydans      S115 0.20890937
      25:              uncultured      S115 0.04531490
      26:           Acinetobacter      S115 0.09677419
      27:  Candidatus_Omnitrophus      S115 0.08986175
      28:             Pseudomonas      S115 0.04685100
      29:          Nitrosarchaeum      S115 0.09831029
      30: Candidatus_Peribacteria      S115 0.00000000
      31:              Nitrospira      S115 0.03917051
      32:           Sulfurifustis      S115 0.01228879
      33:                   Other      S115 0.17204301
      34:         Woesearchaeales      S120 0.12540894
      35:            Sideroxydans      S120 0.00000000
      36:              uncultured      S120 0.15158124
      37:           Acinetobacter      S120 0.00000000
      38:  Candidatus_Omnitrophus      S120 0.06434024
      39:             Pseudomonas      S120 0.05779716
      40:          Nitrosarchaeum      S120 0.00000000
      41: Candidatus_Peribacteria      S120 0.05125409
      42:              Nitrospira      S120 0.01853871
      43:           Sulfurifustis      S120 0.03707743
      44:                   Other      S120 0.20828790
                            Genus SAMPLE_ID      value
                           <fctr>    <char>      <num>

---

    Code
      res$palette
    Output
              Woesearchaeales            Sideroxydans              uncultured 
                    "#A50026"               "#D73027"               "#F46D43" 
                Acinetobacter  Candidatus_Omnitrophus             Pseudomonas 
                    "#FDAE61"               "#FEE090"               "#E0F3F8" 
               Nitrosarchaeum Candidatus_Peribacteria              Nitrospira 
                    "#ABD9E9"               "#74ADD1"               "#4575B4" 
                Sulfurifustis                   Other 
                    "#313695"             "lightgrey" 

---

    Code
      res$data
    Output
                               Genus SAMPLE_ID       value
                              <fctr>    <char>       <num>
       1:            Woesearchaeales      S100 0.020220588
       2:               Sideroxydans      S100 0.192095588
       3:                 uncultured      S100 0.000000000
       4:              Acinetobacter      S100 0.167279412
       5:     Candidatus_Omnitrophus      S100 0.051470588
       6:                Pseudomonas      S100 0.090992647
       7:             Nitrosarchaeum      S100 0.063419118
       8:    Candidatus_Peribacteria      S100 0.000000000
       9:                 Nitrospira      S100 0.030330882
      10:              Sulfurifustis      S100 0.000000000
      11:              Aquabacterium      S100 0.034007353
      12:                      WOR-1      S100 0.001838235
      13:  Candidatus_Falkowbacteria      S100 0.000000000
      14: Candidatus_Kerfeldbacteria      S100 0.000000000
      15:              Parcubacteria      S100 0.000000000
      16:                      Other      S100 0.092830882
      17:            Woesearchaeales      S103 0.288659794
      18:               Sideroxydans      S103 0.000000000
      19:                 uncultured      S103 0.117525773
      20:              Acinetobacter      S103 0.022680412
      21:     Candidatus_Omnitrophus      S103 0.063917526
      22:                Pseudomonas      S103 0.000000000
      23:             Nitrosarchaeum      S103 0.000000000
      24:    Candidatus_Peribacteria      S103 0.041237113
      25:                 Nitrospira      S103 0.000000000
      26:              Sulfurifustis      S103 0.034020619
      27:              Aquabacterium      S103 0.000000000
      28:                      WOR-1      S103 0.026804124
      29:  Candidatus_Falkowbacteria      S103 0.036082474
      30: Candidatus_Kerfeldbacteria      S103 0.032989691
      31:              Parcubacteria      S103 0.017525773
      32:                      Other      S103 0.092783505
      33:            Woesearchaeales      S115 0.031490015
      34:               Sideroxydans      S115 0.208909370
      35:                 uncultured      S115 0.045314900
      36:              Acinetobacter      S115 0.096774194
      37:     Candidatus_Omnitrophus      S115 0.089861751
      38:                Pseudomonas      S115 0.046850998
      39:             Nitrosarchaeum      S115 0.098310292
      40:    Candidatus_Peribacteria      S115 0.000000000
      41:                 Nitrospira      S115 0.039170507
      42:              Sulfurifustis      S115 0.012288786
      43:              Aquabacterium      S115 0.025345622
      44:                      WOR-1      S115 0.000000000
      45:  Candidatus_Falkowbacteria      S115 0.000000000
      46: Candidatus_Kerfeldbacteria      S115 0.000000000
      47:              Parcubacteria      S115 0.001536098
      48:                      Other      S115 0.145161290
      49:            Woesearchaeales      S120 0.125408942
      50:               Sideroxydans      S120 0.000000000
      51:                 uncultured      S120 0.151581243
      52:              Acinetobacter      S120 0.000000000
      53:     Candidatus_Omnitrophus      S120 0.064340240
      54:                Pseudomonas      S120 0.057797165
      55:             Nitrosarchaeum      S120 0.000000000
      56:    Candidatus_Peribacteria      S120 0.051254089
      57:                 Nitrospira      S120 0.018538713
      58:              Sulfurifustis      S120 0.037077426
      59:              Aquabacterium      S120 0.000000000
      60:                      WOR-1      S120 0.028353326
      61:  Candidatus_Falkowbacteria      S120 0.020719738
      62: Candidatus_Kerfeldbacteria      S120 0.022900763
      63:              Parcubacteria      S120 0.027262814
      64:                      Other      S120 0.109051254
                               Genus SAMPLE_ID       value
                              <fctr>    <char>       <num>

---

    Code
      res$palette
    Output
                 Woesearchaeales               Sideroxydans 
                       "#000000"                  "#004949" 
                      uncultured              Acinetobacter 
                       "#009292"                  "#ff6db6" 
          Candidatus_Omnitrophus                Pseudomonas 
                       "#ffb6db"                  "#490092" 
                  Nitrosarchaeum    Candidatus_Peribacteria 
                       "#006ddb"                  "#b66dff" 
                      Nitrospira              Sulfurifustis 
                       "#6db6ff"                  "#b6dbff" 
                   Aquabacterium                      WOR-1 
                       "#920000"                  "#924900" 
       Candidatus_Falkowbacteria Candidatus_Kerfeldbacteria 
                       "#db6d00"                  "#24ff24" 
                   Parcubacteria                      Other 
                       "#ffff6d"                "lightgrey" 

