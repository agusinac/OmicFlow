# Testing Log2 Foldchanges

    Code
      dfe_ungrouped$data[, .SD, .SDcols = colnames(dfe_ungrouped$data)[!grepl(
        "pvalue", colnames(dfe_ungrouped$data))]]
    Output
                                                       Genus
                                                      <char>
       1:                                    Woesearchaeales
       2:                         Candidatus_Kerfeldbacteria
       3:                                    Rokubacteriales
       4:                                         Nitrospira
       5:                        Candidatus_Buchananbacteria
       6:                                      Sulfurifustis
       7:                             Candidatus_Omnitrophus
       8:                                            AKYH767
       9:                                         Babeliales
      10:                                        Subgroup_22
      11:                                  Sediminibacterium
      12:                                             DTB120
      13:                                              WOR-1
      14:                                     Berkelbacteria
      15:                                  Saccharimonadales
      16:                       Candidatus_Magasanikbacteria
      17:                                     Nitrosarchaeum
      18:                         Candidatus_Methanoperedens
      19:                                    Gracilibacteria
      20:                                             MVP-88
      21:                            Candidatus_Peribacteria
      22:                                              GAL15
      23:                Deep_Sea_Euryarchaeotic_Group(DSEG)
      24:                                      Parcubacteria
      25:                                      Lacunisphaera
      26:                                           Coxiella
      27:                          Candidatus_Kaiserbacteria
      28:                                        Pseudomonas
      29:                                      Acinetobacter
      30:                                    Polycyclovorans
      31:                                            Nevskia
      32:                                Vicinamibacteraceae
      33:                                  Nitrosopumilaceae
      34:                                   Nitrosotaleaceae
      35:                                    Marine_Group_II
      36:                            Candidatus_Azambacteria
      37:                     Methylobacterium-Methylorubrum
      38: Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
      39:                                            Hoeflea
      40:                                       Ochrobactrum
      41:                          Candidatus_Falkowbacteria
      42:                          Candidatus_Nomurabacteria
      43:                       Candidatus_Peregrinibacteria
      44:                                       Sphingopyxis
      45:                                       Sphingomonas
      46:                                    Pseudarcobacter
      47:                                         Legionella
      48:                                             MBNT15
      49:                                               WWE3
      50:                                              BSV26
      51:                                        Rhodococcus
      52:                                      Aquabacterium
      53:                                       Sideroxydans
      54:                                             KD4-96
      55:                                              OLB12
      56:                                         Lineage_IV
      57:                          Candidatus_Spechtbacteria
      58:                         Candidatus_Roizmanbacteria
      59:                                               TM7a
      60:                                               MND1
      61:                                     Anaerobacillus
      62:                                      Cutibacterium
      63:                                          0319-6G20
                                                       Genus
                                                      <char>
          Log2FC_male_vs_female_in_all         abun
                                 <num>        <num>
       1:                   -1.1364442 0.1164448349
       2:                   -0.5266196 0.0139726135
       3:                   -0.4714780 0.0097595251
       4:                    0.3191694 0.0220100256
       5:                    8.5188498 0.0013631407
       6:                   -0.3207641 0.0208467079
       7:                   -0.4090910 0.0673975263
       8:                    8.5188498 0.0013631407
       9:                  -10.3465137 0.0003840246
      10:                    6.5928504 0.0051799346
      11:                   -7.3465137 0.0030721966
      12:                    2.6212130 0.0041202026
      13:                    0.1716904 0.0142489213
      14:                    2.0906983 0.0030296901
      15:                    9.5025003 0.0006893382
      16:                   -0.6440919 0.0090828295
      17:                   -0.6324247 0.0404323524
      18:                   -9.7615512 0.0005760369
      19:                  -10.3465137 0.0003840246
      20:                   -9.7615512 0.0005760369
      21:                    0.3137238 0.0231228007
      22:                   -9.3465137 0.0007680492
      23:                   -9.7615512 0.0005760369
      24:                    0.5162447 0.0115811713
      25:                   -8.5391588 0.0013440860
      26:                    9.8407779 0.0005452563
      27:                   -9.9218409 0.0005154639
      28:                    1.6671240 0.0489102025
      29:                    0.4857974 0.0716835044
      30:                    8.7655347 0.0011488971
      31:                   -6.4885327 0.0055683564
      32:                   -9.3465137 0.0007680492
      33:                   -7.1144860 0.0036082474
      34:                   -7.0985862 0.0036482335
      35:                   -8.0245856 0.0019201229
      36:                   -8.0245856 0.0019201229
      37:                    9.0874628 0.0009191176
      38:                   -7.8870821 0.0021121352
      39:                    8.2558154 0.0016357688
      40:                    8.0874628 0.0018382353
      41:                   -0.8002925 0.0142005531
      42:                  -10.3465137 0.0003840246
      43:                    9.5025003 0.0006893382
      44:                    7.9175378 0.0020680147
      45:                    7.7655347 0.0022977941
      46:                    6.6951454 0.0048253676
      47:                   -9.3465137 0.0007680492
      48:                   -1.2408651 0.0073356643
      49:                   -9.3368784 0.0007731959
      50:                    8.2801079 0.0016084559
      51:                    8.5025003 0.0013786765
      52:                    0.4241101 0.0148382438
      53:                   -0.1210528 0.1002512396
      54:                    8.7655347 0.0011488971
      55:                   -7.5391588 0.0026881720
      56:                    9.2558154 0.0008178844
      57:                   -8.3368784 0.0015463918
      58:                   -9.3368784 0.0007731959
      59:                   -7.6460740 0.0024961598
      60:                    7.9175378 0.0020680147
      61:                   -6.7026575 0.0048003072
      62:                    7.6280312 0.0025275735
      63:                    9.8407779 0.0005452563
          Log2FC_male_vs_female_in_all         abun
                                 <num>        <num>

---

    Code
      dfe_grouped$data[, .SD, .SDcols = colnames(dfe_grouped$data)[!grepl("pvalue",
        colnames(dfe_grouped$data))]]
    Output
                                                       Genus
                                                      <char>
       1:                                    Woesearchaeales
       2:                         Candidatus_Kerfeldbacteria
       3:                                    Rokubacteriales
       4:                                         Nitrospira
       5:                        Candidatus_Buchananbacteria
       6:                                      Sulfurifustis
       7:                             Candidatus_Omnitrophus
       8:                                            AKYH767
       9:                                         Babeliales
      10:                                        Subgroup_22
      11:                                  Sediminibacterium
      12:                                             DTB120
      13:                                              WOR-1
      14:                                     Berkelbacteria
      15:                                  Saccharimonadales
      16:                       Candidatus_Magasanikbacteria
      17:                                     Nitrosarchaeum
      18:                         Candidatus_Methanoperedens
      19:                                    Gracilibacteria
      20:                                             MVP-88
      21:                            Candidatus_Peribacteria
      22:                                              GAL15
      23:                Deep_Sea_Euryarchaeotic_Group(DSEG)
      24:                                      Parcubacteria
      25:                                      Lacunisphaera
      26:                                           Coxiella
      27:                          Candidatus_Kaiserbacteria
      28:                                        Pseudomonas
      29:                                      Acinetobacter
      30:                                    Polycyclovorans
      31:                                            Nevskia
      32:                                Vicinamibacteraceae
      33:                                  Nitrosopumilaceae
      34:                                   Nitrosotaleaceae
      35:                                    Marine_Group_II
      36:                            Candidatus_Azambacteria
      37:                     Methylobacterium-Methylorubrum
      38: Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
      39:                                            Hoeflea
      40:                                       Ochrobactrum
      41:                          Candidatus_Falkowbacteria
      42:                          Candidatus_Nomurabacteria
      43:                       Candidatus_Peregrinibacteria
      44:                                       Sphingopyxis
      45:                                       Sphingomonas
      46:                                    Pseudarcobacter
      47:                                         Legionella
      48:                                             MBNT15
      49:                                               WWE3
      50:                                              BSV26
      51:                                        Rhodococcus
      52:                                      Aquabacterium
      53:                                       Sideroxydans
      54:                                             KD4-96
      55:                                              OLB12
      56:                                         Lineage_IV
      57:                          Candidatus_Spechtbacteria
      58:                         Candidatus_Roizmanbacteria
      59:                                               TM7a
      60:                                               MND1
      61:                                     Anaerobacillus
      62:                                      Cutibacterium
      63:                                          0319-6G20
                                                       Genus
                                                      <char>
          Log2FC_tumor_vs_healthy_in_male Log2FC_tumor_vs_healthy_in_female
                                    <num>                             <num>
       1:                      -3.0013382                        -3.0013382
       2:                      -5.1612543                        -5.1612543
       3:                      -5.6789733                        -5.6789733
       4:                       1.9065007                         1.9065007
       5:                      -8.5188498                        -8.5188498
       6:                      -2.5324674                        -2.5324674
       7:                       0.1400454                         0.1400454
       8:                      -8.5188498                        -8.5188498
       9:                      10.3465137                        10.3465137
      10:                      -6.5928504                        -6.5928504
      11:                       7.3465137                         7.3465137
      12:                      -2.6212130                        -2.6212130
      13:                      -4.9071624                        -4.9071624
      14:                      -2.0906983                        -2.0906983
      15:                       9.5025003                         9.5025003
      16:                      -4.5015467                        -4.5015467
      17:                       3.6283460                         3.6283460
      18:                       9.7615512                         9.7615512
      19:                      10.3465137                        10.3465137
      20:                       9.7615512                         9.7615512
      21:                      -4.4345400                        -4.4345400
      22:                       9.3465137                         9.3465137
      23:                       9.7615512                         9.7615512
      24:                      -4.8657887                        -4.8657887
      25:                       8.5391588                         8.5391588
      26:                      -9.8407779                        -9.8407779
      27:                      -9.9218409                        -9.9218409
      28:                       1.2539621                         1.2539621
      29:                       3.5413121                         3.5413121
      30:                       8.7655347                         8.7655347
      31:                       6.4885327                         6.4885327
      32:                       9.3465137                         9.3465137
      33:                      -7.1144860                        -7.1144860
      34:                       7.0985862                         7.0985862
      35:                       8.0245856                         8.0245856
      36:                       8.0245856                         8.0245856
      37:                       9.0874628                         9.0874628
      38:                       7.8870821                         7.8870821
      39:                      -8.2558154                        -8.2558154
      40:                       8.0874628                         8.0874628
      41:                      -5.1379091                        -5.1379091
      42:                      10.3465137                        10.3465137
      43:                       9.5025003                         9.5025003
      44:                       7.9175378                         7.9175378
      45:                       7.7655347                         7.7655347
      46:                       6.6951454                         6.6951454
      47:                       9.3465137                         9.3465137
      48:                      -6.0908567                        -6.0908567
      49:                      -9.3368784                        -9.3368784
      50:                       8.2801079                         8.2801079
      51:                       8.5025003                         8.5025003
      52:                       5.0745358                         5.0745358
      53:                       2.3183080                         2.3183080
      54:                       8.7655347                         8.7655347
      55:                       7.5391588                         7.5391588
      56:                      -9.2558154                        -9.2558154
      57:                      -8.3368784                        -8.3368784
      58:                      -9.3368784                        -9.3368784
      59:                       7.6460740                         7.6460740
      60:                       7.9175378                         7.9175378
      61:                       6.7026575                         6.7026575
      62:                       7.6280312                         7.6280312
      63:                      -9.8407779                        -9.8407779
          Log2FC_tumor_vs_healthy_in_male Log2FC_tumor_vs_healthy_in_female
                                    <num>                             <num>
                  abun
                 <num>
       1: 0.1164448349
       2: 0.0139726135
       3: 0.0097595251
       4: 0.0220100256
       5: 0.0013631407
       6: 0.0208467079
       7: 0.0673975263
       8: 0.0013631407
       9: 0.0003840246
      10: 0.0051799346
      11: 0.0030721966
      12: 0.0041202026
      13: 0.0142489213
      14: 0.0030296901
      15: 0.0006893382
      16: 0.0090828295
      17: 0.0404323524
      18: 0.0005760369
      19: 0.0003840246
      20: 0.0005760369
      21: 0.0231228007
      22: 0.0007680492
      23: 0.0005760369
      24: 0.0115811713
      25: 0.0013440860
      26: 0.0005452563
      27: 0.0005154639
      28: 0.0489102025
      29: 0.0716835044
      30: 0.0011488971
      31: 0.0055683564
      32: 0.0007680492
      33: 0.0036082474
      34: 0.0036482335
      35: 0.0019201229
      36: 0.0019201229
      37: 0.0009191176
      38: 0.0021121352
      39: 0.0016357688
      40: 0.0018382353
      41: 0.0142005531
      42: 0.0003840246
      43: 0.0006893382
      44: 0.0020680147
      45: 0.0022977941
      46: 0.0048253676
      47: 0.0007680492
      48: 0.0073356643
      49: 0.0007731959
      50: 0.0016084559
      51: 0.0013786765
      52: 0.0148382438
      53: 0.1002512396
      54: 0.0011488971
      55: 0.0026881720
      56: 0.0008178844
      57: 0.0015463918
      58: 0.0007731959
      59: 0.0024961598
      60: 0.0020680147
      61: 0.0048003072
      62: 0.0025275735
      63: 0.0005452563
                  abun
                 <num>

