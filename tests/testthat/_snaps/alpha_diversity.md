# Testing Alpha diversity

    Code
      res_shannon$data
    Output
               V1 CONTRAST_sex
            <num>       <char>
      1: 3.403898         male
      2: 3.776849       female
      3: 3.682609       female
      4: 3.686005         male

---

    Code
      res_shannon$stats
    Output
        .y. group1 group2 n1 n2 statistic     p p.adj p.adj.signif y.position
      1  V1 female   male  2  2         3 0.667 0.667           ns    3.78792
              groups xmin xmax
      1 female, male    1    2

---

    Code
      res_invsimpson$data
    Output
               V1 CONTRAST_sex
            <num>       <char>
      1: 22.91679         male
      2: 34.51324       female
      3: 27.09552       female
      4: 29.54116         male

---

    Code
      res_invsimpson$stats
    Output
        .y. group1 group2 n1 n2 statistic     p p.adj p.adj.signif y.position
      1  V1 female   male  2  2         3 0.667 0.667           ns   35.10964
              groups xmin xmax
      1 female, male    1    2

---

    Code
      res_simpson$data
    Output
                V1 CONTRAST_sex
             <num>       <char>
      1: 0.9563639         male
      2: 0.9710256       female
      3: 0.9630935       female
      4: 0.9661489         male

---

    Code
      res_simpson$stats
    Output
        .y. group1 group2 n1 n2 statistic     p p.adj p.adj.signif y.position
      1  V1 female   male  2  2         3 0.667 0.667           ns     0.9716
              groups xmin xmax
      1 female, male    1    2

