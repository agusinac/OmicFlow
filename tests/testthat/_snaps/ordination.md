# Testing UniFrac ordination

    Code
      res$anova_data
    Output
                 pairs Df  SumsOfSqs   F.Model         R2   p.value     p.adj
      1 male vs female  1 0.01153988 0.2067667 0.09369667 0.6666667 0.6666667

---

    Code
      res$dist
    Output
                S100      S103      S115      S120
      S100 0.0000000 0.3722812 0.1069508 0.3381820
      S103 0.3722812 0.0000000 0.3299658 0.1645516
      S115 0.1069508 0.3299658 0.0000000 0.3038010
      S120 0.3381820 0.1645516 0.3038010 0.0000000

---

    Code
      res$pcs
    Output
                PC1          PC2          PC3 groups samples
              <num>        <num>        <num> <char>  <char>
      1:  0.1823439  0.006040795  0.044981932   male       1
      2: -0.1794024 -0.074588820  0.009885969 female       2
      3:  0.1397769 -0.016636934 -0.050476140 female       3
      4: -0.1427184  0.085184959 -0.004391762   male       4

