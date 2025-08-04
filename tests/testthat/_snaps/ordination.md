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
                S100      S103      S115
      S103 0.3722812                    
      S115 0.1069508 0.3299658          
      S120 0.3381820 0.1645516 0.3038010

---

    Code
      res$pcs
    Output
      Call: vegan::wcmdscale(d = distmat, k = 15, eig = TRUE)
      
            Inertia Rank
      Total  0.1232    3
      
      Results have 4 points, 3 axes
      
      Eigenvalues:
      [1] 0.10534 0.01313 0.00469
      
      Weights: Constant
      

