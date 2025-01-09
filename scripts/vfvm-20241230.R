library(magrittr)
library(data.table)
library(coin)

sessionInfo()
# R version 4.4.2 (2024-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: Asia/Taipei
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] coin_1.4-3        survival_3.7-0    data.table_1.16.2 magrittr_2.0.3   
# 
# loaded via a namespace (and not attached):
#  [1] codetools_0.2-20  multcomp_1.4-26   Matrix_1.7-1      lattice_0.22-6    TH.data_1.1-2    
#  [6] splines_4.4.2     zoo_1.8-12        matrixStats_1.4.1 modeltools_0.2-23 libcoin_1.0-10   
# [11] parallel_4.4.2    stats4_4.4.2      mvtnorm_1.3-1     sandwich_3.1-1    grid_4.4.2       
# [16] compiler_4.4.2    rstudioapi_0.17.0 tools_4.4.2       MASS_7.3-61      

dd <- fread("../data/Table-S6.csv", stringsAsFactors = T)
dd %>%
  setkey(Species, SandBurialDuration, RecoveryDuration) %>%
  .[, .(p =
          wilcox_test(fvfm ~ Treatment, distribution = "exact") %>%
          pvalue %>%
          as.double), by = .(Species, SandBurialDuration, RecoveryDuration)] %>%
  .[, p.adj.fdr := p.adjust(p, "fdr"), by = .(Species)] %>%
  .[, p.adj.fdr.rank :=
      cut(
        p.adj.fdr,
        c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
        labels = c("****", "***", "**", "*", "NS")
      )] %>%
  print
# Key: <Species, SandBurialDuration, RecoveryDuration>
#                            Species SandBurialDuration RecoveryDuration           p   p.adj.fdr p.adj.fdr.rank
#                             <fctr>             <fctr>            <int>       <num>       <num>         <fctr>
#  1: Hs-TAR-NTCR (Harveylithon sp.)                 1m                0 0.028571429 0.048351648              *
#  2: Hs-TAR-NTCR (Harveylithon sp.)                 1m                1 0.028571429 0.048351648              *
#  3: Hs-TAR-NTCR (Harveylithon sp.)                 1m                2 0.028571429 0.048351648              *
#  4: Hs-TAR-NTCR (Harveylithon sp.)                 1m                3 0.028571429 0.048351648              *
#  5: Hs-TAR-NTCR (Harveylithon sp.)                 1m                4 0.028571429 0.048351648              *
#  6: Hs-TAR-NTCR (Harveylithon sp.)                 1m                5 0.028571429 0.048351648              *
#  7: Hs-TAR-NTCR (Harveylithon sp.)                 1m                6 0.342857143 0.359183673             NS
#  8: Hs-TAR-NTCR (Harveylithon sp.)                 1m                7 0.057142857 0.083809524             NS
#  9: Hs-TAR-NTCR (Harveylithon sp.)                 1m               14 0.142857143 0.184873950             NS
# 10: Hs-TAR-NTCR (Harveylithon sp.)                 1m               21 0.028571429 0.048351648              *
# 11: Hs-TAR-NTCR (Harveylithon sp.)                 1m               28 0.028571429 0.048351648              *
# 12: Hs-TAR-NTCR (Harveylithon sp.)                 3m                0 0.028571429 0.048351648              *
# 13: Hs-TAR-NTCR (Harveylithon sp.)                 3m                1 0.028571429 0.048351648              *
# 14: Hs-TAR-NTCR (Harveylithon sp.)                 3m                2 0.028571429 0.048351648              *
# 15: Hs-TAR-NTCR (Harveylithon sp.)                 3m                3 0.028571429 0.048351648              *
# 16: Hs-TAR-NTCR (Harveylithon sp.)                 3m                4 0.200000000 0.244444444             NS
# 17: Hs-TAR-NTCR (Harveylithon sp.)                 3m                5 0.485714286 0.485714286             NS
# 18: Hs-TAR-NTCR (Harveylithon sp.)                 3m                6 0.028571429 0.048351648              *
# 19: Hs-TAR-NTCR (Harveylithon sp.)                 3m                7 0.314285714 0.345714286             NS
# 20: Hs-TAR-NTCR (Harveylithon sp.)                 3m               14 0.285714286 0.330827068             NS
# 21: Hs-TAR-NTCR (Harveylithon sp.)                 3m               21 0.057142857 0.083809524             NS
# 22: Hs-TAR-NTCR (Harveylithon sp.)                 3m               28 0.114285714 0.157142857             NS
# 23:          Lo-NTCR (L. okamurae)                 1m                0 0.015873016 0.024943311              *
# 24:          Lo-NTCR (L. okamurae)                 1m                1 0.031746032 0.043650794              *
# 25:          Lo-NTCR (L. okamurae)                 1m                2 0.015873016 0.024943311              *
# 26:          Lo-NTCR (L. okamurae)                 1m                3 0.015873016 0.024943311              *
# 27:          Lo-NTCR (L. okamurae)                 1m                4 0.015873016 0.024943311              *
# 28:          Lo-NTCR (L. okamurae)                 1m                5 0.015873016 0.024943311              *
# 29:          Lo-NTCR (L. okamurae)                 1m                6 0.373015873 0.410317460             NS
# 30:          Lo-NTCR (L. okamurae)                 1m                7 0.111111111 0.128654971             NS
# 31:          Lo-NTCR (L. okamurae)                 1m               14 0.015873016 0.024943311              *
# 32:          Lo-NTCR (L. okamurae)                 1m               21 0.015873016 0.024943311              *
# 33:          Lo-NTCR (L. okamurae)                 1m               28 0.111111111 0.128654971             NS
# 34:          Lo-NTCR (L. okamurae)                 3m                0 0.015873016 0.024943311              *
# 35:          Lo-NTCR (L. okamurae)                 3m                1 0.015873016 0.024943311              *
# 36:          Lo-NTCR (L. okamurae)                 3m                2 0.015873016 0.024943311              *
# 37:          Lo-NTCR (L. okamurae)                 3m                3 0.015873016 0.024943311              *
# 38:          Lo-NTCR (L. okamurae)                 3m                4 0.015873016 0.024943311              *
# 39:          Lo-NTCR (L. okamurae)                 3m                5 0.015873016 0.024943311              *
# 40:          Lo-NTCR (L. okamurae)                 3m                6 0.063492063 0.082166200             NS
# 41:          Lo-NTCR (L. okamurae)                 3m                7 0.031746032 0.043650794              *
# 42:          Lo-NTCR (L. okamurae)                 3m               14 0.015873016 0.024943311              *
# 43:          Lo-NTCR (L. okamurae)                 3m               21 1.000000000 1.000000000             NS
# 44:          Lo-NTCR (L. okamurae)                 3m               28 0.730158730 0.764928193             NS
# 45:    Se-TAR-NTCR (S. erythraeum)                 1m                0 0.004040404 0.008888889             **
# 46:    Se-TAR-NTCR (S. erythraeum)                 1m                1 0.004040404 0.008888889             **
# 47:    Se-TAR-NTCR (S. erythraeum)                 1m                2 0.008080808 0.013675214              *
# 48:    Se-TAR-NTCR (S. erythraeum)                 1m                3 0.008080808 0.013675214              *
# 49:    Se-TAR-NTCR (S. erythraeum)                 1m                4 0.004040404 0.008888889             **
# 50:    Se-TAR-NTCR (S. erythraeum)                 1m                5 0.028282828 0.041481481              *
# 51:    Se-TAR-NTCR (S. erythraeum)                 1m                6 0.008080808 0.013675214              *
# 52:    Se-TAR-NTCR (S. erythraeum)                 1m                7 0.153535354 0.187654321             NS
# 53:    Se-TAR-NTCR (S. erythraeum)                 1m               14 1.000000000 1.000000000             NS
# 54:    Se-TAR-NTCR (S. erythraeum)                 1m               21 0.002020202 0.008888889             **
# 55:    Se-TAR-NTCR (S. erythraeum)                 1m               28 0.004040404 0.008888889             **
# 56:    Se-TAR-NTCR (S. erythraeum)                 3m                0 0.004040404 0.008888889             **
# 57:    Se-TAR-NTCR (S. erythraeum)                 3m                1 0.004040404 0.008888889             **
# 58:    Se-TAR-NTCR (S. erythraeum)                 3m                2 0.004040404 0.008888889             **
# 59:    Se-TAR-NTCR (S. erythraeum)                 3m                3 0.004040404 0.008888889             **
# 60:    Se-TAR-NTCR (S. erythraeum)                 3m                4 0.004040404 0.008888889             **
# 61:    Se-TAR-NTCR (S. erythraeum)                 3m                5 0.072727273 0.100000000             NS
# 62:    Se-TAR-NTCR (S. erythraeum)                 3m                6 0.569696970 0.626666667             NS
# 63:    Se-TAR-NTCR (S. erythraeum)                 3m                7 0.028282828 0.041481481              *
# 64:    Se-TAR-NTCR (S. erythraeum)                 3m               14 0.303030303 0.350877193             NS
# 65:    Se-TAR-NTCR (S. erythraeum)                 3m               21 0.775757576 0.812698413             NS
# 66:    Se-TAR-NTCR (S. erythraeum)                 3m               28 0.153535354 0.187654321             NS
#                            Species SandBurialDuration RecoveryDuration           p   p.adj.fdr p.adj.fdr.rank

