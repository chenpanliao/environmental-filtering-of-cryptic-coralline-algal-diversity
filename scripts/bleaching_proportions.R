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
#  [1] cli_3.6.3         knitr_1.48        rlang_1.1.4       TH.data_1.1-2     xfun_0.48         zoo_1.8-12       
#  [7] htmltools_0.5.8.1 libcoin_1.0-10    stats4_4.4.2      rmarkdown_2.28    modeltools_0.2-23 grid_4.4.2       
# [13] evaluate_1.0.1    MASS_7.3-61       fastmap_1.2.0     yaml_2.3.10       mvtnorm_1.3-1     compiler_4.4.2   
# [19] multcomp_1.4-26   codetools_0.2-20  sandwich_3.1-1    rstudioapi_0.17.0 lattice_0.22-6    digest_0.6.37    
# [25] parallel_4.4.2    splines_4.4.2     Matrix_1.7-1      tools_4.4.2       matrixStats_1.4.1

dd3 <- fread("../data/Table-S5.csv", stringsAsFactors = T)
# p.type:
#   p1 = control-1month vs buried-1month
#   p2 = control-3month vs buried-3month
#   p3 = control-1month vs control-3month
#   p4 = buried-1month vs buried-3month
dd3[SandBurialDuration == "1 month", 
    .(p1 = wilcox_test(Difference ~ Treatment, distribution = "exact") %>% pvalue %>% as.double), by = .(Species)] %>%
  merge(
    dd3[SandBurialDuration == "3 months", 
        .(p2 = wilcox_test(Difference ~ Treatment, distribution = "exact") %>% pvalue %>% as.double), 
        by = .(Species)]) %>%
  merge(
    dd3[Treatment == "Control", 
        .(p3 = wilcox_test(Difference ~ SandBurialDuration, distribution = "exact") %>% pvalue %>% as.double),
        by = .(Species)]) %>%
  merge(
    dd3[Treatment == "Buried", 
        .(p4 = wilcox_test(Difference ~ SandBurialDuration, distribution = "exact") %>% pvalue %>% as.double),
        by = .(Species)]) %>%
  melt(id.vars = 1,
       value.name = "p",
       variable.name = "p.type") %>%
  setkey(Species, p.type) %>% 
  .[, p.adjust.fdr := p.adjust(p, "fdr"), by = .(Species)] %>%
  .[, p.adjust.fdr.text := cut(p.adjust.fdr,
                               c(0, 0.001, 0.01, 0.05, 1),
                               labels = c("***", "**", "*", "NS"))] %>% 
  print

# Key: <Species, p.type>
#                            Species p.type           p p.adjust.fdr p.adjust.fdr.text
#                             <fctr> <fctr>       <num>        <num>            <fctr>
#  1:       Hc-NTCR (H. catarinense)     p1 0.002020202  0.004040404                **
#  2:       Hc-NTCR (H. catarinense)     p2 0.002020202  0.004040404                **
#  3:       Hc-NTCR (H. catarinense)     p3 0.028571429  0.038095238                 *
#  4:       Hc-NTCR (H. catarinense)     p4 0.538306138  0.538306138                NS
#  5: Hs-TAR-NTCR (Harveylithon sp.)     p1 0.742857143  0.742857143                NS
#  6: Hs-TAR-NTCR (Harveylithon sp.)     p2 0.200000000  0.304761905                NS
#  7: Hs-TAR-NTCR (Harveylithon sp.)     p3 0.228571429  0.304761905                NS
#  8: Hs-TAR-NTCR (Harveylithon sp.)     p4 0.142857143  0.304761905                NS
#  9:     Ls-NTCR (Lithophyllum sp.)     p1 1.000000000  1.000000000                NS
# 10:     Ls-NTCR (Lithophyllum sp.)     p2 0.004761905  0.012987013                 *
# 11:     Ls-NTCR (Lithophyllum sp.)     p3 1.000000000  1.000000000                NS
# 12:     Ls-NTCR (Lithophyllum sp.)     p4 0.006493506  0.012987013                 *
# 13:  Ss-TAR-NTCR (Sporolithon sp.)     p1 0.523232323  0.697643098                NS
# 14:  Ss-TAR-NTCR (Sporolithon sp.)     p2 0.731313131  0.731313131                NS
# 15:  Ss-TAR-NTCR (Sporolithon sp.)     p3 0.114285714  0.274436674                NS
# 16:  Ss-TAR-NTCR (Sporolithon sp.)     p4 0.137218337  0.274436674                NS
