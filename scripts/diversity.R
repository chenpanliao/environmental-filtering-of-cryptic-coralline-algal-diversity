library(magrittr)
library(data.table)
library(iNEXT)
library(vegan)
library(ggrepel)
library(ggpubr)
library(dendextend)
library(ggtree)
library(ape)
library(adespatial)
sessionInfo()
# R version 4.4.2 (2024-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: Asia/Taipei
# tzcode source: internal
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] adespatial_0.3-24 ape_5.8           ggtree_3.12.0     dendextend_1.18.1 ggpubr_0.6.0      ggrepel_0.9.6     ggplot2_3.5.1    
#  [8] vegan_2.6-8       lattice_0.22-6    permute_0.9-7     iNEXT_3.0.1       coin_1.4-3        survival_3.7-0    data.table_1.16.2
# [15] magrittr_2.0.3   
# 
# loaded via a namespace (and not attached):
#   [1] libcoin_1.0-10      RColorBrewer_1.1-3  wk_0.9.4            rstudioapi_0.17.0   jsonlite_1.8.9      TH.data_1.1-2      
#   [7] modeltools_0.2-23   farver_2.1.2        fs_1.6.4            adegraphics_1.0-21  vctrs_0.6.5         spdep_1.3-6        
#  [13] rstatix_0.7.2       htmltools_0.5.8.1   progress_1.2.3      broom_1.0.7         s2_1.1.7            Formula_1.2-5      
#  [19] gridGraphics_0.5-1  adegenet_2.1.10     spData_2.3.3        KernSmooth_2.23-24  adephylo_1.1-16     plyr_1.8.9         
#  [25] sandwich_3.1-1      zoo_1.8-12          uuid_1.2-1          igraph_2.0.3        mime_0.12           lifecycle_1.0.4    
#  [31] pkgconfig_2.0.3     Matrix_1.7-1        R6_2.5.1            fastmap_1.2.0       shiny_1.9.1         digest_0.6.37      
#  [37] aplot_0.2.3         colorspace_2.1-1    patchwork_1.3.0     phylobase_0.8.12    fansi_1.0.6         httr_1.4.7         
#  [43] abind_1.4-8         mgcv_1.9-1          compiler_4.4.2      proxy_0.4-27        withr_3.0.1         backports_1.5.0    
#  [49] carData_3.0-5       viridis_0.6.5       DBI_1.2.3           ggsignif_0.6.4      MASS_7.3-61         classInt_0.4-10    
#  [55] tools_4.4.2         units_0.8-5         rncl_0.8.7          httpuv_1.6.15       glue_1.8.0          nlme_3.1-166       
#  [61] promises_1.3.0      grid_4.4.2          sf_1.0-18           cluster_2.1.6       reshape2_1.4.4      ade4_1.7-22        
#  [67] generics_0.1.3      seqinr_4.2-36       gtable_0.3.5        class_7.3-22        tidyr_1.3.1         hms_1.1.3          
#  [73] sp_2.1-4            xml2_1.3.6          car_3.1-3           utf8_1.2.4          pillar_1.9.0        stringr_1.5.1      
#  [79] yulab.utils_0.1.7   later_1.3.2         splines_4.4.2       dplyr_1.1.4         treeio_1.28.0       deldir_2.0-4       
#  [85] tidyselect_1.2.1    gridExtra_2.3       stats4_4.4.2        matrixStats_1.4.1   stringi_1.8.4       boot_1.3-31        
#  [91] lazyeval_0.2.2      ggfun_0.1.6         codetools_0.2-20    interp_1.1-6        tibble_3.2.1        ggplotify_0.1.2    
#  [97] cli_3.6.3           xtable_1.8-4        munsell_0.5.1       Rcpp_1.0.13         png_0.1-8           XML_3.99-0.17      
# [103] parallel_4.4.2      RNeXML_2.4.11       prettyunits_1.2.0   latticeExtra_0.6-30 jpeg_0.1-10         viridisLite_0.4.2  
# [109] mvtnorm_1.3-1       tidytree_0.4.6      scales_1.3.0        e1071_1.7-16        purrr_1.0.2         crayon_1.5.3       

#### Data ####
d0 <- fread("../data/mOTU-composition-8sites-season.csv")
dd.8sites <- 
  melt(d0, id.vars = c("area", "site", "season"), variable.name = "mOTU", value.name = "count") %>% 
  .[, .(countSum = sum(count)), by = .(area, site, mOTU)] %>% 
  dcast(area + site ~ mOTU, value.var = "countSum") %>% 
  setkey(area, site)
dd.8sites.Y <- dd.8sites[, -1:-2] %>% as.data.frame %>% set_rownames(dd.8sites$site) 

#### Benthic composition of functional groups differed between the TAR (represented by Datan G2) and the NTCR (represented by Shimen) ####
d.coverage <- 
  fread("../data/Table-S2.csv") %>% 
  setkey(site, season, replicate) %>% 
  .[, Others := 100 - CCA - FM - TM - ACA - Zoantharia]
myPerm <- how(plots = Plots(d.coverage$season), nperm = 4999)
adonis2(d.coverage[, CCA:Others] %>% as.matrix ~ site, data = d.coverage, 
        permutations = myPerm,
        by = "margin")
#          Df SumOfSqs      R2      F Pr(>F)    
# site      1   1.3345 0.53663 18.529  2e-04 ***
# Residual 16   1.1523 0.46337                  
# Total    17   2.4867 1.00000                  
    

#### Differences in mOTU composition between the TAR and the NTCR; PERMANOVA ####
adonis2(dd.8sites.Y %>% decostand("hellinger") %>% vegdist ~ area, 
        data = dd.8sites, 
        permutations = 4999)
#          Df SumOfSqs      R2      F Pr(>F)  
# Model     1  0.56707 0.39536 3.9233 0.0364 *
# Residual  6  0.86724 0.60464                
# Total     7  1.43431 1.00000                
        

#### Rarefaction curves at Datan G2, Shimen, TAR and NTCR ####

# "Shimen", "Datan G2" only
dd21 <- 
  dd.8sites[site %in% c("Shimen", "Datan G2")] %>% 
  melt(
    id.vars = "site",
    measure.vars = grep("mOTU[0-9]{3}", names(.)),
    variable.name = "mOTU",
    value.name = "count"
  ) %>% 
  .[, .(countSum = sum(count)), by = .(site, mOTU)] %>% 
  dcast(site ~ mOTU, value.var = "countSum")
dd21.Y <- dd21[, -1] %>% as.data.frame %>% set_rownames(dd21$site) 
reCurve21 <-
  lapply(seq_len(nrow(dd21.Y)), function(i) {dd21.Y[i,] %>% unlist %>% sort(decreasing = T)} %>% {.[.>0]}) %>% 
  set_names(dd21$site) %>% 
  iNEXT(q = c(0,1,2), datatype = "abundance", endpoint = 500, knots = 100, nboot = 500)

# pooling TAR and NTCR
dd22 <- 
  dd.8sites %>%
  melt(
    id.vars = "area",
    measure.vars = grep("mOTU[0-9]{3}", names(.)),
    variable.name = "mOTU",
    value.name = "count"
  ) %>% 
  .[, .(countSum = sum(count)), by = .(area, mOTU)] %>% 
  dcast(area ~ mOTU, value.var = "countSum") %>% {
    # exclude all zero mOTUs
    foo <- .
    cbind(foo[, .(area)], foo[, which(colSums(foo[, -1]) > 0) %>% names, with = F])
  }
dd22.Y <- dd22[, -1] %>% as.data.frame %>% set_rownames(dd22$area) 
reCurve22 <-
  lapply(seq_len(nrow(dd22.Y)), function(i) {dd22.Y[i,] %>% unlist %>% sort(decreasing = T)} %>% {.[.>0]}) %>% 
  set_names(dd22$area) %>% 
  iNEXT(q = c(0,1,2), datatype = "abundance", endpoint = 500, knots = 250, nboot = 500)

reCurve.fig <- 
  rbind(
    reCurve21$iNextEst$size_based %>% as.data.table, 
    reCurve22$iNextEst$size_based %>% as.data.table
  ) %>% 
  .[, Method := factor(Method, levels = c("Rarefaction", "Observed", "Extrapolation"))] %>% 
  .[, Order.q := factor(Order.q)] %>% 
  .[, Assemblage := factor(Assemblage, 
                           levels = c("TAR (Taoyuan)", "Datan G2", "NTCR (New Taipei)", "Shimen"), 
                           labels = c("TAR", "Datan G2", "NTCR", "Shimen"))]

reCurve.fig[Order.q == "0" & Method != "Extrapolation"] %>%
  ggplot(aes(m, qD)) +
  theme_pubr(8, border = T) +
  theme(legend.text = element_text(size = 8), strip.text = element_text(size = 8)) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, fill = Assemblage, group = Assemblage), alpha = 0.2) +
  geom_line(
    data = reCurve.fig[Order.q == "0" & Method == "Rarefaction"],
    aes(
      linewidth = Assemblage,
      group = interaction(Assemblage, Method), color = Assemblage)
  ) +
  geom_point(
    data = reCurve.fig[Order.q == 0 & Method == "Observed"],
    aes(color = Assemblage, shape = Assemblage)
  ) +
  scale_fill_manual("Assemblage", values = c("#FE6100", "#FFB000", "#785EF0", "#648FFF")) +
  scale_color_manual("Assemblage", values = c("#FE6100", "#FFB000", "#785EF0", "#648FFF")) +
  scale_linetype_manual(values = 1:2) +
  scale_linewidth_manual(values = c(1, 0.5, 1, 0.5)) +
  scale_shape_manual(values = c(15:18)) +
  scale_y_continuous("Richness") +
  scale_x_continuous("Sample size") +
  guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2)) 

#### relative abundant ####
fig.rankAbundance.data <- 
  rbind(
    dd21 %>% copy %>% setnames("site", "Assemblage"),
    dd22 %>% copy %>% setnames("area", "Assemblage")
  ) %>% 
  .[, Assemblage := factor(Assemblage, 
                           levels = c("TAR (Taoyuan)", "Datan G2", "NTCR (New Taipei)", "Shimen"), 
                           labels = c("TAR", "Datan G2", "NTCR", "Shimen"))]
fig.rankAbundance.data.labels <-
  fig.rankAbundance.data[, -1] %>% as.matrix %>% set_rownames(fig.rankAbundance.data$Assemblage) %>% 
  decostand("total") %>% 
  as.data.table(keep.rownames = "Assemblage") %>% 
  melt(id.vars = "Assemblage", variable.name = "mOTU", value = "Rel. abundance") %>% 
  .[`Rel. abundance` > 0] %>% 
  .[, mOTU := gsub("^mOTU", "", mOTU)] %>% 
  .[, mOTU := gsub("023", "(023)", mOTU)] %>% 
  .[, mOTU := gsub("014", "(014)", mOTU)] %>% 
  .[, mOTU := gsub("066", "(066)", mOTU)] %>% 
  .[, mOTU := gsub("139", "(139)", mOTU)] %>% 
  .[, Assemblage := factor(Assemblage, 
                           levels = c("TAR", "Datan G2", "NTCR", "Shimen"), 
                           labels = c("TAR", "Datan G2", "NTCR", "Shimen"))] %>% 
  .[order(Assemblage, `Rel. abundance`, decreasing = T)] %>% 
  .[, `Ranked species number` := 1:.N, by = .(Assemblage)]

fig.rankAbundance.data.labels %>% 
  ggplot(aes(`Ranked species number`, `Rel. abundance`)) +
  geom_point(aes(color = Assemblage, shape = Assemblage), size = 1) +
  geom_line(aes(color = Assemblage, group = Assemblage, linewidth = Assemblage)) +
  theme_pubr(8, border = T) +
  theme(legend.text = element_text(size = 8), strip.text = element_text(size = 8)) +
  scale_color_manual("Area/Site", values = c("#FE6100", "#FFB000", "#785EF0", "#648FFF")) +
  scale_linetype_manual("Area/Site", values = 1:2) +
  scale_shape_manual("Area/Site", values = c(15:18)) +
  scale_linewidth_manual("Area/Site", values = c(1, 0.5, 1, 0.5)) +
  geom_text_repel(
    data = fig.rankAbundance.data.labels[`Ranked species number` <= 2 | grepl("^[^0-9]", mOTU)],
    aes(label = mOTU, color = Assemblage), 
    size = 8*25.4/72,
    show.legend = F
  ) +
  scale_x_continuous(limits = c(1, 7), breaks = 1:20, expand = expansion(c(0.15, 0.05))) +
  scale_y_continuous(breaks = seq(0,1,0.1), expand = expansion(c(0.15, 0.15))) +
  guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2)) +
  theme(legend.direction = "vertical")

#### heatmap ####
dd.8sites
dd.8sites.Y
sample.tree <- 
  hclust(dd.8sites.Y %>% decostand("hellinger") %>% vegdist("bray"), method = "complete") %>% 
  dendextend::rotate(order = c("Cianshueiwan", "Shimen", "Yongxin", "Datan G1", "Datan G2", "Baiyu", "Baosheng", "Yongan") %>% rev)
f.sample.tree <-
  (ggtree(sample.tree, right = F, ladderize = F) %<+% dd.8sites[, .(site = site)]) + 
  geom_tiplab(offset = 0, align = T, hjust = 0, vjust = 0.5, size = 8*25.4/72, color = 1) +
  coord_flip() +
  scale_x_continuous(transform = "reverse")
sample.order <- sample.tree$labels[sample.tree$order]

# mOTUs by total count
totalCount <- dd.8sites.Y %>% colSums %>% sort(decreasing = T) %>% {. / sum(.)} %T>% print
#     mOTU023     mOTU139     mOTU014     mOTU075     mOTU116     mOTU025     mOTU055     mOTU097     mOTU104     mOTU018 
# 0.318037975 0.183544304 0.075949367 0.063291139 0.044303797 0.034810127 0.028481013 0.025316456 0.022151899 0.018987342 
#     mOTU008     mOTU027     mOTU066     mOTU114     mOTU117     mOTU103     mOTU107     mOTU065     mOTU098     mOTU056 
# 0.017405063 0.015822785 0.014240506 0.011075949 0.011075949 0.009493671 0.009493671 0.007911392 0.007911392 0.006329114 
#     mOTU058     mOTU073     mOTU100     mOTU072     mOTU110     mOTU142     mOTU001     mOTU002     mOTU011     mOTU040 
# 0.006329114 0.006329114 0.006329114 0.004746835 0.004746835 0.004746835 0.003164557 0.003164557 0.003164557 0.003164557 
#     mOTU053     mOTU089     mOTU106     mOTU112     mOTU147     mOTU151     mOTU015     mOTU105     mOTU108     mOTU111 
# 0.003164557 0.003164557 0.003164557 0.003164557 0.003164557 0.003164557 0.001582278 0.001582278 0.001582278 0.001582278 
#     mOTU146     mOTU150 
# 0.001582278 0.001582278 
mOTUs.names <- colnames(dd.8sites.Y)
mOTU.used <- names(totalCount)[1:which(names(totalCount) == "mOTU066")]
mOTU.used.not <- mOTUs.names[!(mOTUs.names %in% mOTU.used)]

Y1 <- cbind(dd.8sites.Y[, mOTU.used], Others = dd.8sites.Y[, !(mOTUs.names %in% mOTU.used)] %>% rowSums)
Y1.relative <- Y1 %>% decostand("total")
d1 <- cbind(dd.8sites[, !(names(dd.8sites) %in% mOTUs.names), with = F], Y1)
d1.relative <- cbind(dd.8sites[, !(names(dd.8sites) %in% mOTUs.names), with = F], Y1.relative)


phy <- read.tree("../data/mOTU.only.removeSampleCode.nwk")
phy.drop <- phy %>% drop.tip(phy$tip.label[!(phy$tip.label %in% mOTU.used)]) %T>% plot
phy.drop$tip.label <- gsub("^mOTU", "", phy.drop$tip.label)
phy.drop %>% plot
f.mOTU.tree <-
  ggtree(phy.drop, right = F, ladderize = F) + 
  geom_tiplab(offset = 0.01, align = T, hjust = -0.5, vjust = 0.5, angle = 90, size = 2, color = 8) +
  # coord_flip() +
  scale_y_discrete(expand = expansion(c(0.02, 0.08)))
f.mOTU.tree
mOTU.order <- phy.drop$tip.label


d1.long.count <- 
  melt(d1, 
       id.vars = c("site"), 
       measure.vars = c(mOTU.used, "Others"),
       variable.name = "mOTU", value.name = "count") %>% 
  .[, mOTU := factor(mOTU, levels = c(mOTU.used, "Others"), labels = gsub("^mOTU", "", c(mOTU.used, "Others (29 mOTUs)")))]
d1.long.relative <- 
  melt(d1.relative, 
       id.vars = c("site"), 
       measure.vars = c(mOTU.used, "Others"),
       variable.name = "mOTU", value.name = "count.relative") %>% 
  .[, mOTU := factor(mOTU, levels = c(mOTU.used, "Others"), labels = gsub("^mOTU", "", c(mOTU.used, "Others (29 mOTUs)")))]
d1.long <- merge(d1.long.count, d1.long.relative)

colors <- gray(seq(1, 0, length.out = 100))
f.heatmap <- 
  ggplot(d1.long, aes(site, mOTU)) +
  geom_tile(aes(fill = count.relative), color = gray(0.7)) +
  geom_text(data = d1.long[count > 0], aes(label = count), 
            color = gray(1), size = 6.5*25.4/72, family = "mono", hjust = 1, vjust = 0.5, position = position_nudge(x = 0.4)) +
  theme_pubr(8, border = T, legend = "bottom") +
  scale_y_discrete(
    limits = c("Others (29 mOTUs)", mOTU.order %>% rev), 
    expand = expansion(c(0,0)),
    labels = 
      c("Others (29 mOTUs)", mOTU.order %>% rev) %>% 
      gsub("139", "Sporolithon sp. (139)", .) %>% 
      gsub("023", "Harveylithon sp. (023)", .) %>% 
      gsub("014", "H. catarinense (014)", .) %>% 
      gsub("066", "Lithophyllum sp. (066)", .) 
    ) +
  scale_x_discrete("Site", limits = sample.order, expand = expansion(c(0,0)), position = "top") +
  scale_fill_gradientn(
    "Relative abundance", 
    transform = "sqrt",
    breaks = c(0,0.01,0.05, seq(0.1,1,0.1)),
    colours = colors
  ) +
  theme(legend.key.width = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.caption = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.text.x.top = element_text(size = 8, hjust = 0, vjust = 0.5)
        ) +
  theme(plot.title = element_text(size = 8))

deeptime::ggarrange2(
  f.mOTU.tree,
  f.heatmap,
  f.sample.tree,
  layout = matrix(c(0, 3, 1, 2), nrow = 2, byrow = T),
  byrow = T,
  nrow = 2,
  ncol = 2,
  heights = c(1, 5),
  widths = c(1, 3)
)
