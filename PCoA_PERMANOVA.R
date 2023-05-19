library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)
library(ggpubr)

# Read in phyloseq object
ps = readRDS("ps.Exp3.wcastings")

trtpalette = c("#fcb18c","#B7003F","#2d105a") #Amynthas pressure palette

#### PCoA, Arboretum/Gallistel ####
ps.A = subset_samples(ps,sample_data(ps)$site == "Arb")

# Prune away zero's
ps.A = prune_taxa(taxa_sums(ps.A) > 0, ps.A)
ps.A
sample_data(ps.A)

# Hellinger transformation
ps.A = transform_sample_counts(ps.A, function(x) (x / sum(x))^0.5 )

# PCoA, Arboretum
ps.PCoA = ordinate(ps.A, method="PCoA", distance="bray")
p.A = plot_ordination(ps.A, ps.PCoA, color = "trt", shape = "fraction")
p.A = p.A + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.A = p.A + scale_color_manual(values = trtpalette, 
                               labels = c("Control", "Low", "High"),
                               name = expression(paste(italic("Amynthas"), " pressure")))
p.A = p.A + scale_shape_manual(values = c(1, 17, 14, 4), labels = c("Bulk soil", "Free microaggregate",
                                                                   "Occluded microagg.", "Worm casting"),
                               name = "Fraction")
p.A # 800 x 225 for MS


######################## #
#### PCoA, Lakeshore ####
######################## #

ps.L = subset_samples(ps,sample_data(ps)$site == "Lak")

# Prune away zero's
ps.L = prune_taxa(taxa_sums(ps.L) > 0, ps.L)
ps.L
sample_data(ps.L)

# Hellinger transformation
ps.L = transform_sample_counts(ps.L, function(x) (x / sum(x))^0.5 )

# PCoA, Lakeshore
ps.PCoA = ordinate(ps.L, method="PCoA", distance="bray")
p.L = plot_ordination(ps.L, ps.PCoA, color = "trt", shape = "fraction")
p.L = p.L + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.L = p.L + scale_color_manual(values = trtpalette[-2],
                               labels = c("Control", "High"),
                               name = expression(paste(italic("Amynthas"), " pressure")))
p.L = p.L + scale_shape_manual(values = c(1, 17, 14, 4), labels = c("Bulk soil", "Free microaggregate",
                                                                 "Occluded microagg.", "Worm casting"),
                               name = "Fraction")
p.L # 800x225


######################### #
#### Arrange PCoA figs ####
######################### #

ggarrange(p.A, p.L,
          labels = c("A", "B"),
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "right")

#ggsave("Figures/PCoA.worm.tiff", width=7.25, height=4.5, units = "in", device='tiff', dpi=400)



######################## #
#### Statistics, Gallistel ####
######################## #
# First, run PERMANOVA

# # Create veganotu function
# veganotu = function(physeq) {
#   require("vegan")
#   OTU = otu_table(physeq)
#   if (taxa_are_rows(OTU)) {
#     OTU = t(OTU)
#   }
#   return(as(OTU, "matrix"))
# }

ps.A = subset_samples(ps.A, sample_data(ps.A)$fraction != "casting")
ps.A.veg = veganotu(ps.A)
ps.A.df = data.frame(sample_data(ps.A))
DistVar.A = vegdist(ps.A.veg, method = "bray")
adonis2(DistVar.A ~ trt*fraction, data = ps.A.df, method = "bray")

# adonis2(formula = DistVar.A ~ trt * fraction, data = ps.A.df, method = "bray")
#                 Df SumOfSqs   R2       F    Pr(>F)    
# trt            2   3.8401 0.17794 14.0383  0.001 ***
# fraction       2   0.2485 0.01152  0.9085  0.549    
# trt:fraction   4   0.2587 0.01199  0.4729  1.000    
# Residual     126  17.2334 0.79855                   
# Total        134  21.5807 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# # Create the pairwise adonis function (JUST BELOW)

## Run the pairwise adonis function for pairwise PERMANOVA
pairwise.adonis(DistVar.A, ps.A.df$trt)
#             pairs       Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
#   1   LowWorm vs NoWorm  1  1.705227 13.58682 0.1337459   0.001      0.001  **
#   2 LowWorm vs HighWorm  1  1.377285 10.19366 0.1038118   0.001      0.001  **
#   3  NoWorm vs HighWorm  1  2.677676 18.78045 0.1758791   0.001      0.001  **

# # Create the pairwise adonis function
# pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='BH',reduce=NULL,perm=999)
# {
#   co <- combn(unique(as.character(factors)),2)
#   pairs <- c()
#   Df <- c()
#   SumsOfSqs <- c()
#   F.Model <- c()
#   R2 <- c()
#   p.value <- c()
# 
#   for(elem in 1:ncol(co)){
#     if(inherits(x, 'dist')){
#       x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
#                       factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
#     }
# 
#     else  (
#       if (sim.function == 'daisy'){
#         x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
#       }
#       else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
#     )
# 
#     ad <- adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])],
#                  permutations = perm);
#     pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
#     Df <- c(Df,ad$aov.tab[1,1])
#     SumsOfSqs <- c(SumsOfSqs, ad$aov.tab[1,2])
#     F.Model <- c(F.Model,ad$aov.tab[1,4]);
#     R2 <- c(R2,ad$aov.tab[1,5]);
#     p.value <- c(p.value,ad$aov.tab[1,6])
#   }
#   p.adjusted <- p.adjust(p.value,method=p.adjust.m)
# 
#   sig = c(rep('',length(p.adjusted)))
#   sig[p.adjusted <= 0.05] <-'.'
#   sig[p.adjusted <= 0.01] <-'*'
#   sig[p.adjusted <= 0.001] <-'**'
#   sig[p.adjusted <= 0.0001] <-'***'
#   pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
# 
#   if(!is.null(reduce)){
#     pairw.res <- subset (pairw.res, grepl(reduce,pairs))
#     pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
# 
#     sig = c(rep('',length(pairw.res$p.adjusted)))
#     sig[pairw.res$p.adjusted <= 0.1] <-'.'
#     sig[pairw.res$p.adjusted <= 0.05] <-'*'
#     sig[pairw.res$p.adjusted <= 0.01] <-'**'
#     sig[pairw.res$p.adjusted <= 0.001] <-'***'
#     pairw.res <- data.frame(pairw.res[,1:7],sig)
#   }
#   class(pairw.res) <- c("pwadonis", "data.frame")
#   return(pairw.res)
# }

## PERMDISP--Homogeneity of multivariate dispersions. 
#  A test to confirm assumption of homogeneity of variance prior to PERMANOVA.
# Must be done for factors separately
# For fraction:
betadisp.A = betadisper(DistVar.A, ps.A.df$fraction)
# For worm trt:
betadisp.A = betadisper(DistVar.A, ps.A.df$trt)

betadisp.A
plot(betadisp.A)
anova(betadisp.A) # NS for ps.A.df$fraction

# Significant for WORM trt:
# Response: Distances
#             Df  Sum Sq   Mean Sq F value   Pr(>F)   
# Groups      2 0.03967 0.0198369  5.2142 0.006614 **
#   Residuals 132 0.50218 0.0038044                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1    

permutest(betadisp.A)
TukeyHSD(betadisp.A)

# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#                     diff          lwr        upr     p adj
# LowWorm-NoWorm   -0.02008271 -0.050906035 0.01074062 0.2736547
# HighWorm-NoWorm   0.02189562 -0.008927702 0.05271895 0.2151826
# HighWorm-LowWorm  0.04197833  0.011155007 0.07280166 0.0044522

boxplot(betadisp.A, xlab = "", las = 2, cex.axis = 0.8)


######################## #
#### Statistics, Lakeshore ####
######################## #
# First, run PERMANOVA
ps.L = subset_samples(ps.L, sample_data(ps.L)$fraction != "casting")
ps.L.veg = veganotu(ps.L)
ps.L.df = data.frame(sample_data(ps.L))
DistVar.L = vegdist(ps.L.veg, method = "bray")
adonis2(DistVar.L ~ trt*fraction, data = ps.L.df, method = "bray")

# adonis2(formula = DistVar.L ~ trt * fraction, data = ps.L.df, method = "bray")
#               Df SumOfSqs     R2      F   Pr(>F)    
# trt           1   2.3397 0.16306 16.7495  0.001 ***
# fraction      2   0.1635 0.01140  0.5854  0.922    
# trt:fraction  2   0.1118 0.00779  0.4002  0.999    
# Residual     84  11.7337 0.81775                   
# Total        89  14.3487 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Run the pairwise adonis function for pairwise PERMANOVA
pairwise.adonis(DistVar.L, ps.L.df$trt)

#               pairs   Df SumsOfSqs  F.Model       R2  p.value   p.adjusted sig
# 1 NoWorm vs HighWorm  1  2.339683 17.14479 0.1630589   0.001      0.001  **


## PERMDISP--Homogeneity of multivariate dispersions. 
#  A test to confirm assumption of homogeneity of variance prior to PERMANOVA.
# Must be done for factors separately
# For fraction:
betadisp.L = betadisper(DistVar.L, ps.L.df$fraction)
# For worm trt:
betadisp.L = betadisper(DistVar.L, ps.L.df$trt)

betadisp.L
plot(betadisp.L)
anova(betadisp.L) # NS for ps.L.df$fraction

#significant in trt
# Response: Distances
#             Df   Sum Sq   Mean Sq     F value   Pr(>F)   
# Groups      1   0.025712  0.0257116  8.4444     0.004631 **
# Residuals   88  0.267945  0.0030448                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

permutest(betadisp.L)
TukeyHSD(betadisp.L)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
# diff         lwr         upr     p adj
# HighWorm-NoWorm -0.0338044 -0.05692245 -0.01068635 0.0046313
boxplot(betadisp.L, xlab = "", las = 2, cex.axis = 0.8)



#### Option to add pH or soil moisture to PCoA

#read in metadata that has soil pH
metadata = read.csv(file ="sample_metadata_Exp3_with_pH_moisture.csv", header=TRUE, stringsAsFactors=TRUE)
#metadata


#merge metadata df with sample_data to put the metadata in the correct order, clean things up
df2 = merge(data.frame(sample_data(ps)),metadata,by="SampleID")
head(df2)
colnames(df2)
df3 = df2[,c(1:7, 14:17)]
colnames(df3)
colnames(df3) = c("SampleID", "sample", "fraction", "site", "plot", "core", "trt",
                  "MixComm", "Blank", "pH", "moisture")
colnames(df3)
head(df3)
row.names(df3) = df3$SampleID
sample_data(ps) = sample_data(df3)
str(sample_data(ps))
sample_data(ps)

#### PCoA, Arboretum ####
ps.A = subset_samples(ps,sample_data(ps)$site == "Arb")
ps.A = subset_samples(ps.A,sample_data(ps.A)$fraction == "fresh")

# Prune away zero's
ps.A = prune_taxa(taxa_sums(ps.A) > 0, ps.A)

# Hellinger transformation
ps.A = transform_sample_counts(ps.A, function(x) (x / sum(x))^0.5 )

# PCoA, Arboretum
ps.PCoA = ordinate(ps.A, method="PCoA", distance="bray")

# Soil pH
p.A = plot_ordination(ps.A, ps.PCoA, color = "pH", shape = "trt") +
  geom_point(size = 3)
p.A = p.A + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.A = p.A + scale_shape_manual(values = c(16, 15, 17), labels = c("Control", "Low", "High"),
                               name = expression(paste(italic("Amynthas"), " pressure")))
p.A = p.A + scale_color_gradient2(low = "orange", mid = "darkorange2", midpoint = 6, high = "darkorange4")
p.A # 800 x 225 for MS

# Soil moisture
p.A = plot_ordination(ps.A, ps.PCoA, color = "moisture", shape = "trt") +
  geom_point(size = 3)
p.A = p.A + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.A = p.A + scale_shape_manual(values = c(16, 15, 17), labels = c("Control", "Low", "High"),
                               name = expression(paste(italic("Amynthas"), " pressure")))
p.A = p.A + scale_color_gradient2(low = "lightblue1", mid = "skyblue4", high = "midnightblue", midpoint = 0.35, name = "Soil moisture")
p.A # 800 x 225 for MS

######################## #
#### PCoA, Lakeshore ####
######################## #

ps.L = subset_samples(ps,sample_data(ps)$site == "Lak")
ps.L = subset_samples(ps.L,sample_data(ps.L)$fraction == "fresh")

# Prune away zero's
ps.L = prune_taxa(taxa_sums(ps.L) > 0, ps.L)

# Do a Hellinger transformation
ps.L = transform_sample_counts(ps.L, function(x) (x / sum(x))^0.5 )

# PCoA, pH, Lakeshore
ps.PCoA = ordinate(ps.L, method="PCoA", distance="bray")
p.L = plot_ordination(ps.L, ps.PCoA, color = "pH", shape = "trt") +
  geom_point(size = 3)
p.L = p.L + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.L = p.L + scale_shape_manual(values = c(16, 17), labels = c("Control", "High"),
                               name = expression(paste(italic("Amynthas"), " pressure")))
p.L = p.L + scale_color_gradient2(low = "orange", mid = "darkorange2", midpoint = 6, high = "darkorange4")
p.L # 800x225


# PCoA, moisture, Lakeshore
ps.PCoA = ordinate(ps.L, method="PCoA", distance="bray")
p.L = plot_ordination(ps.L, ps.PCoA, color = "moisture", shape = "trt") +
  geom_point(size = 3)
p.L = p.L + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.L = p.L + scale_shape_manual(values = c(16, 17), labels = c("Control", "High"),
                               name = expression(paste(italic("Amynthas"), " pressure")))
p.L = p.L + scale_color_gradient2(low = "lightblue1", mid = "skyblue4", high = "midnightblue", midpoint = 0.35, name = "Soil moisture")
p.L


######################## #
#### Statistics, Gallistel ####
######################## #
# First, run PERMANOVA

# # Create veganotu function
# veganotu = function(physeq) {
#   require("vegan")
#   OTU = otu_table(physeq)
#   if (taxa_are_rows(OTU)) {
#     OTU = t(OTU)
#   }
#   return(as(OTU, "matrix"))
# }

ps.A.veg = veganotu(ps.A)
ps.A.df = data.frame(sample_data(ps.A))
DistVar.A = vegdist(ps.A.veg, method = "bray")
#DistVar.A
adonis2(DistVar.A ~ moisture * trt, data = ps.A.df, method = "bray", na.action = na.omit)


# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
#
# adonis2(formula = DistVar.A ~ pH * trt, data = ps.A.df, method = "bray", na.action = na.omit)
#             Df SumOfSqs      R2       F Pr(>F)    
#   pH        1   1.6301 0.22802 16.6127  0.001 ***
#   trt       2   1.1753 0.16440  5.9888  0.001 ***
#   pH:trt    2   0.6147 0.08599  3.1322  0.001 ***
#   Residual 38   3.7286 0.52159                   
#   Total    43   7.1487 1.00000                   
# ---
# adonis2(formula = DistVar.A ~ moisture * trt, data = ps.A.df, method = "bray", na.action = na.omit)
#                 Df SumOfSqs      R2      F Pr(>F)    
#   moisture      1   0.7450 0.10546 6.4549  0.001 ***
#   trt           2   1.2637 0.17888 5.4744  0.001 ***
#   moisture:trt  2   0.6700 0.09484 2.9026  0.001 ***
#   Residual     38   4.3860 0.62082                  
#   Total        43   7.0648 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


######################## #
#### Statistics, Lakeshore ####
######################## #
# First, run PERMANOVA
ps.L.veg = veganotu(ps.L)
ps.L.df = data.frame(sample_data(ps.L))
DistVar.L = vegdist(ps.L.veg, method = "bray")
adonis2(DistVar.L ~ moisture * trt, data = ps.L.df, method = "bray")

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = DistVar.L ~ pH * trt, data = ps.L.df, method = "bray")
#             Df SumOfSqs      R2      F Pr(>F)   
#   pH        1   0.3740 0.07854 2.9479  0.016 * 
#   trt       1   0.7884 0.16556 6.2146  0.002 **
#   pH:trt    1   0.3011 0.06323 2.3734  0.027 * 
#   Residual 26   3.2986 0.69267                 
#   Total    29   4.7621 1.00000  
# 
# adonis2(formula = DistVar.L ~ moisture * trt, data = ps.L.df, method = "bray")
#                  Df SumOfSqs      R2      F Pr(>F)    
#   moisture      1   0.4109 0.08628 3.6906  0.006 ** 
#   trt           1   0.7896 0.16580 7.0920  0.001 ***
#   moisture:trt  1   0.6670 0.14006 5.9908  0.001 ***
#   Residual     26   2.8947 0.60785                  
#   Total        29   4.7621 1.00000  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


