# Load required packages
library(ggplot2)
library(phyloseq)
library(dplyr)
library(picante)
library(DESeq2)


# Read in phyloseq object
ps = readRDS("ps.Exp3")
ps

### Prune to just one site
ps.site = subset_samples(ps,sample_data(ps)$site == "Lak")

# Prune away zero's
ps.site = prune_taxa(taxa_sums(ps.site) > 0, ps.site)
ps.site


# FYI: For Faith's PD, you do not need to normalize to relative abundances as this metric is based on presence/absence

# # That said, I did check how rarefied OTU table works--results are virtually identical.
# # # Set OTU table for rarefaction
# otu_data = otu_table(ps.site)
# otu_data[1:5, 1:5]
# # 
# # rarefied = rrarefy(otu_data, min(rowSums(otu_data)))
# # rarefied[1:5, 1:5]
# # otu_data = t(rarefied)
# # head(otu_data)
# 

# read in the phylogeny
# Rooted and unrooted trees seem to produce nearly identical results. Must modify (include.root) option.
phylo = ggtree::read.tree("tree_Exp3_fewer_quotes.nwk") # unrooted tree
#phylo = ggtree::read.tree("tree_rooted_Exp3.nwk") # rooted tree

### Notes on rooted vs unrooted trees from pd R documentation:
# If the root is to be included in all calculations of PD (include.root=TRUE), 
# the TREE MUST BE ROOTED. Single-species samples will be assigned a PD value equal
# to the distance from the root to the present.
#
# If the root is not included in all calculations by default (include.root=FALSE), 
# the tree need not rooted, but in the case of single-species samples the PD will 
# be equal to NA and a warning will be issued.

# ## Fix tip labels, IF they contain extra set of quotes...this seemed to be an issue with the unrooted tree!
# library(stringr)
# library(phylotools)
# tips = phylo$tip.label
# tips[1]
# new.tips =str_sub(tips, 2, -2)
# new.tips[1]
# df = data.frame(tips, new.tips)
# phylo2 = sub.taxa.label(phylo, df) #this step takes a while...
# write.tree(phylo2, "tree_Exp3_fewer_quotes.nwk")

head(phylo)
phylo.site = prune.sample(otu_table(ps.site), phylo)
phy_tree(ps.site)=phylo.site

# Run Faith's phylogenetic diversity for this site
Faiths = pd(otu_table(ps.site), phylo.site, include.root=FALSE)
Faiths
Faiths$SampleID = row.names(Faiths)

# Add sample metadata
meta = read.csv("sample_metadata_Exp3.csv", header=TRUE)
meta
meta = meta [,-c(8,9)]
head(meta)

Faiths = merge(Faiths, meta, by = "SampleID")

# ANOVA

ano = aov(PD ~ trt*fraction, data = Faiths)
summary(ano)
model.tables(ano, "means", se=TRUE)

par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

### Arboretum
#               Df Sum Sq Mean Sq F value   Pr(>F)    
# trt            2   1982   991.0   8.214 0.000448 ***
# fraction       2    223   111.7   0.925 0.399086    
# trt:fraction   4    271    67.8   0.562 0.690827    
# Residuals    123  14840   120.6                                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



### Lakeshore 
#               Df Sum Sq Mean Sq F value Pr(>F)    
# trt           1   1004  1004.0  13.603 0.0004 ***
# fraction      2     38    19.0   0.257 0.7740    
# trt:fraction  2     52    26.2   0.354 0.7027    
# Residuals    84   6200    73.8      

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

### Arboretum
# $trt
#                     diff        lwr       upr     p adj
# LowWorm-HighWorm  9.217511  3.7238869 14.711134 0.0003417
# NoWorm-HighWorm   6.175635  0.5847716 11.766498 0.0265469
# NoWorm-LowWorm   -3.041876 -8.6327390  2.548988 0.4030598

### Lakeshore
# $trt
#                   diff      lwr      upr     p adj
# NoWorm-HighWorm 6.680089 3.078389 10.28179 0.0003996


### Option to save 
# write.csv(Faiths,"Derived_data/FaithsPD_WORM_Arb.csv")
# write.csv(Faiths,"Derived_data/FaithsPD_WORM_Lak.csv")


#### Load data, if not continuing from above ####
Faiths = read.csv(file="Derived_data/FaithsPD_WORM_Lak.csv", row.names = 1)
# head(Faiths)

#### Plot Estimates ####

## Very basic plot
p = ggplot(Faiths,aes(y=PD,x=fraction,color=trt))
p = p + geom_point()
p


## Add SD bars
PDplot = Faiths %>%
    group_by(trt, fraction) %>%
    dplyr::summarize(PD_mean=mean(PD),
                     PD_sd=sd(PD))

p = ggplot(PDplot,aes(y=PD_mean,x=fraction,color=trt))
p = p + geom_point() + geom_errorbar(aes(ymin = PD_mean - PD_sd, ymax = PD_mean + PD_sd))
#p = p + ylim(2000,4500)
p

### Create a combined 'Faiths', with multiple sites, for facet graphing, below.
Faiths.1 = read.csv(file="Derived_data/FaithsPD_WORM_Arb.csv", row.names = 1)
Faiths.2 = read.csv(file="Derived_data/FaithsPD_WORM_Lak.csv", row.names = 1)

Faiths = rbind(Faiths.1, Faiths.2)

# Find mean and error, for graphing
PDplot = Faiths %>%
  group_by(site, trt, fraction) %>%
  dplyr::summarize(PD_mean=mean(PD),
                   PD_se=sd(PD)/sqrt(n()))
head(PDplot)

# Re-name variables and put them in order
PDplot$trt = recode_factor(PDplot$trt, "NoWorm" = "Control", "LowWorm" = "Low","HighWorm" = "High")
PDplot$trt = ordered(PDplot$trt, levels = c("Control", "Low", "High"))
PDplot$fraction = recode_factor(PDplot$fraction,
                                      "fresh" = "Bulk \nsoil",
                                      "freem" = "Free \nmicroagg.",
                                      "occm" = "Occluded \nmicroagg.")
PDplot$fraction = ordered(PDplot$fraction, levels=c("Bulk \nsoil", "Free \nmicroagg.", "Occluded \nmicroagg."))
PDplot$site = recode_factor(PDplot$site,Arb = "Gallistel", Lak = "Lakeshore")
PDplot$site = ordered(PDplot$site, levels=c("Gallistel", "Lakeshore"))

trtpalette = c("#fcb18c","#B7003F","#2d105a") #Worm sites
AmynthasTitle = expression(paste(italic("Amynthas"), " pressure"))
dodge = position_dodge(0.5)

# Plot final richness estimates with errors: ±1.96SE should represent 95% confidence intervals
p = ggplot(PDplot,aes(y=PD_mean,x=fraction,color=trt))
p = p + geom_point(size=2, position = dodge) + geom_errorbar(aes(ymin=PD_mean-1.96*PD_se, ymax=PD_mean+1.96*PD_se), position=dodge, width=0.15)
p = p + facet_grid(~site)
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylim(90,130) + ylab("Faith's phylogenetic diversity")
p = p + xlab("Soil fraction")
p = p + scale_color_manual(AmynthasTitle, values=trtpalette)
#p = p + theme(legend.title = element_blank()) #+ theme(legend.position = c(0.15, 0.15), legend.box = "horizontal")
#p = p + theme(axis.text.x = element_text(angle = 45))
p # 650 x 350

#ggsave("Figures/PD.WORM.tiff", width=6.5, height=3.5, units = "in", device='tiff', dpi=400)


# Plot final estimates across treatments only!
head(Faiths)
# Re-name variables and put them in order
Faiths$trt = recode_factor(Faiths$trt, "NoWorm" = "Control", "LowWorm" = "Low","HighWorm" = "High")
Faiths$trt = ordered(Faiths$trt, levels = c("Control", "Low", "High"))
Faiths$site = recode_factor(Faiths$site,Arb = "Gallistel", Lak = "Lakeshore")
Faiths$site = ordered(Faiths$site, levels=c("Gallistel", "Lakeshore"))

p = ggplot(Faiths,aes(y=PD, x=trt, color=trt)) +
  geom_boxplot()
p = p + facet_grid(~site)
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylim(50,150) + ylab("Faith's phylogenetic diversity")
p = p + xlab(AmynthasTitle)
p = p + scale_color_manual(AmynthasTitle, values=trtpalette)
#p = p + theme(legend.title = element_blank()) #+ theme(legend.position = c(0.15, 0.15), legend.box = "horizontal")
#p = p + theme(axis.text.x = element_text(angle = 45))
p # 650 x 350

ggsave("Figures/PD.WORM_bytrtonly.tiff", width=6.5, height=3.5, units = "in", device='tiff', dpi=400)
