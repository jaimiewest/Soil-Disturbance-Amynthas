library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)
library(ggpubr)

# Read in phyloseq object
ps = readRDS("ps.Exp3")

######################## #
#### Subset for Galistel/Arboretum ####
######################## #
ps.A = subset_samples(ps,sample_data(ps)$site == "Arb")

# Prune away zero's
ps.A = prune_taxa(taxa_sums(ps.A) > 0, ps.A)
ps.A
sample_data(ps.A)

# Hellinger transformation
ps.A = transform_sample_counts(ps.A, function(x) (x / sum(x))^0.5 )


######################## #
#### Subset for Lakeshore ####
######################## #

ps.L = subset_samples(ps,sample_data(ps)$site == "Lak")

# Prune away zero's
ps.L = prune_taxa(taxa_sums(ps.L) > 0, ps.L)
ps.L
sample_data(ps.L)

# Hellinger transformation
ps.L = transform_sample_counts(ps.L, function(x) (x / sum(x))^0.5 )


######################## #
#### Distance to Spatial Median/Centroid, Gallistel ####
######################## #

# Must do this for each fraction separately--freem, occm, and fresh, as follows.

# ## Create veganotu function
# veganotu = function(physeq) {
#   require("vegan")
#   OTU = otu_table(physeq)
#   if (taxa_are_rows(OTU)) {
#     OTU = t(OTU)
#   }
#   return(as(OTU, "matrix"))
# }


### Free microaggregates, Gallistel
ps.A.fr = subset_samples(ps.A, fraction == "freem")
DistVar.A.fr = vegdist(veganotu(ps.A.fr), method = "bray")
ps.A.fr.df = data.frame(sample_data(ps.A.fr))

betadisp.A.fr = betadisper(DistVar.A.fr, ps.A.fr.df$trt)

trt = data.frame(betadisp.A.fr$group)
distances = data.frame(betadisp.A.fr$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Arb"
data$fraction= "freem"
data.Arb.fr = data

### Occluded microaggregates, Gallistel
ps.A.oc = subset_samples(ps.A, fraction == "occm")
DistVar.A.oc = vegdist(veganotu(ps.A.oc), method = "bray")
ps.A.oc.df = data.frame(sample_data(ps.A.oc))

betadisp.A.oc = betadisper(DistVar.A.oc, ps.A.oc.df$trt)

trt = data.frame(betadisp.A.oc$group)
distances = data.frame(betadisp.A.oc$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Arb"
data$fraction= "occm"
data.Arb.oc = data

### Bulk soil ("fresh"), Gallistel
ps.A.bu = subset_samples(ps.A, fraction == "fresh")
DistVar.A.bu = vegdist(veganotu(ps.A.bu), method = "bray")
ps.A.bu.df = data.frame(sample_data(ps.A.bu))

betadisp.A.bu = betadisper(DistVar.A.bu, ps.A.bu.df$trt)

trt = data.frame(betadisp.A.bu$group)
distances = data.frame(betadisp.A.bu$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Arb"
data$fraction= "fresh"
data.Arb.bu = data


######################## #
#### Distance to Centroid, Lakeshore ####
######################## #

### Free microaggregates, Lakeshore
ps.L.fr = subset_samples(ps.L, fraction == "freem")
DistVar.L.fr = vegdist(veganotu(ps.L.fr), method = "bray")
ps.L.fr.df = data.frame(sample_data(ps.L.fr))

betadisp.L.fr = betadisper(DistVar.L.fr, ps.L.fr.df$trt)

trt = data.frame(betadisp.L.fr$group)
distances = data.frame(betadisp.L.fr$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Lak"
data$fraction= "freem"
data.Lak.fr = data

### Occluded microaggregates, Lakeshore
ps.L.oc = subset_samples(ps.L, fraction == "occm")
DistVar.L.oc = vegdist(veganotu(ps.L.oc), method = "bray")
ps.L.oc.df = data.frame(sample_data(ps.L.oc))

betadisp.L.oc = betadisper(DistVar.L.oc, ps.L.oc.df$trt)

trt = data.frame(betadisp.L.oc$group)
distances = data.frame(betadisp.L.oc$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Lak"
data$fraction= "occm"
data.Lak.oc = data

### Bulk soil, Lakeshore
ps.L.bu = subset_samples(ps.L, fraction == "fresh")
DistVar.L.bu = vegdist(veganotu(ps.L.bu), method = "bray")
ps.L.bu.df = data.frame(sample_data(ps.L.bu))

betadisp.L.bu = betadisper(DistVar.L.bu, ps.L.bu.df$trt)

trt = data.frame(betadisp.L.bu$group)
distances = data.frame(betadisp.L.bu$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Lak"
data$fraction= "fresh"
data.Lak.bu = data


################################## #
#### Boxplots of Distance to Centroid
################################## #

# Combine the dataframes
data.bytrt = rbind(data.Arb.fr, data.Lak.fr, data.Arb.oc, data.Lak.oc, data.Arb.bu, data.Lak.bu)
data.bytrt

# Re-name variables and put them in order
data.bytrt$trt = recode_factor(data.bytrt$trt, "NoWorm" = "Control", "LowWorm" = "Low","HighWorm" = "High")
data.bytrt$trt = ordered(data.bytrt$trt, levels = c("Control", "Low", "High"))
data.bytrt$fraction = recode_factor(data.bytrt$fraction,
                                "fresh" = "Bulk \nsoil",
                                "freem" = "Free \nmicroagg.",
                                "occm" = "Occluded \nmicroagg.")
data.bytrt$fraction = ordered(data.bytrt$fraction, levels=c("Bulk \nsoil", "Free \nmicroagg.", "Occluded \nmicroagg."))
data.bytrt$site = ordered(data.bytrt$site, levels=c("Arb", "Lak"))

trtpalette = c("#fcb18c","#B7003F","#2d105a") # Gallistel site
trtpalette = c("#fcb18c","#2d105a") # Lakeshore site
AmynthasTitle = expression(paste(italic("Amynthas"), " pressure"))

disp.trt = ggplot(subset(data.bytrt, site =="Arb"), aes(x = trt, y = distances, color = trt ))
disp.trt = disp.trt + geom_boxplot()
disp.trt = disp.trt + facet_wrap(~ site, labeller = as_labeller(c("Arb" = "Gallistel, between-plot scale \n(i.e., treatment scale)",
                                                                  "Lak" = "Lakeshore, between-plot scale \n(i.e., treatment scale)"))) +
  labs(x = "", y = "Distance to spatial median", fill = "") +
  scale_color_manual(AmynthasTitle, values = trtpalette) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background = element_rect(fill="white")) +
  ylim(0.132, 0.575)
disp.trt # 550 x 370

# Re-run ggplot above separately for each site, and save


### Use ANOVA to find significant differences in distance to centroid, by site
dat.Arb = subset(data.bytrt, site =="Arb")
dat.Lak = subset(data.bytrt, site =="Lak")

ano = aov(distances ~ trt*fraction, data = dat.Lak)
summary(ano)
model.tables(ano, "means")
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))
tuk = TukeyHSD(ano, conf.level = 0.95)
tuk
cld <- multcompLetters4(ano, tuk)
cld
