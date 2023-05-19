library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)
library(ggpubr)


# Read in phyloseq object
ps = readRDS("ps.Exp3")

######################## #
#### Subset for Gallistel ####
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
#### Distance to Centroid, WITHIN EACH PLOT, Gallistel ####
######################## #

# Must do this for each soil fraction--freem, occm, and fresh, as follows

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

betadisp.A.fr = betadisper(DistVar.A.fr, ps.A.fr.df$plot)
betadisp.A.fr
plot = data.frame(betadisp.A.fr$group)
distances = data.frame(betadisp.A.fr$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Arb"
data$fraction= "freem"
data.Arb.fr = data
data.Arb.fr

### Occluded microaggregates, Gallistel
ps.A.oc = subset_samples(ps.A, fraction == "occm")
DistVar.A.oc = vegdist(veganotu(ps.A.oc), method = "bray")
ps.A.oc.df = data.frame(sample_data(ps.A.oc))

betadisp.A.oc = betadisper(DistVar.A.oc, ps.A.oc.df$plot)

plot = data.frame(betadisp.A.oc$group)
distances = data.frame(betadisp.A.oc$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Arb"
data$fraction= "occm"
data.Arb.oc = data

### Bulk soil, Gallistel
ps.A.bu = subset_samples(ps.A, fraction == "fresh")
DistVar.A.bu = vegdist(veganotu(ps.A.bu), method = "bray")
ps.A.bu.df = data.frame(sample_data(ps.A.bu))

betadisp.A.bu = betadisper(DistVar.A.bu, ps.A.bu.df$plot)

plot = data.frame(betadisp.A.bu$group)
distances = data.frame(betadisp.A.bu$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Arb"
data$fraction= "fresh"
data.Arb.bu = data


######################## #
#### Distance to Centroid, WITHIN PLOT, Lakeshore ####
######################## #


### Free microaggregates, Lakeshore
ps.L.fr = subset_samples(ps.L, fraction == "freem")
DistVar.L.fr = vegdist(veganotu(ps.L.fr), method = "bray")
ps.L.fr.df = data.frame(sample_data(ps.L.fr))

betadisp.L.fr = betadisper(DistVar.L.fr, ps.L.fr.df$plot)

plot = data.frame(betadisp.L.fr$group)
distances = data.frame(betadisp.L.fr$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Lak"
data$fraction= "freem"
data.Lak.fr = data

### Occluded microaggregates, Lakeshore
ps.L.oc = subset_samples(ps.L, fraction == "occm")
DistVar.L.oc = vegdist(veganotu(ps.L.oc), method = "bray")
ps.L.oc.df = data.frame(sample_data(ps.L.oc))

betadisp.L.oc = betadisper(DistVar.L.oc, ps.L.oc.df$plot)

plot = data.frame(betadisp.L.oc$group)
distances = data.frame(betadisp.L.oc$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Lak"
data$fraction= "occm"
data.Lak.oc = data

### Bulk soil, Lakeshore
ps.L.bu = subset_samples(ps.L, fraction == "fresh")
DistVar.L.bu = vegdist(veganotu(ps.L.bu), method = "bray")
ps.L.bu.df = data.frame(sample_data(ps.L.bu))

betadisp.L.bu = betadisper(DistVar.L.bu, ps.L.bu.df$plot)

plot = data.frame(betadisp.L.bu$group)
distances = data.frame(betadisp.L.bu$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Lak"
data$fraction= "fresh"
data.Lak.bu = data


################################## #
#### Boxplots of Distance to Centroid, WITHIN PLOT
################################## #

# Combine the dataframes
data.byplot = rbind(data.Arb.fr, data.Lak.fr, data.Arb.oc, data.Lak.oc, data.Arb.bu, data.Lak.bu)
head(data.byplot)
data.byplot$SampleID = rownames(data.byplot)

# Add Amynthas pressure treatment, based on plot numbers
meta = read.csv("sample_metadata_Exp3.csv", header=TRUE)
head(meta)
meta = meta [,-c(3,4,5,8,9)]
head(meta)

data.byplot = merge(data.byplot, meta, by = "SampleID")
head(data.byplot)

# Re-name variables and put them in order
data.byplot$trt = recode_factor(data.byplot$trt, "NoWorm" = "Control", "LowWorm" = "Low","HighWorm" = "High")
data.byplot$trt = ordered(data.byplot$trt, levels = c("Control", "Low", "High"))
data.byplot$fraction = recode_factor(data.byplot$fraction,
                                    "fresh" = "Bulk \nsoil",
                                    "freem" = "Free \nmicroagg.",
                                    "occm" = "Occluded \nmicroagg.")
data.byplot$fraction = ordered(data.byplot$fraction, levels=c("Bulk \nsoil", "Free \nmicroagg.", "Occluded \nmicroagg."))
data.byplot$site = ordered(data.byplot$site, levels=c("Arb", "Lak"))

trtpalette = c("#fcb18c","#B7003F","#2d105a") # Gallistel site
trtpalette = c("#fcb18c","#2d105a") # Lakeshore site
AmynthasTitle = expression(paste(italic("Amynthas"), " pressure"))

disp.plot = ggplot(subset(data.byplot, site=="Lak"), aes(x = trt, y = distances, color = trt ))
disp.plot = disp.plot + geom_boxplot()
disp.plot = disp.plot + facet_wrap(~ site, labeller = as_labeller(c("Arb" = "Gallistel, within-plot scale \n ",
                                                                    "Lak" = "Lakeshore, within-plot scale \n "))) +
  labs(x = "", y = " ", fill = "") +
  scale_color_manual(AmynthasTitle, values = trtpalette) +
  scale_fill_manual(name=NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(strip.background = element_rect(fill="white")) +
  ylim(0.132,0.575)
disp.plot # 550 x 370

#disp.plot.Arb = disp.plot
#disp.plot.Lak = disp.plot

### Use ANOVA to find significant differences in distance to MEDIAN, by site
dat.Arb = subset(data.byplot, site =="Arb")
dat.Lak = subset(data.byplot, site =="Lak")

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
