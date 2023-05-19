library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(rcartocolor)


TCTN = read.csv(file ="TC_TN_forR.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA")
head(TCTN)

# Subset for site(s)
TCTN = subset(TCTN, site=="Arb" | site == "Lak")

# Find mean
TCTNmean = TCTN %>%
  group_by(site, trt, fraction) %>% 
  dplyr::summarize(n=n(),
                   "Total carbon, percent" = mean(C_pct, na.rm=TRUE),
                   "Total nitrogen, percent" = mean(N_pct, na.rm=TRUE),
                   "C:N ratio" = mean(C_N_ratio, na.rm=TRUE)) %>%
  gather("element", "mean", - c(site, trt, fraction, n), factor_key=TRUE)  
TCTNmean

# Find SE
TCTNSE = TCTN %>%
  group_by(site, trt, fraction) %>% 
  dplyr::summarize(n=n(),
                   "Total carbon, percent" = sd(C_pct, na.rm=TRUE)/sqrt(n()),
                   "Total nitrogen, percent" = sd(N_pct, na.rm=TRUE)/sqrt(n()),
                   "C:N ratio" = sd(C_N_ratio, na.rm=TRUE)/sqrt(n())) %>%
  gather("element", "SE", - c(site, trt, fraction, n), factor_key=TRUE)  
TCTNSE

# Join mean proportion and SE data together
TCTNmeanSE = merge(TCTNmean,TCTNSE, by=c("site","trt", "fraction", "n", "element"))      
head(TCTNmeanSE)

TCTNmeanSE = TCTNmeanSE %>%
  mutate(occluded = case_when(fraction == "occm" ~ "occluded",
                              fraction == "fresh" | fraction == "freem" | fraction == "mac" ~ "not"))


# Rename and reorder fractions for graphing
TCTNmeanSE$trt = recode_factor(TCTNmeanSE$trt, NoWorm = "Control", 
                              LowWorm = "Low", HighWorm = "High")
TCTNmeanSE$trt = ordered(TCTNmeanSE$trt, levels=c("Control", "Low", "High"))

TCTNmeanSE$fraction = recode_factor(TCTNmeanSE$fraction, fresh = "Bulk \nsoil",
                                    mac = "Macroagregate", freem = "Free \nmicroaggregate",
                                    occm = "Occluded \nmicroaggregate")
TCTNmeanSE$fraction = ordered(TCTNmeanSE$fraction, levels=c("Bulk \nsoil", "Macroagregate",
                                                            "Free \nmicroaggregate",
                                                            "Occluded \nmicroaggregate"))

trtpalette = c("#fcb18c","#B7003F","#2d105a") #Worm sites
TCTNmeanSE$site = ordered(TCTNmeanSE$site, levels=c("Arb", "Lak"))
site.labs=c("Arb" = "Gallistel", "Lak" = "Lakeshore")
AmynthasTitle = expression(paste(italic("Amynthas"), " pressure"))

### Graph: Total C Percent, BAR graph
c = ggplot(subset(TCTNmeanSE, element == "Total carbon, percent"), aes(x = fraction, y = mean, fill=trt, pattern=occluded)) +
  geom_bar_pattern(stat = "identity", color="black", position=position_dodge(),
                   pattern_colour = "#666666",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.03) +
  guides(fill = guide_legend(override.aes=list(pattern="none"))) +
  scale_pattern_manual(values=c(occluded="circle", not="none"), guide = "none")
c = c + geom_errorbar(aes(ymin=mean-1.96*SE,ymax=mean+1.96*SE), position=position_dodge(0.9),width=0.15)
c = c + facet_grid(~site, labeller=labeller(site=site.labs))
c = c + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
c = c + ylab("Total carbon (%)")  + ylim(-1, 17.5)
c = c + xlab("Soil fraction")
c = c + scale_fill_manual(AmynthasTitle, values=trtpalette)
c = c + theme(legend.position = c(0.65, 0.78), legend.box = "horizontal")
c = c + theme(legend.key = element_rect(color = "black"))
c # 800 x 400 for SI


# ### Graph: Total N Percent, BAR graph
n = ggplot(subset(TCTNmeanSE, element == "Total nitrogen, percent"), aes(x = fraction, y = mean, fill=trt, pattern=occluded)) +
  geom_bar_pattern(stat = "identity", color="black", position=position_dodge(),
                   pattern_colour = "#666666",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.03) +
  guides(fill = guide_legend(override.aes=list(pattern="none"))) +
  scale_pattern_manual(values=c(occluded="circle", not="none"), guide = "none")
n = n + geom_errorbar(aes(ymin=mean-1.96*SE,ymax=mean+1.96*SE), position=position_dodge(0.9),width=0.15)
n = n + facet_grid(~site, labeller=labeller(site=site.labs))
n = n + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
n = n + ylab("Total nitrogen (%)")  + ylim(-0.05, 1)
n = n + xlab("Soil fraction")
n = n + scale_fill_manual(AmynthasTitle, values=trtpalette)
n = n + theme(legend.position = c(0.62, 0.81), legend.box = "horizontal")
n = n + theme(legend.key = element_rect(color = "black"))
n # 800 x 400 for SI

setwd("~/Box\ Sync/Research/0_Comm_Assemb_Exp_3/Soil analyses/TC TN")
ggsave("Npct.WORM.tiff", width=8, height=4, units = "in", device='tiff', dpi=400)



# Statistics, treatment differences within fraction
# Must update site for different sites
# and must update the fraction in the aov statement (e.g., mac_C, freem_C, occm_C, fresh)
dat = subset(TCTN, site == "Lak" & fraction == "occm")
ano = aov(N_pct ~ trt, data = dat)
summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

cld = multcompLetters4(ano, tuk)
cld



head(TCTN)
### Graph: boxplots, C:N ratio

# Order fractions for graphing
TCTN$trt = recode_factor(TCTN$trt, NoWorm = "Control", 
                               LowWorm = "Low", HighWorm = "High")
TCTN$trt = ordered(TCTN$trt, levels=c("Control", "Low", "High"))

TCTN$fraction = recode_factor(TCTN$fraction, fresh = "Bulk \nsoil",
                                    mac = "Macroagregate", freem = "Free \nmicroaggregate",
                                    occm = "Occluded \nmicroaggregate")
TCTN$fraction = ordered(TCTN$fraction, levels=c("Bulk \nsoil", "Macroagregate",
                                                            "Free \nmicroaggregate",
                                                            "Occluded \nmicroaggregate"))
trtpalette = c("#fcb18c","#B7003F","#2d105a") #Worm sites
TCTN$site = ordered(TCTN$site, levels=c("Arb", "Lak"))

p = ggplot(subset(TCTN, site=="Arb" | site=="Lak"), aes(x = fraction, y = C_N_ratio, color=trt)) +
   geom_boxplot()
p = p + facet_grid(~site, labeller=labeller(site=site.labs))
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab("C:N ratio") # + ylim(1250,2500)
p = p + xlab("Soil fraction")
p = p + scale_color_manual(AmynthasTitle, values=trtpalette)
p = p + theme(legend.position = c(0.1, 0.78), legend.box = "horizontal")
p # 800 x 350 for SI

# this may help to produce raw numbers for a table:
subset(TCTNmeanSE, element=="Total nitrogen, percent" & trt=="High")
TCTN

#### ANOVA, by site, for table ####
# Arlington, TC
# Total C was not significantly affected by Amynthas at Arb (occm and mac were p <0.1)
dat = subset(TCTN, site == "Lak")
ano = aov(N_pct ~ trt * fraction, data = dat)
summary(ano)
model.tables(ano, "means", se=TRUE)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

#library(multcompView)
cld = multcompLetters4(ano, tuk)
cld



#### Graph of total C per g of soil by fraction
# need to create agg mean

# Get carbon data
TCTN = read.csv(file ="TC_TN_forR.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA")
head(TCTN)
TCTN.C = TCTN[, -c(1,8,9,11)]
TCTN.C = spread(TCTN.C, fraction, C_pct)
head(TCTN.C)

# Get aggregate proportion data
agg = read.csv("aggregate_fractions_proportions.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA", row.names=1)
head(agg)
agg = agg[, -c(13,14)]

TCTN.C.agg = merge(TCTN.C,agg, by=c("sample", "site","trt"))      
head(TCTNagg)
TCTN.C.agg = subset(TCTN.C.agg, site == "Arb" | site == "Lak")
TCTN.C.agg = subset(TCTN.C.agg, sample != "235")

TCTN.C.agg$freem_C = TCTN.C.agg$freem.x * TCTN.C.agg$freem.y *10
TCTN.C.agg$occm_C = TCTN.C.agg$occm * TCTN.C.agg$occm.Mac* TCTN.C.agg$Mac *10
TCTN.C.agg$mac_C = TCTN.C.agg$mac *  TCTN.C.agg$Mac *10

# Option to add in estimates for silt + clay fractions
TCTN.C.agg$SC_C = (TCTN.C.agg$fresh*10-TCTN.C.agg$mac_C-TCTN.C.agg$freem_C)
TCTN.C.agg$occSC_C = (TCTN.C.agg$mac_C-TCTN.C.agg$occm_C)

levels(TCTNmeanSE$site)
# Find proportion of each aggregate
TCTN.C.aggmean = TCTN.C.agg %>%
  group_by(site, trt) %>% 
  dplyr::summarize(n=n(),
                   "freem" = mean(freem_C, na.rm=TRUE),
                   "mac" = mean(mac_C, na.rm = TRUE),
                   "occm" = mean(occm_C, na.rm = TRUE),
                   "fresh" = mean((fresh*10), na.rm = TRUE),
                   "SC" = mean(SC_C, na.rm=TRUE),
                   "occSC" = mean(occSC_C, na.rm=TRUE)) %>%
  gather("fraction", "agg_fraction_C", - c(site, trt, n), factor_key=TRUE)  
TCTN.C.aggmean
str(TCTN.C.aggmean)

# Find SE
TCTN.C.aggSE = subset(TCTN.C.agg) %>%
  group_by(site, trt) %>% 
  dplyr::summarize(n=n(),
                   "freem" = sd(freem_C, na.rm=TRUE)/sqrt(n()),
                   "mac" = sd(mac_C, na.rm = TRUE)/sqrt(n),
                   "occm" = sd(occm_C, na.rm=TRUE)/sqrt(n()),
                   "fresh" = sd((fresh*10), na.rm=TRUE)/sqrt(n()),
                   "SC" = sd(SC_C, na.rm=TRUE)/sqrt(n()),
                   "occSC" = sd(occSC_C, na.rm=TRUE)/sqrt(n())
  ) %>%
  gather("fraction", "SE", - c(site, trt, n), factor_key=TRUE)  
TCTN.C.aggSE
str(TCTNmean)

# Join mean proportion and SE data together
aggmeanSE = merge(TCTN.C.aggmean,TCTN.C.aggSE, by=c("site","trt", "n", "fraction"))      
str(aggmeanSE)

# Find the proportion of total C (in bulk soil) in each fraction
aggmeanSE = aggmeanSE %>%
  group_by(site, trt) %>%
  mutate(prop = agg_fraction_C / agg_fraction_C[match('fresh', fraction)])
aggmeanSE$prop.SE = aggmeanSE$SE*aggmeanSE$prop/aggmeanSE$agg_fraction_C

aggmeanSE = aggmeanSE %>%
  mutate(occluded = case_when(fraction == "occm" | fraction == "occSC" ~ "occluded",
                              fraction == "fresh" | fraction == "freem" | fraction == "mac" | fraction == "SC"~ "not"))
str(aggmeanSE)


#option to save
#write.csv(aggmeanSE,"Carbon_in_aggregate_fractions_WORM.csv")

#aggmeanSE = read.csv("Carbon_in_aggregate_fractions_WORM.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA", row.names=1)


### Graph: amount of C in each soil fraction--FOR MANUSCRIPT FIGURE
# Order fractions for graphing
aggmeanSE$fraction = ordered(aggmeanSE$fraction, levels=c("fresh", "mac", "freem", "SC", "occm", "occSC"))
aggmeanSE$trt = ordered(aggmeanSE$trt, levels=c("NoWorm", "LowWorm", "HighWorm"))
site.labs=c("Arb" = "Gallistel", "Lak" = "Lakeshore")
aggmeanSE$trt = recode_factor(aggmeanSE$trt, NoWorm = "Control",
                              LowWorm = "Low",
                              HighWorm  ="High")
AmynthasTitle = expression(paste(italic("Amynthas"), " pressure"))
trtpalette = c("#fcb18c","#B7003F","#2d105a") #Worm sites

p = ggplot(aggmeanSE, aes(x = fraction, y = agg_fraction_C, fill=trt, pattern = occluded)) +
  geom_bar_pattern(stat = "identity", color="black", position=position_dodge(),
                   pattern_colour = "#666666",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.03) +
  guides(fill = guide_legend(override.aes=list(pattern="none"))) +
  scale_pattern_manual(values=c(occluded="circle", not="none"), guide = "none") +
  geom_errorbar(data=subset(aggmeanSE, fraction!="SC" & fraction!="occSC"), #no error bars for the estimates of SC fraction
                aes(ymin=agg_fraction_C-1.96*SE,ymax=agg_fraction_C+1.96*SE), position=position_dodge(0.9),width=0.15)
p = p + facet_grid(~site, labeller=labeller(site=site.labs))
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab("mg C (in fraction) per g bulk soil") #+ ylim(-1,40)
p = p + xlab("Soil fraction")
p = p + scale_x_discrete(labels=c("fresh"="Bulk \nsoil","mac"= "Macro.", "freem"="Free \nmicro.", "SC"="Silt + clay \n(estimate)",
                                  "occm"="Occluded \nmicro.", "occSC"="Occluded \nsilt + clay \n& POM \n(estimate)"))
p = p + scale_fill_manual(AmynthasTitle, values=trtpalette)
p = p + theme(legend.position = c(0.9, 0.75), legend.text=element_text(size=rel(0.75)), legend.title = element_text(size=rel(0.85)))
p = p + theme(legend.key = element_rect(color = "black"))
p

ggsave("CinFraction.WORM.tiff", width=7.5, height=3.5, units = "in", device='tiff', dpi=400)

head(TCTN.C.agg)
# Statistics, treatment differences within fraction
# Must update site for different sites
# and must update the fraction in the aov statement (fresh, mac_C, freem_C, occm_C)
dat = subset(TCTN.C.agg, site == "Lak")
ano = aov(fresh ~ trt, data = dat)
summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

cld = multcompLetters4(ano, tuk)
cld



#### Graph of total N per g of soil by fraction
# need to create agg mean

# Get nitrogen data
TCTN = read.csv(file ="TC_TN_forR.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA")
head(TCTN)
TCTN.N = TCTN[, -c(1,8,10,11)]
TCTN.N = spread(TCTN.N, fraction, N_pct)
head(TCTN.N)

# Get aggregate proportion data
agg = read.csv("aggregate_fractions_proportions.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA", row.names=1)
head(agg)
agg = agg[, -c(13,14)]

TCTN.N.agg = merge(TCTN.N,agg, by=c("sample", "site","trt"))      
head(TCTNagg)
TCTN.N.agg = subset(TCTN.N.agg, site == "Arb" | site == "Lak")
TCTN.N.agg = subset(TCTN.N.agg, sample != "235")
head(TCTN.N.agg)



TCTN.N.agg$freem_N = TCTN.N.agg$freem.x * TCTN.N.agg$freem.y *10
TCTN.N.agg$occm_N = TCTN.N.agg$occm * TCTN.N.agg$occm.Mac* TCTN.N.agg$Mac *10
TCTN.N.agg$mac_N = TCTN.N.agg$mac *  TCTN.N.agg$Mac *10

# Option to add in estimates for silt + clay fractions
TCTN.N.agg$SC_N = (TCTN.N.agg$fresh*10-TCTN.N.agg$mac_N-TCTN.N.agg$freem_N)
TCTN.N.agg$occSC_N = (TCTN.N.agg$mac_N-TCTN.N.agg$occm_N)

head(TCTN.N.agg)

levels(TCTNmeanSE$site)
# Find proportion of each aggregate
TCTN.N.aggmean = TCTN.N.agg %>%
  group_by(site, trt) %>% 
  dplyr::summarize(n=n(),
                   "freem" = mean(freem_N, na.rm=TRUE),
                   "mac" = mean(mac_N, na.rm = TRUE),
                   "occm" = mean(occm_N, na.rm = TRUE),
                   "fresh" = mean((fresh*10), na.rm = TRUE),
                   "SC" = mean(SC_N, na.rm=TRUE),
                   "occSC" = mean(occSC_N, na.rm=TRUE)) %>%
  gather("fraction", "agg_fraction_N", - c(site, trt, n), factor_key=TRUE)  
head(TCTN.N.aggmean)
str(TCTN.N.aggmean)

# Find SE
TCTN.N.aggSE = subset(TCTN.N.agg) %>%
  group_by(site, trt) %>% 
  dplyr::summarize(n=n(),
                   "freem" = sd(freem_N, na.rm=TRUE)/sqrt(n()),
                   "mac" = sd(mac_N, na.rm = TRUE)/sqrt(n),
                   "occm" = sd(occm_N, na.rm=TRUE)/sqrt(n()),
                   "fresh" = sd((fresh*10), na.rm=TRUE)/sqrt(n()),
                   "SC" = sd(SC_N, na.rm=TRUE)/sqrt(n()),
                   "occSC" = sd(occSC_N, na.rm=TRUE)/sqrt(n())
  ) %>%
  gather("fraction", "SE", - c(site, trt, n), factor_key=TRUE)  
TCTN.N.aggSE
str(TCTNmean)

# Join mean proportion and SE data together
aggmeanSE = merge(TCTN.N.aggmean,TCTN.N.aggSE, by=c("site","trt", "n", "fraction"))      
str(aggmeanSE)

# Find the proportion of total N (in bulk soil) in each fraction
aggmeanSE = aggmeanSE %>%
  group_by(site, trt) %>%
  mutate(prop = agg_fraction_N / agg_fraction_N[match('fresh', fraction)])
aggmeanSE$prop.SE = aggmeanSE$SE*aggmeanSE$prop/aggmeanSE$agg_fraction_N

aggmeanSE = aggmeanSE %>%
  mutate(occluded = case_when(fraction == "occm" | fraction == "occSC" ~ "occluded",
                              fraction == "fresh" | fraction == "freem" | fraction == "mac" | fraction == "SC"~ "not"))
str(aggmeanSE)


#option to save
#write.csv(aggmeanSE,"Nitrogen_in_aggregate_fractions_WORM.csv")

#aggmeanSE = read.csv("Nitrogen_in_aggregate_fractions_WORM.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA", row.names=1)


### Graph: amount of N in each soil fraction--FOR MANUSCRIPT FIGURE
# Order fractions for graphing
aggmeanSE$fraction = ordered(aggmeanSE$fraction, levels=c("fresh", "mac", "freem", "SC", "occm", "occSC"))
aggmeanSE$trt = ordered(aggmeanSE$trt, levels=c("NoWorm", "LowWorm", "HighWorm"))
site.labs=c("Arb" = "Gallistel", "Lak" = "Lakeshore")
aggmeanSE$trt = recode_factor(aggmeanSE$trt, NoWorm = "Control",
                              LowWorm = "Low",
                              HighWorm  ="High")
AmynthasTitle = expression(paste(italic("Amynthas"), " pressure"))
trtpalette = c("#fcb18c","#B7003F","#2d105a") #Worm sites

p = ggplot(aggmeanSE, aes(x = fraction, y = agg_fraction_N, fill=trt, pattern = occluded)) +
  geom_bar_pattern(stat = "identity", color="black", position=position_dodge(),
                   pattern_colour = "#666666",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.03) +
  guides(fill = guide_legend(override.aes=list(pattern="none"))) +
  scale_pattern_manual(values=c(occluded="circle", not="none"), guide = "none") +
  geom_errorbar(data=subset(aggmeanSE, fraction!="SC" & fraction!="occSC"), #no error bars for the estimates of SC fraction
                aes(ymin=agg_fraction_N-1.96*SE,ymax=agg_fraction_N+1.96*SE), position=position_dodge(0.9),width=0.15)
p = p + facet_grid(~site, labeller=labeller(site=site.labs))
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab("mg N (in fraction) per g bulk soil") #+ ylim(-1,40)
p = p + xlab("Soil fraction")
p = p + scale_x_discrete(labels=c("fresh"="Bulk \nsoil","mac"= "Macro.", "freem"="Free \nmicro.", "SC"="Silt + clay \n(estimate)",
                                  "occm"="Occluded \nmicro.", "occSC"="Occluded \nsilt + clay \n& POM \n(estimate)"))
p = p + scale_fill_manual(AmynthasTitle, values=trtpalette)
p = p + theme(legend.position = c(0.9, 0.75), legend.text=element_text(size=rel(0.75)), legend.title = element_text(size=rel(0.85)))
p = p + theme(legend.key = element_rect(color = "black"))
p

ggsave("N_in_Fraction.WORM.tiff", width=7.5, height=3.5, units = "in", device='tiff', dpi=400)

head(TCTN.N.agg)
# Statistics, treatment differences within fraction
# Must update site for different sites
# and must update the fraction in the aov statement (mac_C, freem_C, occm_C, fresh)
dat = subset(TCTN.N.agg, site == "Arb")
ano = aov(mac_N ~ trt, data = dat)
summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

cld = multcompLetters4(ano, tuk)
cld



### Graph: proportion of C in fraction
head(aggmeanSE)
p = ggplot(aggmeanSE, aes(x = fraction, y = prop, fill=trt, pattern = occluded)) +
  geom_bar_pattern(stat = "identity", color="black", position=position_dodge(),
                   pattern_colour = "#666666",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.05,
                   pattern_spacing = 0.03) +
  guides(fill = guide_legend(override.aes=list(pattern="none"))) +
  scale_pattern_manual(values=c(occluded="circle", not="none"), guide = "none") +  geom_errorbar(data=subset(aggmeanSE, fraction!="SC" & fraction!="occSC"), #no error bars for the estimates of SC fraction
                aes(ymin=prop-1.96*prop.SE,ymax=prop+1.96*prop.SE), position=position_dodge(0.9),width=0.15)
p = p + facet_grid(~site, labeller=labeller(site=site.labs))
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab("Proportion of total N in soil fraction") #+ ylim(-1,40)
p = p + xlab("Soil fraction")
p = p + scale_x_discrete(labels=c("fresh"="Bulk \nsoil","mac"= "Macro.", "freem"="Free \nmicro.", "SC"="Silt + clay \n(estimate)",
                                  "occm"="Occluded \nmicro.", "occSC"="Occluded \nsilt + clay \n& POM \n(estimate)"))
p = p + scale_fill_manual(AmynthasTitle, values=trtpalette)
p = p + theme(legend.position = c(0.9, 0.75), legend.text=element_text(size=rel(0.75)), legend.title = element_text(size=rel(0.85)))
p = p + theme(legend.key = element_rect(color = "black"))
#p = p + theme(axis.text.x = element_text(angle = 45))
p # 750 x 350 for SI


ggsave("NinFraction.Proportion of whole.WORM.tiff", width=7.5, height=3.5, units = "in", device='tiff', dpi=400)

