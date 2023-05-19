library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)

CO2 = read.csv(file ="MicroResp_CO2_forR.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA")
row.names(CO2) = CO2$sample

head(CO2)

# Choose sites (or other factors) of interest
CO2=subset(CO2, site == "Arb" | site == "Lak")

# Find mean
CO2mean = CO2%>%
  group_by(site, trt) %>% 
  dplyr::summarize(n=n(),
                   "CO2 rate, per g soil" = mean(mean_CO2rate_ugC_perg_perh, na.rm=TRUE),
                   "CO2 rate, per g soil C" = mean(mean_CO2rate_ugC_per_gC_perh, na.rm=TRUE)) %>%
  gather("element", "mean", - c(site, trt, n), factor_key=TRUE)  
CO2mean

# Find SE
CO2SE = CO2%>%
  group_by(site, trt) %>% 
  dplyr::summarize(n=n(),
                   "CO2 rate, per g soil" = sd(mean_CO2rate_ugC_perg_perh, na.rm=TRUE)/sqrt(n()),
                   "CO2 rate, per g soil C" = sd(mean_CO2rate_ugC_per_gC_perh, na.rm=TRUE)/sqrt(n())) %>%
  gather("element", "SE", - c(site, trt, n), factor_key=TRUE)  
CO2SE

# Join mean proportion and SE data together
CO2meanSE = merge(CO2mean,CO2SE, by=c("site","trt", "n", "element"))      
CO2meanSE


# Re-name and order factors for graphing
CO2$trt = recode_factor(CO2$trt, NoWorm = "Control",
                              LowWorm = "Low",
                              HighWorm  ="High")
CO2$trt = ordered(CO2$trt, levels = c("Control", "Low", "High"))
CO2$site = ordered(CO2$site, levels=c("Arb", "Lak"))
site.labs=c("Arb" = "Gallistel", "Lak" = "Lakeshore")
trtpalette = c("#fcb18c","#B7003F","#2d105a") #Worm sites
AmynthasTitle = expression(paste(italic("Amynthas"), " pressure"))

### Graph: CO2 respiration rate--BOXPLOTS -- Per g soil basis
p = ggplot(CO2, aes(x=trt, y=mean_CO2rate_ugC_perg_perh, color=trt)) +
  geom_boxplot() #+
  #geom_jitter(alpha=0.3)
p = p + facet_grid(~site, labeller=labeller(site=site.labs), scale="free", space="free_x")
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab(expression(CO[2]*", \u03BCg C/g soil/hour")) # + ylim(1250,2500)
p = p + xlab(element_blank())
p = p + scale_color_manual(values=trtpalette)
p = p + theme(legend.title = element_blank()) + theme(legend.position = "none")#c(0.15, 0.15), legend.box = "horizontal")
p = p + scale_x_discrete(labels=c("Control"="Control","Low"=expression(paste("Low ", italic("Amynthas"))), "High"=expression(paste("High ", italic("Amynthas")))))
p = p + theme(legend.title= element_blank())
p # 650 x 330 for MS. wider to accommodate x axis labels

ggsave("CO2.pergsoil.worm_NOjitter.tiff", width=6.5, height=3.3, units = "in", device='tiff', dpi=400)


### Graph: CO2 respiration rate--BOXPLOTS -- Per g soil C basis
c = ggplot(CO2, aes(x=trt, y=mean_CO2rate_ugC_per_gC_perh, color=trt)) +
  geom_boxplot() #+
  #geom_jitter(alpha=0.3)
c = c + facet_grid(~site, labeller=labeller(site=site.labs), scale ="free", space="free_x")
c = c + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
c = c + ylab(expression(CO[2]*", \u03BCg C/g soil carbon/hour")) # + ylim(1250,2500)
c = c + xlab(element_blank())
c = c + scale_color_manual(values=trtpalette)
c = c + theme(legend.title = element_blank()) + theme(legend.position = "none")#c(0.15, 0.15), legend.box = "horizontal")
c = c + scale_x_discrete(labels=c("Control"="Control","Low"=expression(paste("Low ", italic("Amynthas"))), "High"=expression(paste("High ", italic("Amynthas")))))
c = c + theme(legend.title= element_blank())
c # 650 x 330 for MS. wider to accommodate x axis labels

ggsave("CO2.pergC.worm_NOjitter.tiff", width=6.5, height=3.3, units = "in", device='tiff', dpi=400)




#### ANOVA, by site ####
head(CO2)
dat = subset(CO2, site == "Lak")
ano = aov(mean_CO2rate_ugC_per_gC_perh ~ trt, data = dat)
ano = aov(mean_CO2rate_ugC_perg_perh ~ trt, data = dat)
summary(ano)
model.tables(ano, "means", se=TRUE)

par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

cld = multcompLetters4(ano, tuk)
cld



### Graph: CO2 respiration rate--by moisture
### on a "per g Carbon" basis
resp = ggplot(CO2, aes(x = soil_moisture*100, y = mean_CO2rate_ugC_per_gC_perh, color = trt))
resp = resp + geom_point() +
  geom_smooth(method = "lm", formula = y~x, se = FALSE) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~','~")),
           show.legend=FALSE)
resp = resp + facet_grid(~site, labeller = labeller(site =site.labs))
resp = resp + scale_color_manual(values = trtpalette, name = AmynthasTitle)
resp = resp + labs(x = "Soil moisture, percent")
resp = resp + ylab(expression(CO[2]*", \u03BCg C/g soil C/hour")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
resp

#ggsave("Respiration_vs_moisture.tiff", width=6.8, height=3.6, units = "in", device='tiff', dpi=400)

### on a "per g SOIL" basis
resp2 = ggplot(CO2, aes(x=soil_moisture*100, y=mean_CO2rate_ugC_perg_perh, color=trt))
resp2 = resp2 + geom_point() +
  geom_smooth(method = "lm", formula = y~x, se = FALSE) +
  stat_cor(aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~','~")),
           show.legend=FALSE)
resp2 = resp2 + facet_grid(~site, labeller = labeller(site =site.labs))
resp2 = resp2 + scale_color_manual(values = trtpalette, name = "Treatment",
                                 labels =c("Control" = expression(paste("Control")),
                                           "Low" = expression(paste("Low ", italic("Amynthas"))),
                                           "High" = expression(paste("High ", italic("Amynthas")))))
resp2 = resp2 + labs(x = "Soil moisture, percent")
resp2 = resp2 + ylab(expression(CO[2]*", \u03BCg C/g soil/hour")) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
resp2

