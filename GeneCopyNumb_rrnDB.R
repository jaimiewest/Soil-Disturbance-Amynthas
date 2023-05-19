library(Biostrings)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(multcompView)



##### We want to see how the mean predicted 16S copy number corresponds to Amynthas worm pressure

# Adding RDP Classifier data
# From the rrnDB website: "Estimate is an on-line interface to the RDP Classifier tool, including adjustment of 
# relative abundance of taxons based on 16S gene copy number data from rrnDB."
# "Estimate runs RDP Classifier version 2.13 using 16S training set #18 (or most recent) incorporating current rrnDB copy number data. 
# The necessary training files were re-created following the documented use of RDP Classifier's 'train' command, 
# replacing the copy number file from the original training set with one derived from the most recent downloadable 
# pan-taxa statistics."

# 1. https://rrndb.umms.med.umich.edu/estimate/run_classifier
# 2. Navigate to the “Download” page
# 3. Download “rrnDB-5.8_pantaxa_stats_RDP.tsv”, or the most recent version of this
# 4. Navigate back to rrnDB "Estimate" page and upload the fasta file for the full dataset (confidence cutoff 0.8)
# 5. Wait several minutes for it to run
# 6. Download the Classification assignment file (“dna-sequences.tsv”)


# Load in sequence classifications file
RDP = read.csv("rrnDB/dna-sequences_Exp3_12Jan2023.tsv",header=FALSE,sep=";")
head(RDP)
# We want to extract the genus they assigned for each OTU.
# Create function to extract genus, if present
GenusGenerator <- function(taxon) {
  Genus = strsplit(gsub(".*family\\s*|genus*", "", taxon),split="\t")[[1]][2]
  return(Genus)
}

# Extract the genus
RDP$GenusRDP = sapply(RDP$V1,GenusGenerator)
head(RDP)

# Might as well pull out OTU ID to be certain
OTUGenerator <- function(taxon) {
  OTU = strsplit(paste(taxon),split="\t")[[1]][1]
  return(OTU)
}
RDP$OTU = sapply(RDP$V1,OTUGenerator)
head(RDP)

# Ok, we've got what we need.
# Trim it down
RDP = RDP[,c(2,3)]

# Can now pull data from RRNDB.
# Reading in the rrnDB v5.8 file
rrnDB = read.csv("rrnDB/rrnDB-5.8_pantaxa_stats_RDP.tsv",sep="\t")
head(rrnDB)

# Creating a list of genera in the DB
rrnDBGenera = as.character(rrnDB[rrnDB$rank=="genus",]$name)

# Matching up genus name with mean predicted copy number
for (i in 1:length(RDP$GenusRDP)){
  GenusRDP = paste(RDP$GenusRDP[i])
  CopyNum = ifelse(GenusRDP %in% rrnDBGenera, rrnDB[rrnDB$name==GenusRDP,9],"")
  RDP$CopyNum[i] = CopyNum
}
tail(RDP)


# Bring in ps object and normalize to relative abundances
ps.full <- readRDS("ps.Exp3")
ps.norm <- transform_sample_counts(ps.full, function(x) x / sum(x))

# Work with melted phyloseq object
mdf = psmelt(ps.norm)


# Add the rrnDB copy number data to the melted phyloseq object
mdf = plyr::join(mdf,RDP,by="OTU")
mdf$CopyNum = as.numeric(mdf$CopyNum)
mdf$Abundance = as.numeric(mdf$Abundance)

head(mdf)

# From Nemergut et al. (2016) - OTU data were then normalized (standardized) for copy number 
# by dividing by copy number. For each sample, we calculated the community aggregated trait value 
# (weighted mean) by taking the product of the estimated operon copy number and the relative abundance 
# for each OTU, and summing this value across all OTUs in a sample. 

# So, first, we divide abundance by copy number
# Then, we re-calculate the relative abundanace, now adjusted for copy number
# The risk there, is, for any organisms without assigned copy numbers, they are excluded from this calculation.
# However, I think we have a pretty good fraction of the community with copy numbers
# To check:

d = mdf %>%
  dplyr::group_by(site, trt, fraction, Sample)%>%
  dplyr::filter(is.na(CopyNum))%>%
  dplyr::summarize(NoCopyNum = sum(Abundance))
hist(d$NoCopyNum)
d
# Not too bad, but a wide range



# Calculating weighted mean copy numbers:
df = mdf %>%
  dplyr::filter(!is.na(CopyNum))%>%
  dplyr::mutate(WtAbund = Abundance/CopyNum)%>%
  dplyr::group_by(site, trt, fraction, Sample)%>%
  dplyr::mutate(AdjAbund = WtAbund/sum(WtAbund))%>%
  dplyr::mutate(WtCopyNum = CopyNum*AdjAbund)%>%
  dplyr::summarize(WtMeanCopyNum = sum(WtCopyNum,na.rm=TRUE))
hist(df$WtMeanCopyNum)

#write.csv(df, 'Derived_data/Weighted_Mean_16S_Gene_Copy_number_EXP3_Jan2023.csv')

#df <- read.csv(file = 'Derived_data/Weighted_Mean_16S_Gene_Copy_number_EXP3_Jan2023.csv', row.names = 1)
head(df)

# Subset for Amynthas pressure sites
df = subset(df, site == "Arb" | site == "Lak")
# Rename Sample column as Sample ID
df$SampleID = df$Sample
df=df[-4]


# Rename and put the variables in order
df$trt = recode_factor(df$trt, NoWorm = "Control",LowWorm = "Low",HighWorm = "High")
df$trt = ordered(df$trt, levels = c("Control", "Low", "High"))
df$fraction = recode_factor(df$fraction, fresh = "Bulk \nsoil",
                              freem = "Free \nmicro.",
                              occm = "Occluded \nmicro.")
df$fraction = ordered(df$fraction, levels=c("Bulk \nsoil", "Free \nmicro.",
                                            "Occluded \nmicro."))
df$site = recode_factor(df$site, Arb = "Gallistel", Lak = "Lakeshore")
df$site = ordered(df$site, levels=c("Gallistel", "Lakeshore"))
AmynthasTitle = expression(paste(italic("Amynthas"), " pressure"))
trtpalette = c("#fcb18c","#B7003F","#2d105a") #Worm sites

# Add plot and sample metadata
meta = read.csv(file = 'sample_metadata_Exp3.csv', row.names = NULL)
head(meta)
meta = meta[-c(3,4,7,8,9)]

df = merge(df, meta, by="SampleID")

#### Worm disturbance sites: Graph of gene copy numbers, within plot and fraction ####
p.4 = ggplot(subset(df, (site == "Gallistel" | site == "Lakeshore")),
             aes(x=fraction, y=WtMeanCopyNum, color=trt))
p.4 = p.4 + geom_boxplot(alpha = 0.7) #+ geom_jitter(alpha = 0.1) 
p.4 = p.4 + scale_y_continuous(limits = c(1.2, 1.7),
                               breaks = c(1.2,1.3,1.4,1.5, 1.6, 1.7)) +
  labs(x="Soil fraction", y="Weighted mean predicted\n16S rRNA gene copy number", title=NULL)
p.4 = p.4 + facet_grid(~site)
p.4 = p.4 + scale_color_manual(AmynthasTitle, values=trtpalette)
p.4 = p.4 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
#p.4 = p.4 + theme(legend.position=c(0.75,0.8))
p.4 #750 x 350

ggsave("Figures/GeneCopyNumb.Amynthas.tiff", width=7.5, height=3.5, units = "in", device='tiff', dpi=400)

### Try graphing with means, SE
# Find mean
head(df)
df.mean = df%>%
  group_by(site, trt, fraction) %>% 
  dplyr::summarize(n=n(),
                   "Copies" = mean(WtMeanCopyNum, na.rm=TRUE)) %>%
  gather("element", "mean", - c(site, trt, fraction, n), factor_key=TRUE)  
df.mean

# Find SE
df.SE = df%>%
  group_by(site, trt, fraction) %>% 
  dplyr::summarize(n=n(),
                   "Copies" = sd(WtMeanCopyNum, na.rm=TRUE)/sqrt(n())) %>%
  gather("element", "SE", - c(site, trt, fraction, n), factor_key=TRUE)  
df.SE

# Join mean proportion and SE data together
df.meanSE = merge(df.mean, df.SE, by=c("site","trt", "fraction", "n", "element"))      
head(df.meanSE)

p = ggplot(df.meanSE, aes(x=factor(fraction), y=mean, color=trt)) +
  geom_point(size=3, position = dodge) +
  geom_errorbar(aes(ymin=mean-1.96*SE,ymax=mean+1.96*SE), position = dodge, width=0.2)
p = p + facet_grid(~site)
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white")) 
p = p + ylab("Weighted mean predicted\n16S rRNA gene copy number") + ylim(1.3,1.55)
p = p + xlab("Soil fraction")
p = p + scale_color_manual(AmynthasTitle, values=trtpalette)
p = p + theme(legend.position=c(0.815,0.845), legend.text=element_text(size=rel(0.75)), legend.title = element_text(size=rel(0.85)))
p = p + guides(color = guide_legend(byrow = TRUE)) + theme(legend.spacing.y = unit(0.0, 'cm'))
p # 450 x 450

ggsave("Figures/GeneCopyNumb.meanSE.Amynthas.tiff", width=4.5, height=4.5, units = "in", device='tiff', dpi=400)




### Summarize the results

# Calculating means across treatments:
head(df)
df.mean = df %>%
  #dplyr::filter(!is.na(CopyNum))%>%
  dplyr::group_by(site, trt, fraction)%>%
  dplyr::summarize(mean_WtMeanCopyNum=mean(WtMeanCopyNum))
df.mean



# Statistics, repeat for both sites
dat = subset(df, site == "Lakeshore")
ano = aov(WtMeanCopyNum ~ trt*fraction, data = dat)
summary(ano)
model.tables(ano, "means")

par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

cld = multcompLetters4(ano, tuk)
cld
