library(Biostrings)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(multcompView) # to generate compact letter display (using aov and tukey's)

##### We want to see how the weighted mean predicted 16S copy number corresponds to tillage treatment and soil fraction

#### Adding RDP Classifier data ####

# From the rrnDB website: "Estimate is an on-line interface to the RDP Classifier tool, including adjustment of 
# relative abundance of taxons based on 16S gene copy number data from rrnDB."
# "Estimate runs RDP Classifier version 2.12 using 16S training set #16 incorporating current rrnDB copy number data. 
# The necessary training files were re-created following the documented use of RDP Classifier's 'train' command, 
# replacing the copy number file from the original training set with one derived from the most recent downloadable 
# pan-taxa statistics."

# 1. https://rrndb.umms.med.umich.edu/estimate/run_classifier
# 2. Navigate to the “Download” page
# 3. Download “rrnDB-5.7_pantaxa_stats_RDP.tsv”, or the most recent version of this
# 4. Navigate back to rrnDB "Estimate" page and upload the fasta file for the full dataset (confidence cutoff 0.8)
# 5. Wait several minutes for it to run
# 6. Download the Classification assignment file (“dna-sequences.tsv”)


# Load in sequence classifications file
RDP = read.csv("rrnDB/dna-sequences_Exp3.tsv",header=FALSE,sep=";")
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
# Reading in the rrnDB file
rrnDB = read.csv("rrnDB/rrnDB-5.7_pantaxa_stats_RDP.tsv",sep="\t")
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

# Working with melted phyloseq object
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
# To check:

d = mdf %>%
  dplyr::group_by(site, trt, fraction, Sample)%>%
  dplyr::filter(is.na(CopyNum))%>%
  dplyr::summarize(NoCopyNum = sum(Abundance))
hist(d$NoCopyNum)
d
# Not too bad

# Check out distribution of NoCopyNum across treatments.
# create a df that excludes the trt columns (for the grey background full dataset)
# it would be good to check out individual sites, but this gives us an idea.
head(dat.notrt)
dat.notrt=d[, -2]
dat.notrt=dat.notrt[, -2]
p = ggplot(d, aes(x=NoCopyNum, fill=trt)) +
  geom_histogram(data = dat.notrt, fill = "grey", alpha = 0.5) +
  geom_histogram(color = "black") +
  facet_wrap(~trt) +
  theme_bw()
p = p + scale_color_manual(values=palette)
p = p + xlab("Proportion of relative abundance for which OTUs lack a predicted copy number")
p = p + ylab("Count, number of samples")
p

# Tillage communities have a higher proportion of copy numbers assigned than No-till communities

# Calculating weighted mean copy numbers:
df = mdf %>%
  dplyr::filter(!is.na(CopyNum))%>%
  dplyr::mutate(WtAbund = Abundance/CopyNum)%>%
  dplyr::group_by(site, trt, fraction, Sample)%>%
  dplyr::mutate(AdjAbund = WtAbund/sum(WtAbund))%>%
  dplyr::mutate(WtCopyNum = CopyNum*AdjAbund)%>%
  dplyr::summarize(WtMeanCopyNum = sum(WtCopyNum,na.rm=TRUE))
hist(df$WtMeanCopyNum)

# Option to save
#write.csv(df, 'Derived_data/Weighted_Mean_16S_Gene_Copy_number.csv')
#df = read.csv(file = 'Derived_data/Weighted_Mean_16S_Gene_Copy_number.csv', row.names = 1)

head(df)

# Rename Sample column as Sample ID
df$SampleID = df$Sample
df=df[-4]

# Put the variables in order
df$trt = recode_factor(df$trt, NoTill = "No-tillage", ConvTill = "Tillage")
df$trt = ordered(df$trt, levels = c("No-tillage", "Tillage"))
trtpalette = c("#CC9933", "#117733")
df$fraction = ordered(df$fraction, levels=c("fresh", "freem", "occm"))
df$site = ordered(df$site, levels=c("Arl", "Lan"))
site.labs=c("Arl" = "Arlington", "Lan" = "Lancaster")
fraction.labs=c("fresh" = "Bulk \nsoil", "freem" = "Free \nmicro", "occm" = "Occluded \nmicro")

# Add plot and sample metadata
meta = read.csv(file = 'sample_metadata_Exp3.csv', row.names = NULL)
head(meta)
meta = meta[-c(3,4,7,8,9)]

df = merge(df, meta, by="SampleID")

#### Tillage sites: Graph of gene copy numbers, within plot and fraction ####
p.4 = ggplot(df ,aes(x = fraction, y = WtMeanCopyNum, color = trt))
p.4 = p.4 + geom_boxplot(alpha = 0.7) #+ geom_jitter(alpha = 0.1) 
p.4 = p.4 + expand_limits(y=c(1.5,2.02)) +
  labs(x="Soil fraction", y="Weighted mean predicted\n16S rRNA gene copy number", title=NULL) +
  scale_x_discrete(labels=fraction.labs)
p.4 = p.4 + facet_grid(~site, labeller = labeller(site=site.labs))
p.4 = p.4 + scale_color_manual(values=trtpalette)
p.4 = p.4 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.4 = p.4 + theme(legend.position=c(.2,.87)) + theme(legend.title=element_blank())
p.4 

ggsave("Figures/GeneCopyNumb.till.tiff", width=4, height=3.65, units = "in", device='tiff', dpi=400)


### Summarize the results

# Calculating means across treatments:
head(df)
df.mean = df %>%
  #dplyr::filter(!is.na(CopyNum))%>%
  dplyr::group_by(site, trt, fraction)%>%
  dplyr::summarize(mean_WtMeanCopyNum = mean(WtMeanCopyNum))
df.mean



# Statistics, ANOVA, Tukey's, generate CLD
# Must manually update and run for each site
dat = subset(df, site == "Lan")
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
