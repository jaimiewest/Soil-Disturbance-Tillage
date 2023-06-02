# Load required packages

install.packages("devtools")
devtools::install_github("adw96/breakaway")
library(breakaway)
# You may need to update Breakaway and then close and reopen R
library(phyloseq)
library(tidyverse)
library(plyr)
library(dplyr) # load plyr, then dplyr in that order
library(ggplot2)

# Read in phyloseq object
ps = readRDS("ps.Exp3")
ps

### Prune to the tillage sites
ps = subset_samples(ps, sample_data(ps)$site == "Arl" | sample_data(ps)$site == "Lan")

# Prune away zero's
ps = prune_taxa(taxa_sums(ps) > 0, ps)
sample_data(ps)

# Collect sample data
SamDat = data.frame(sample_data(ps))

# take a peek at Chao1 estimates
ps2 = subset_samples(ps, sample_data(ps)$site == "Lan")
plot_richness(ps2, x = "fraction", color = 'trt', measures = "Chao1", title = "Chao1 Alpha diversity") + geom_boxplot()



#### Generate Breakaway Richness Estimates ####
#tutorial: https://adw96.github.io/breakaway/articles/breakaway.html

# Set OTU table
otu_data = t(otu_table(ps))
# Set sample data
meta_data = sample_data(ps)
# Had to flip OTU table so rownames match sample data
head(colnames(otu_data) == rownames(meta_data))

# Run Breakaway's frequency table list function
frequencytablelist = build_frequency_count_tables(otu_data)

# Check out one of them (#63)
head(frequencytablelist[[63]])

# Try Breakaway on a couple of samples
breakaway(frequencytablelist[[1]])
breakaway(frequencytablelist[[60]])

# Because no plot pops up, we know that we're dealing with the WLRM
# That's because dada2 won't allow singletons to pass through

# Run the richness estimator (breakaway) on all our samples (lists of frequency tables)
RichEsts = lapply(frequencytablelist,breakaway)
#library(purrr)
# Pull out the estimates, errors, and the model
Estimate = as.matrix(map(RichEsts, "estimate"))
Error = as.matrix(map(RichEsts, "error"))
Model = as.matrix(map(RichEsts, "model"))
df = data.frame(Estimate,Error,Model)

# Add sample ID column, estimate, and error
df$SampleID = row.names(df)
df$Estimate=as.numeric(df$Estimate)
df$Error=as.numeric(df$Error)

# Merge the estimates with the sample data
RichPlot3 = merge(SamDat,df,by="SampleID")
head(RichPlot3)

# # Save Breakaway Richness Estimates
RichPlot3$Model = as.character(RichPlot3$Model) # Model is a list, for some reason..and prevents write.csv
write.csv(RichPlot3,"Derived_data/Richness_Estimates_Tillage.csv")

#### Load Breakaway Richness Estimates, if not continuing from above ####
RichPlot3 = read.csv(file="Derived_data/Richness_Estimates_Tillage.csv", row.names = 1)
head(RichPlot3)
#### Plot Estimates ####

## Plot them a few ways
p = ggplot(RichPlot3,aes(y = Estimate, x = trt, color = fraction))
p = p + geom_point()
p

head(RichPlot3)
## Add SD bars
RichPlot = RichPlot3 %>%
    group_by(trt, fraction) %>%
    dplyr::summarize(Estimate_mean = mean(Estimate),
                     Estimate_sd = sd(Estimate))

p = ggplot(RichPlot, aes(y = Estimate_mean, x = trt, color = fraction))
p = p + geom_point() + geom_errorbar(aes(ymin = Estimate_mean - Estimate_sd, ymax = Estimate_mean + Estimate_sd))
p = p + ylim(2000, 4500)
p


### We want to summarize the estimates across at each site for each treatment and fraction
# First, make a single variable that has all of those elements
RichPlot3$Comp = paste(RichPlot3$site, RichPlot3$trt, RichPlot3$fraction)
RichPlot3$Comp = as.factor(RichPlot3$Comp)
head(RichPlot3)

# Create empty data frame that will hold our output data
RichPlotSumm = data.frame(Comp=levels(RichPlot3$Comp))
RichPlotSumm$Estimate = 0
RichPlotSumm$Error = 0
RichPlotSumm$p = 0
head(RichPlotSumm)

# Run a for loop that goes through each subset of data and makes the estimates using the betta function
for (i in levels(RichPlot3$Comp)){
  d = RichPlot3[RichPlot3$Comp==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSumm[RichPlotSumm$Comp==i,]$Estimate = Betta$table[,1]
  RichPlotSumm[RichPlotSumm$Comp==i,]$Error = Betta$table[,2]
  RichPlotSumm[RichPlotSumm$Comp==i,]$p = Betta$table[,3]
}
head(RichPlotSumm)
# Create a function to extract the coding values from our comparison variable (Comp)
# so we can get the site, trt and fraction info back
substrRight <- function(x, n){
  substr(x, nchar(x)-4, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, 3)
}
substrMid <- function(x, n){
  substr(x, 5, 11)
}
# Pull out variables from the Comp variable so we can use them as plotting variables
# (This won't be perfect because our trt and fraction variable levels have different number of characters)
RichPlotSumm$site = substrLeft(as.character(RichPlotSumm$Comp),1)
RichPlotSumm$trt = substrMid(as.character(RichPlotSumm$Comp),1)
RichPlotSumm$fraction = substrRight(as.character(RichPlotSumm$Comp),1)
RichPlotSumm

write.csv(RichPlotSumm,"Derived_data/Richness_Estimates_forplot_Tillage.csv")

#### Load Breakaway Richness Estimates, if not continuing from above ####
RichPlotSumm = read.csv(file="Derived_data/Richness_Estimates_forplot_Tillage.csv", row.names = 1)

RichPlotSumm

# Re-name variables and put them in order
RichPlotSumm$trt = recode_factor(RichPlotSumm$trt,
                       "NoTill " = "No-tillage", 
                       "ConvTil" = "Tillage")
RichPlotSumm$trt = ordered(RichPlotSumm$trt, levels = c("No-tillage", "Tillage"))

RichPlotSumm$fraction = recode_factor(RichPlotSumm$fraction,
                                      "fresh" = "Bulk \nsoil",
                                      "freem" = "Free \nmicro.",
                                      " occm" = "Occluded \nmicro.")
RichPlotSumm$fraction = ordered(RichPlotSumm$fraction, levels=c("Bulk \nsoil", "Free \nmicro.", "Occluded \nmicro."))

RichPlotSumm$site = recode_factor(RichPlotSumm$site,
                                  Arl = "Arlington",
                                  Lan = "Lancaster")

RichPlotSumm$site = ordered(RichPlotSumm$site, levels=c("Arlington", "Lancaster"))

trtpalette = c("#CC9933", "#117733")
dodge = position_dodge(0.5)


# Plot final richness estimates with errors: ±1.96SE should represent 95% confidence intervals
p = ggplot(RichPlotSumm, aes(y = Estimate, x = fraction, color = trt))
p = p + geom_point(size=2, position = dodge) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), position=dodge, width=0.15)
p = p + facet_grid(~site)
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylim(0,4250) + ylab("Richness estimate")
p = p + xlab("Soil fraction")
p = p + scale_color_manual(values = trtpalette)
p = p + theme(legend.title = element_blank()) #+ theme(legend.position = c(0.15, 0.15), legend.box = "horizontal")
p # 450x300

#ggsave("Figures/richness.tillage.tiff", width=6, height=3, units = "in", device='tiff', dpi=400)


### Statistics
# Arlington
head(RichPlot3)
dat = subset(RichPlot3, site == "Arl")
ano = aov(Estimate ~ trt*fraction, data = dat)
summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

# Lancaster
head(RichPlot3)
dat = subset(RichPlot3, site == "Lan")
ano = aov(Estimate ~ trt*fraction, data = dat)
summary(ano)
model.tables(ano, "means")
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk



