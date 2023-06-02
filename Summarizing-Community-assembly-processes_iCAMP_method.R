library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(magrittr)
library(harrietr)

#### Getting started #### 

### read in data, generated using separate codes for each site, using a high-throughput computing cluster
### Select one site at a time, then repeat for the other site.

# # Lancaster:
# iCAMP.IDs = read.csv("/Volumes/External/CHTC/Exp3/iCAMP/Lancaster/iCAMP.process.CbMNTDiCBraya.csv", header=TRUE, check.names=FALSE, row.names=1)

# Arlington:
iCAMP.IDs = read.csv("/Volumes/External/CHTC/Exp3/iCAMP/Arlington/iCAMP.process.CbMNTDiCBraya.csv", header=TRUE, check.names=FALSE, row.names=1)

head(iCAMP.IDs)
names(iCAMP.IDs) = c("SampleID", "com2", "Variable selection", "Homogeneous selection", "Dispersal limitation", 
                     "Homogenizing dispersal", "Undominated")

# Add sample metadata
meta = read.csv("~/Desktop/Community assembly/Exp3/sample_metadata_Exp3.csv", header=TRUE)
meta = meta [,-c(8,9)]

allcomps = merge(iCAMP.IDs, meta, by = "SampleID")
head(allcomps)
names(allcomps) = c("com1", "SampleID", "Variable selection", "Homogeneous selection", "Dispersal limitation", 
                "Homogenizing dispersal", "Undominated","com1.sample","com1.fraction","com1.site",
                "com1.plot","com1.core", "com1.trt")
allcomps = merge(allcomps, meta, by = "SampleID")
names(allcomps) = c("com2", "com1", "Variable selection", "Homogeneous selection", "Dispersal limitation", 
                "Homogenizing dispersal", "Undominated","com1.sample","com1.fraction","com1.site",
                "com1.plot","com1.core", "com1.trt",
                "com2.sample","com2.fraction","com2.site",
                "com2.plot","com2.core", "com2.trt")

# Make sorting columns
allcomps$sametrt = ifelse(allcomps$com1.trt==allcomps$com2.trt, "same", "no")
allcomps$sameplot = ifelse(allcomps$com1.plot==allcomps$com2.plot, "same", "no")
allcomps$samefraction = ifelse(allcomps$com1.fraction==allcomps$com2.fraction, "same", "no")
head(allcomps)

str(allcomps)
allcomps$com1.plot = as.character(allcomps$com1.plot)
allcomps$com2.plot = as.character(allcomps$com2.plot)
allcomps$com1.core = as.character(allcomps$com1.core)
allcomps$com2.core = as.character(allcomps$com2.core)

# Pull out comparisons within plot
same_plot = allcomps  %>%
  dplyr::filter(samefraction == "same") %>%
  dplyr::filter(sameplot == "same")

head(same_plot)
## ^^ this same_plot df can be used for statistical analysis! at the end of this code

# Find the mean proportion of each process ID, within plot
same_plot.gather = gather(same_plot, key="Process", value="Proportion", 3:7)
same_plot2 = same_plot.gather %>% 
  group_by(com1.trt, com1.fraction, Process, sameplot) %>%
  dplyr::summarise(n=n(),
            "Proportion" = mean(Proportion))
same_plot2

# Find the standard error of each process ID, within plot comparisons
same_plotSE = same_plot.gather %>%
  group_by(com1.trt, com1.fraction, Process, sameplot) %>%
  dplyr::summarise(n=n(),
                   "SE" = sd(Proportion)/sqrt(n()))
same_plotSE

# Join proportion and SE data together
same_plot2 = merge(same_plot2,same_plotSE, by=c("com1.trt", "com1.fraction", "Process",
                                                "sameplot","n"))
same_plot2


# Pull out comparisons between plots (within trt and fraction; excluding within-plot)
different_plot = allcomps %>%
  dplyr::filter(sametrt == "same") %>%
  dplyr::filter(samefraction == "same") %>%
  dplyr::filter(sameplot == "no")

head(different_plot)
## ^^ this different_plot df can be used for statistical analysis! at the end of this code


different_plot.gather = gather(different_plot, key="Process", value="Proportion", 3:7)

# Find the mean proportion of relative abundance for each process ID, by treatment
different_plot2 = different_plot.gather %>% 
  group_by(com1.trt, com1.fraction, Process, sameplot) %>%
  dplyr::summarise(n=n(),
            "Proportion" = mean(Proportion))
different_plot2

# Find the standard error of each process ID, across plot comparisons
different_plotSE = different_plot.gather %>%
  group_by(com1.trt, com1.fraction, Process, sameplot) %>%
  dplyr::summarise(n=n(),
                   "SE" = sd(Proportion)/sqrt(n()))
different_plotSE

# Join proportion and SE data together
different_plot2 = merge(different_plot2, different_plotSE, by=c("com1.trt", "com1.fraction", "Process",
                                                "sameplot","n"))
different_plot2


### Merge summary datasets
Proc.ID = rbind(same_plot2, different_plot2)

## Option to save
#write.csv(Proc.ID,"Derived_data/Community_Assembly_Process_IDs_iCAMP_Arlington.csv")

# Read in data, if not continuing from above
Proc.ID = read.csv(file ="Derived_data/Community_Assembly_Process_IDs_iCAMP_Lancaster.csv",
                           header=TRUE, stringsAsFactors=TRUE, row.names=1)

#### STACKED BAR GRAPH--Within Plot ####
Proc.ID
# Re-name variables and put them in order
Proc.ID$com1.trt = recode_factor(Proc.ID$com1.trt, "NoTill" = "No-tillage", "ConvTill" = "Tillage")
Proc.ID$com1.trt = ordered(Proc.ID$com1.trt, levels = c("No-tillage", "Tillage"))
Proc.ID$com1.fraction = ordered(Proc.ID$com1.fraction, levels = c("fresh", "freem", "occm"))
fraction.labs = c("fresh" = "Bulk soil", "freem" = "Free microaggregate", "occm" = "Occluded microagg.")
Proc.ID$Process = ordered(Proc.ID$Process, levels=c("Undominated","Homogeneous selection",
                                                    "Homogenizing dispersal",
                                                    "Dispersal limitation",
                                                    "Variable selection"))
dodge = position_dodge2(width = 0.3, padding = 0.3)
ID_Palette = c("#CCCCCC", "#99CC99","#336699", "#FFCC66", "#CC3333")

head(Proc.ID)
### Stacked bar, comparisons made within the same plot
p.2 = ggplot(subset(Proc.ID, sameplot=="same"), aes(x=com1.trt, y=Proportion)) +
  geom_bar(position = "stack", stat="identity",aes(fill=Process)) +
  scale_fill_manual(values = ID_Palette) +
  labs(x = "Within-plot sample comparisons",
       y = "Relative influence, Arlington", fill = "Community assembly process") +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  theme_bw() +
  facet_grid(~com1.fraction, labeller=labeller(com1.fraction=fraction.labs)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
p.2 #530 x 330 for ms

###Name plot, depending on site
#Arl.iCAMP.plot = p.2
#Lan.iCAMP.plot = p.2



### Stacked bar, comparisons made between plots/trt level
q.2 = ggplot(subset(Proc.ID, sameplot=="no"), aes(x=com1.trt, y=Proportion)) +
  geom_bar(position = "stack", stat="identity",aes(fill=Process)) +
  scale_fill_manual(values = ID_Palette) +
  labs(x = "Between-plot sample comparisons",
       y = NULL, fill = "Community assembly process") +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  theme_bw() +
  facet_grid(~com1.fraction, labeller=labeller(com1.fraction=fraction.labs)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white"))
q.2 #530 x 330 for ms

###Name plot, depending on site
#Arl.iCAMP.trt = q.2
#Lan.iCAMP.trt = q.2


######################### #
#### Arrange Community Assembly figs ####
######################### #

ggarrange(Arl.iCAMP.plot, Arl.iCAMP.trt,
  Lan.iCAMP.plot, Lan.iCAMP.trt,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom")

ggsave("Figures/Community_Assembly.iCAMP_Tillage.tiff", width=9.5, height=7.5, units = "in", device='tiff', dpi=400)



### Statistics-- identify significant changes to homogenizing dispersal and homogeneous selection
# must repeat for both Arlington and Lancaster data
head(same_plot)

#rename the process ID columns to drop the space
names(same_plot) = c("com2", "com1", "Variable.selection", "Homogeneous.selection", "Dispersal.limitation", 
                    "Homogenizing.dispersal", "Undominated","com1.sample","com1.fraction","com1.site",
                    "com1.plot","com1.core", "com1.trt",
                    "com2.sample","com2.fraction","com2.site",
                    "com2.plot","com2.core", "com2.trt", "sametrt", "sameplot", "samefraction")

# Subset by fraction. Must manually repeat for fresh, freem, occm.
# (and then be sure we are only considering same fraction comparisons)
# must also repeat for same_plot and different_plot datasets

dat = subset(different_plot, com1.fraction=="fresh")

#create model--must repeat manually for each process of interest (e.g., Homogenizing.dispersal, Homogeneous.selection)
ano = aov(Dispersal.limitation ~ com1.trt, data = dat)
summary(ano)
model.tables(ano, "means")
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

tuk = TukeyHSD(ano, conf.level = 0.95)
#tuk

#library(multcompView)
cld = multcompLetters4(ano, tuk)
cld

