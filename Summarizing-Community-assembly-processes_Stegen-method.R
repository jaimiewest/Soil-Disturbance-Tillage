
#library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(magrittr)
library(harrietr)
library(ggpubr) # for ggarrange

#### Getting started #### 


### read in data beta nearest taxon distance (bNTI) and Raup-Crick Bray-Curtis (RCBC) data, 
### generated in separate R codes usign a high throughput computing cluster

### Do one site at a time (i.e. only load in data for one site at a time)
### Arlington
bNTI = read.csv("Derived_data/Community_assembly_Arl/weighted_bNTI_Lan_Hel.csv", header=TRUE, check.names=FALSE, row.names = 1)
RCBC = read.csv("Derived_data/Community_assembly_Arl/RC_results_Lan_Hel.csv", header=TRUE, check.names=FALSE, row.names = 1)

# # Lancaster
# bNTI = read.csv("Derived_data/Community_assembly_Lan/weighted_bNTI_Lan_Hel.csv", header=TRUE, check.names=FALSE, row.names = 1)
# RCBC = read.csv("Derived_data/Community_assembly_Lan/RC_results_Lan_Hel.csv", header=TRUE, check.names=FALSE, row.names = 1)

# Shape up data
#head(bNTI)
bNTI2 = melt_dist(bNTI, dist_name = "bNTI")
head(bNTI2)

names(RCBC) = sub("^X", "", names(RCBC))
RCBC2 = melt_dist(RCBC, dist_name = "RCBC")
head(RCBC2)

# Now, merge bNTI and RCBC dataframes
merged = merge(bNTI2,RCBC2,by=c("iso1", "iso2"), all.y =TRUE)
head(merged)
names(merged) = c("SampleID", "com2", "bNTI", "RCBC")

# Add sample metadata
meta = read.csv("sample_metadata_Exp3.csv", header=TRUE)
meta = meta [,-c(8,9)]
#head(meta)

test = merge(merged, meta, by = "SampleID")
head(test)
names(test) = c("com1", "SampleID", "bNTI", "RCBC","com1.sample","com1.fraction","com1.site",
                "com1.plot","com1.core", "com1.trt")
test = merge(test, meta, by = "SampleID")
names(test) = c("com2", "com1", "bNTI", "RCBC","com1.sample","com1.fraction","com1.site",
                "com1.plot","com1.core", "com1.trt",
                "com2.sample","com2.fraction","com2.site",
                "com2.plot","com2.core", "com2.trt")
head(test)

# Make sorting columns
test$samefraction = ifelse(test$com1.fraction==test$com2.fraction, "same", "no")
test$sameplot = ifelse(test$com1.plot==test$com2.plot, "same", "no")
test$sametrt = ifelse(test$com1.trt==test$com2.trt, "same", "no")
str(test)

test$com1.plot = as.character(test$com1.plot)
test$com2.plot = as.character(test$com2.plot)
test$com1.core = as.character(test$com1.core)
test$com2.core = as.character(test$com2.core)

# Pull out comparisons within plot (and fraction)
same_plot = test %>%
  dplyr::filter(sameplot == "same") %>%
  dplyr::filter(samefraction == "same")

# Pull out comparisons across plots (and within same trt & fraction; excluding within-plot)
different_plot = test %>%
  dplyr::filter(sameplot == "no") %>%
  dplyr::filter(sametrt == "same") %>%
  dplyr::filter(samefraction == "same")

#head(same_plot)

### Assign process IDs; considering selection first, and dispersal only if no evidence for selection

## Within same plot

# create empty column for the process ID
same_plot$Process <- NA

for(r in 1:(nrow(same_plot))) {
  if (same_plot$bNTI[r] < -2) {
    same_plot$Process[r] = "Homogeneous selection"
  } else if (same_plot$bNTI[r] > 2) {
    same_plot$Process[r] = "Variable selection"
  } else if (same_plot$RCBC[r] < -0.95 && same_plot$bNTI[r] >= -2 && same_plot$bNTI[r] <= 2) {
    same_plot$Process[r] = "Homogenizing dispersal"
  } else if (same_plot$RCBC[r] > 0.95 && same_plot$bNTI[r] >= -2 && same_plot$bNTI[r] <= 2) {
    same_plot$Process[r] = "Dispersal limitation"
  } else if (same_plot$RCBC[r] >= -0.95 && same_plot$RCBC[r] <= 0.95 && same_plot$bNTI[r] >= -2 && same_plot$bNTI[r] <= 2) {
    same_plot$Process[r] = "Undominated"
  } else same_plot$Process[r] = NA
}

head(same_plot)

# Find the mean proportion of each process ID, by treatment
ProcID.same.plot = same_plot %>% 
  group_by(com1.trt, com1.fraction, Process, sameplot) %>%
  dplyr::summarise(n=n())
ProcID.same.plot


### ...Within same treatment and fraction, but not same plot

# create empty column for the process ID
different_plot$Process <- NA

for(r in 1:(nrow(different_plot))) {
  if (different_plot$bNTI[r] < -2) {
    different_plot$Process[r] = "Homogeneous selection"
  } else if (different_plot$bNTI[r] > 2) {
    different_plot$Process[r] = "Variable selection"
  } else if (different_plot$RCBC[r] < -0.95 && different_plot$bNTI[r] >= -2 && different_plot$bNTI[r] <= 2) {
    different_plot$Process[r] = "Homogenizing dispersal"
  } else if (different_plot$RCBC[r] > 0.95 && different_plot$bNTI[r] >= -2 && different_plot$bNTI[r] <= 2) {
    different_plot$Process[r] = "Dispersal limitation"
  } else if (different_plot$RCBC[r] >= -0.95 && different_plot$RCBC[r] <= 0.95 && different_plot$bNTI[r] >= -2 && different_plot$bNTI[r] <= 2) {
    different_plot$Process[r] = "Undominated"
  } else different_plot$Process[r] = NA
}

head(different_plot)

# Find the mean proportion of each process ID, within treatments
ProcID.diff.plot = different_plot %>% 
  group_by(com1.trt, com1.fraction, Process, sameplot) %>%
  dplyr::summarise(n=n())
ProcID.diff.plot
ProcID.same.plot

# Make a dataframe specific to this site, before repeating the above code for a different site
Proc.ID = rbind(ProcID.diff.plot, ProcID.same.plot)
Proc.ID

## Option to save
write.csv(Proc.ID,"Derived_data/Community_Assembly_Process_IDs_Arl_summary.csv")
#write.csv(Proc.ID,"Derived_data/Community_Assembly_Process_IDs_Lan_summary.csv")



### I did not graph these data because it was nearly all homogeneous selection, 
### and instead made figures for the data generated using the Ning et al., 2020 iCAMP method