# Load required packages
library(ggplot2)
library(phyloseq)
library(dplyr)
library(picante)
library(DESeq2) #for rarefaction


# Read in phyloseq object
ps = readRDS("ps.Exp3")
ps

### Prune to just one site at a time
ps.site = subset_samples(ps,sample_data(ps)$site == "Lan")

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


# read in the phylogeny
# Rooted and unrooted trees seem to produce nearly identical results. Must modify (include.root) option.
phylo = ggtree::read.tree("tree_Exp3_fewer_quotes.nwk") #unrooted tree
#phylo = ggtree::read.tree("tree_rooted_Exp3.nwk") #rooted tree

### Notes on rooted vs unrooted trees from pd R documentation:
# If the root is to be included in all calculations of PD (include.root=TRUE), 
# the TREE MUST BE ROOTED. Single-species samples will be assigned a PD value equal
# to the distance from the root to the present.
#
# If the root is not included in all calculations by default (include.root=FALSE), 
# the tree need not rooted, but in the case of single-species samples the PD will 
# be equal to NA and a warning will be issued.

# ## Fix tip labels, if they contain extra set of quotes...this seemed to be an issue with the unrooted tree!
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
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

### Arlington
#               Df Sum Sq Mean Sq F value Pr(>F)  
# trt           1     40   40.50   0.484 0.4885  
# fraction      2    514  256.93   3.071 0.0516 .
# trt:fraction  2     26   12.82   0.153 0.8582  
# Residuals    84   7028   83.67                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### Lancaster
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# trt           1      1     0.9   0.011 0.915729    
# fraction      2   1632   816.0  10.657 0.000136 ***
# trt:fraction  2     53    26.6   0.347 0.708304    
# Residuals    51   3905    76.6       

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

### Arlington
# $fraction
#               diff        lwr         upr     p adj
# fresh-freem  4.116514  -1.518484  9.75151209 0.1954120
# occm-freem  -1.545097  -7.180096  4.08990079 0.7904944
# occm-fresh  -5.661611 -11.296609 -0.02661315 0.0486626

### Lancaster

#                 diff        lwr       upr     p adj
# fresh-freem  10.379084   3.525814 17.232354 0.0017279
# occm-freem   -1.741948  -8.595218  5.111322 0.8133807
# occm-fresh  -12.121032 -18.974302 -5.267762 0.0002483

### Option to save 
# write.csv(Faiths,"Derived_data/FaithsPD_tillage_Arl.csv")
# write.csv(Faiths,"Derived_data/FaithsPD_tillage_Lan.csv")


#### Load data, if not continuing from above ####
# Faiths = read.csv(file="Derived_data/FaithsPD_tillage_Lan.csv", row.names = 1)
# head(Faiths)

#### Plot Estimates ####

## Very basic plot
p = ggplot(Faiths, aes(y = PD, x = fraction, color = trt))
p = p + geom_point()
p


## Add SD bars
PDplot = Faiths %>%
    group_by(trt, fraction) %>%
    dplyr::summarize(PD_mean=mean(PD),
                     PD_sd=sd(PD))

p = ggplot(PDplot, aes(y = PD_mean, x = fraction, color = trt))
p = p + geom_point() + geom_errorbar(aes(ymin = PD_mean - PD_sd, ymax = PD_mean + PD_sd))
#p = p + ylim(2000,4500)
p

### Create a combined 'Faiths', with multiple sites, for facet graphing, below.
#Faiths.1 = read.csv(file="Derived_data/FaithsPD_tillage_Arl.csv", row.names = 1)
#Faiths.2 = read.csv(file="Derived_data/FaithsPD_tillage_Lan.csv", row.names = 1)

Faiths = rbind(Faiths.1, Faiths.2)

# Find mean and error, for graphing
PDplot = Faiths %>%
  group_by(site, trt, fraction) %>%
  dplyr::summarize(PD_mean=mean(PD),
                   PD_se=sd(PD)/sqrt(n()))
head(PDplot)

# Re-name variables and put them in order
PDplot$trt = recode_factor(PDplot$trt, "NoTill" = "No-tillage", 
                                 "ConvTill" = "Tillage")
PDplot$trt = ordered(PDplot$trt, levels = c("No-tillage", "Tillage"))
PDplot$fraction = recode_factor(PDplot$fraction,
                                      "fresh" = "Bulk \nsoil",
                                      "freem" = "Free \nmicro.",
                                      "occm" = "Occluded \nmicro.")
PDplot$fraction = ordered(PDplot$fraction, levels=c("Bulk \nsoil", "Free \nmicro.", "Occluded \nmicro."))
PDplot$site = recode_factor(PDplot$site,Arl = "Arlington",
                                  Lan = "Lancaster")
PDplot$site = ordered(PDplot$site, levels=c("Arlington", "Lancaster"))
trtpalette = c("#CC9933", "#117733")
dodge = position_dodge(0.5)

# Plot final richness estimates with errors: ±1.96SE should represent 95% confidence intervals
p = ggplot(PDplot, aes(y = PD_mean, x = fraction, color = trt))
p = p + geom_point(size=2, position = dodge) + geom_errorbar(aes(ymin=PD_mean-1.96*PD_se, ymax=PD_mean+1.96*PD_se), position=dodge, width=0.15)
p = p + facet_grid(~site)
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylim(105,145) + ylab("Faith's phylogenetic diversity")
p = p + xlab("Soil fraction")
p = p + scale_color_manual(values=trtpalette)
p = p + theme(legend.title = element_blank()) 
p # 500x300

#ggsave("Figures/PD.till.tiff", width=6, height=3.5, units = "in", device='tiff', dpi=400)
