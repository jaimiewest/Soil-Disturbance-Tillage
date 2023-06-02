library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)
library(ggpubr)
library(multcompView)

# Read in phyloseq object
ps = readRDS("ps.Exp3")

######################## #
#### Subset for Arlington ####
######################## #
ps.A = subset_samples(ps,sample_data(ps)$site == "Arl")
# Prune away zero's
ps.A = prune_taxa(taxa_sums(ps.A) > 0, ps.A)
# Normalize using Hellinger transformation
ps.A = transform_sample_counts(ps.A, function(x) (x / sum(x))^0.5 )


######################## #
#### Subset for Lancaster ####
######################## #
ps.L = subset_samples(ps,sample_data(ps)$site == "Lan")
# Prune away zero's
ps.L = prune_taxa(taxa_sums(ps.L) > 0, ps.L)
# Normalize using Hellinger transformation
ps.L = transform_sample_counts(ps.L, function(x) (x / sum(x))^0.5 )


######################## #
#### Distance to Spatial Median (centroid), Arlington--between-plot scale ("field scale") ####
######################## #
# Must do this for each fraction separately--freem, occm, and fresh, as follows.

## First, create veganotu function
# veganotu = function(physeq) {
#   require("vegan")
#   OTU = otu_table(physeq)
#   if (taxa_are_rows(OTU)) {
#     OTU = t(OTU)
#   }
#   return(as(OTU, "matrix"))
# }

### Free microaggregate, Arlington
ps.A.fr = subset_samples(ps.A, fraction == "freem")
DistVar.A.fr = vegdist(veganotu(ps.A.fr), method = "bray")
ps.A.fr.df = data.frame(sample_data(ps.A.fr))

betadisp.A.fr = betadisper(DistVar.A.fr, ps.A.fr.df$trt)

trt = data.frame(betadisp.A.fr$group)
distances = data.frame(betadisp.A.fr$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Arl"
data$fraction= "freem"
data.Arl.fr = data

### Occluded microaggregate, Arlington
ps.A.oc = subset_samples(ps.A, fraction == "occm")
DistVar.A.oc = vegdist(veganotu(ps.A.oc), method = "bray")
ps.A.oc.df = data.frame(sample_data(ps.A.oc))

betadisp.A.oc = betadisper(DistVar.A.oc, ps.A.oc.df$trt)

trt = data.frame(betadisp.A.oc$group)
distances = data.frame(betadisp.A.oc$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Arl"
data$fraction= "occm"
data.Arl.oc = data

### Bulk soil ("fresh"), Arlington
ps.A.bu = subset_samples(ps.A, fraction == "fresh")
DistVar.A.bu = vegdist(veganotu(ps.A.bu), method = "bray")
ps.A.bu.df = data.frame(sample_data(ps.A.bu))

betadisp.A.bu = betadisper(DistVar.A.bu, ps.A.bu.df$trt)

trt = data.frame(betadisp.A.bu$group)
distances = data.frame(betadisp.A.bu$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Arl"
data$fraction= "fresh"
data.Arl.bu = data


######################## #
#### Distance to Centroid, Lancaster, between-plot scale ("field scale") ####
######################## #

### Free microaggregate, Lancaster
ps.L.fr = subset_samples(ps.L, fraction == "freem")
DistVar.L.fr = vegdist(veganotu(ps.L.fr), method = "bray")
ps.L.fr.df = data.frame(sample_data(ps.L.fr))

betadisp.L.fr = betadisper(DistVar.L.fr, ps.L.fr.df$trt)

trt = data.frame(betadisp.L.fr$group)
distances = data.frame(betadisp.L.fr$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Lan"
data$fraction= "freem"
data.Lan.fr = data

### Occluded microaggregate, Lancaster
ps.L.oc = subset_samples(ps.L, fraction == "occm")
DistVar.L.oc = vegdist(veganotu(ps.L.oc), method = "bray")
ps.L.oc.df = data.frame(sample_data(ps.L.oc))

betadisp.L.oc = betadisper(DistVar.L.oc, ps.L.oc.df$trt)

trt = data.frame(betadisp.L.oc$group)
distances = data.frame(betadisp.L.oc$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Lan"
data$fraction= "occm"
data.Lan.oc = data

### Bulk soil, Lancaster
ps.L.bu = subset_samples(ps.L, fraction == "fresh")
DistVar.L.bu = vegdist(veganotu(ps.L.bu), method = "bray")
ps.L.bu.df = data.frame(sample_data(ps.L.bu))

betadisp.L.bu = betadisper(DistVar.L.bu, ps.L.bu.df$trt)

trt = data.frame(betadisp.L.bu$group)
distances = data.frame(betadisp.L.bu$distances)
data = cbind(trt,distances)
colnames(data) = c("trt", "distances")

data$site= "Lan"
data$fraction= "fresh"
data.Lan.bu = data

# Combine the dataframes
data.bytrt = rbind(data.Arl.fr, data.Lan.fr, data.Arl.oc, data.Lan.oc, data.Arl.bu, data.Lan.bu)
head(data.bytrt)

### Option to save 
#write.csv(data.bytrt,"Derived_data/Exp3_Till_BC_Distance_to_centroid_within_trt.csv")
head(data.bytrt)
#### Load data, if not continuing from above ####
data.bytrt = read.csv(file="Derived_data/Exp3_Till_BC_Distance_to_centroid_within_trt.csv", row.names=1)
head(data.bytrt)

################################## #
#### Boxplots of Distance to Centroid, between-plot scale
################################## #

# Re-name variables and put them in order
data.bytrt$trt = recode_factor(data.bytrt$trt, "NoTill" = "No-tillage", "ConvTill" = "Tillage")
data.bytrt$trt = ordered(data.bytrt$trt, levels = c("No-tillage", "Tillage"))
data.bytrt$fraction = recode_factor(data.bytrt$fraction, "fresh" = "Bulk \nsoil", "freem" = "Free \nmicroagg.", "occm" = "Occluded \nmicroagg.")
data.bytrt$fraction = ordered(data.bytrt$fraction, levels=c("Bulk \nsoil", "Free \nmicroagg.", "Occluded \nmicroagg."))
data.bytrt$site = ordered(data.bytrt$site, levels=c("Arl", "Lan"))
trtpalette = c("#CC9933", "#117733")

disp.trt = ggplot(subset(data.bytrt, site =="Arl"), aes(x = fraction, y = distances, color = trt ))
disp.trt = disp.trt + geom_boxplot()
disp.trt = disp.trt + facet_wrap(~ site, labeller = as_labeller(c("Arl" = "Arlington, between-plot scale \n(i.e., treatment scale)",
                                                                  "Lan" = "Lancaster, between-plot scale \n(i.e., treatment scale)"))) +
  labs(x = "", y = "Distance to spatial median", fill = "") +
  scale_color_manual(values = trtpalette,
                     name = "Treatment") +
  scale_fill_manual(name=NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background = element_rect(fill="white"))+
  ylim(0.14, 0.345)
disp.trt # 550 x 370

## Re-run ggplot above separately for each site
## Save the Arlington between-plot figure:
#disp.trt.Arl = disp.trt
## Save the Lancaster between-plot figure:
#disp.trt.Lan = disp.trt


######################## #
#### Calculate Distance to Centroid, WITHIN PLOT-scale, Arlington ####
######################## #
# Again, must do this for each soil fraction--freem, occm, and fresh, as follows
# (veganotu function was created above)

### Free microaggregates, Arlington
ps.A.fr = subset_samples(ps.A, fraction == "freem")
DistVar.A.fr = vegdist(veganotu(ps.A.fr), method = "bray")
ps.A.fr.df = data.frame(sample_data(ps.A.fr))

betadisp.A.fr = betadisper(DistVar.A.fr, ps.A.fr.df$plot)
betadisp.A.fr
plot = data.frame(betadisp.A.fr$group)
distances = data.frame(betadisp.A.fr$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Arl"
data$fraction= "freem"
data.Arl.fr = data

### Occluded microaggregates, Arlington
ps.A.oc = subset_samples(ps.A, fraction == "occm")
DistVar.A.oc = vegdist(veganotu(ps.A.oc), method = "bray")
ps.A.oc.df = data.frame(sample_data(ps.A.oc))

betadisp.A.oc = betadisper(DistVar.A.oc, ps.A.oc.df$plot)

plot = data.frame(betadisp.A.oc$group)
distances = data.frame(betadisp.A.oc$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Arl"
data$fraction= "occm"
data.Arl.oc = data

### Bulk soil, Arlington
ps.A.bu = subset_samples(ps.A, fraction == "fresh")
DistVar.A.bu = vegdist(veganotu(ps.A.bu), method = "bray")
ps.A.bu.df = data.frame(sample_data(ps.A.bu))

betadisp.A.bu = betadisper(DistVar.A.bu, ps.A.bu.df$plot)

plot = data.frame(betadisp.A.bu$group)
distances = data.frame(betadisp.A.bu$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Arl"
data$fraction= "fresh"
data.Arl.bu = data


######################## #
#### Distance to Centroid, WITHIN PLOT, Lancaster ####
######################## #


### Free microaggregates, Lancaster
ps.L.fr = subset_samples(ps.L, fraction == "freem")
DistVar.L.fr = vegdist(veganotu(ps.L.fr), method = "bray")
ps.L.fr.df = data.frame(sample_data(ps.L.fr))

betadisp.L.fr = betadisper(DistVar.L.fr, ps.L.fr.df$plot)

plot = data.frame(betadisp.L.fr$group)
distances = data.frame(betadisp.L.fr$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Lan"
data$fraction= "freem"
data.Lan.fr = data

### Occluded microaggregate, Lancaster
ps.L.oc = subset_samples(ps.L, fraction == "occm")
DistVar.L.oc = vegdist(veganotu(ps.L.oc), method = "bray")
ps.L.oc.df = data.frame(sample_data(ps.L.oc))

betadisp.L.oc = betadisper(DistVar.L.oc, ps.L.oc.df$plot)

plot = data.frame(betadisp.L.oc$group)
distances = data.frame(betadisp.L.oc$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Lan"
data$fraction= "occm"
data.Lan.oc = data

### Bulk soil, Lancaster
ps.L.bu = subset_samples(ps.L, fraction == "fresh")
DistVar.L.bu = vegdist(veganotu(ps.L.bu), method = "bray")
ps.L.bu.df = data.frame(sample_data(ps.L.bu))

betadisp.L.bu = betadisper(DistVar.L.bu, ps.L.bu.df$plot)

plot = data.frame(betadisp.L.bu$group)
distances = data.frame(betadisp.L.bu$distances)
data = cbind(plot,distances)
colnames(data) = c("plot", "distances")

data$site= "Lan"
data$fraction= "fresh"
data.Lan.bu = data

################################## #
#### Boxplots of Distance to Centroid, WITHIN PLOT
################################## #

# Combine the dataframes
data.byplot = rbind(data.Arl.fr, data.Lan.fr, data.Arl.oc, data.Lan.oc, data.Arl.bu, data.Lan.bu)
head(data.byplot)
data.byplot$SampleID = rownames(data.byplot)

# Add tillage treatment, based on plot numbers
meta = read.csv("sample_metadata_Exp3.csv", header=TRUE)
head(meta)
meta = meta [,-c(3,4,5,8,9)]
head(meta)

data.byplot = merge(data.byplot, meta, by = "SampleID")
head(data.byplot)

### Option to save 
#write.csv(data.byplot,"Derived_data/Exp3_Till_BC_Distance_to_centroid_within_plot.csv")

#### Load data, if not continuing from above ####
#data.byplot = read.csv(file="Derived_data/Exp3_Till_BC_Distance_to_centroid_within_plot.csv", row.names=1)
head(data.byplot)

# Re-name variables and put them in order
data.byplot$trt = recode_factor(data.byplot$trt, "NoTill" = "No-tillage", "ConvTill" = "Tillage")
data.byplot$trt = ordered(data.byplot$trt, levels = c("No-tillage", "Tillage"))
data.byplot$fraction = recode_factor(data.byplot$fraction, "fresh" = "Bulk \nsoil", "freem" = "Free \nmicroagg.", "occm" = "Occluded \nmicroagg.")
data.byplot$fraction = ordered(data.byplot$fraction, levels=c("Bulk \nsoil", "Free \nmicroagg.", "Occluded \nmicroagg."))
data.byplot$site = ordered(data.byplot$site, levels=c("Arl", "Lan"))
trtpalette = c("#CC9933", "#117733")

disp.plot = ggplot(subset(data.byplot, site=="Arl"), aes(x = fraction, y = distances, color = trt ))
disp.plot = disp.plot + geom_boxplot()
disp.plot = disp.plot + facet_wrap(~ site, labeller = as_labeller(c("Arl" = "Arlington, within-plot scale \n ",
                                                                    "Lan" = "Lancaster, within-plot scale \n "))) +
  labs(x = "", y = " ", fill = "") +
  scale_color_manual(values = trtpalette, name = "Treatment") +
  scale_fill_manual(name=NULL) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(strip.background = element_rect(fill="white"))+
  ylim(0.14,0.345)
disp.plot # 550 x 370

## Re-run ggplot above separately for each site
## Save the Arlington within-plot figure:
#disp.plot.Arl = disp.plot
## Save the Lancaster within-plot figure:
#disp.plot.Lan = disp.plot


######################## #
#### Distance to Centroid, WITHIN EACH CORE, Arlington ####
######################## #

# (veganotu function was created above)

DistVar.A = vegdist(veganotu(ps.A), method = "bray")
ps.A.df = data.frame(sample_data(ps.A))

betadisp.A = betadisper(DistVar.A, ps.A.df$sample)
betadisp.A
samplecore = data.frame(betadisp.A$group)
distances = data.frame(betadisp.A$distances)
data = cbind(samplecore,distances)
colnames(data) = c("samplecore", "distances")

data$site= "Arl"
data.Arl = data

######################## #
#### Distance to Centroid, WITHIN EACH CORE, Lancaster ####
######################## #

### Lancaster
DistVar.L = vegdist(veganotu(ps.L), method = "bray")
ps.L.df = data.frame(sample_data(ps.L))

betadisp.L = betadisper(DistVar.L, ps.L.df$sample)
betadisp.L
samplecore = data.frame(betadisp.L$group)
distances = data.frame(betadisp.L$distances)
data = cbind(samplecore,distances)
colnames(data) = c("samplecore", "distances")

data$site= "Lan"
data.Lan = data

# Combine the dataframes
data.bycore = rbind(data.Arl, data.Lan)
head(data.bycore)
data.bycore$SampleID = rownames(data.bycore)

# Add tillage treatment, based on plot numbers
meta = read.csv("sample_metadata_Exp3.csv", header=TRUE)
head(meta)
meta = meta [,-c(4,8,9)]

data.bycore = merge(data.bycore, meta, by = "SampleID")

### Option to save 
#write.csv(data.bycore,"Derived_data/Exp3_Till_BC_Distance_to_centroid_within_CORE.csv")

#### Load data, if not continuing from above ####
#data.bycore = read.csv(file="Derived_data/Exp3_Till_BC_Distance_to_centroid_within_CORE.csv", row.names=1)
head(data.bycore)

################################## #
#### Boxplots of Distance to Centroid, WITHIN CORE
################################## #

# Re-name variables and put them in order
data.bycore$trt = recode_factor(data.bycore$trt, "NoTill" = "No-tillage", "ConvTill" = "Tillage")
data.bycore$trt = ordered(data.bycore$trt, levels = c("No-tillage", "Tillage"))
data.bycore$site = ordered(data.bycore$site, levels=c("Arl", "Lan"))
trtpalette = c("#CC9933", "#117733")

disp.core = ggplot(subset(data.bycore, fraction!="fresh" & site == "Lan"), aes(x = trt, y = distances, color = trt ))
disp.core = disp.core + geom_boxplot()
disp.core = disp.core + facet_wrap(~ site, labeller = as_labeller(c("Arl" = "Arlington, within-soil core scale \n(microagg. fractions w/in core)", "Lan" = "Lancaster, within-soil core scale \n(microagg. fractions w/in core)"))) +
  labs(x = "", y = " ", fill = "") +
  scale_color_manual(values = trtpalette,
                     name = "Treatment") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(strip.background = element_rect(fill="white")) +
  ylim(0.14,0.345)
disp.core # 550 x 370

## Re-run ggplot above separately for each site
## Save the Arlington within-core figure:
#disp.core.Arl = disp.core
## Save the Lancaster within-core figure:
#disp.core.Lan = disp.core



######################### #
#### Arrange ALL the Figures together--for Manuscript ####
######################### #

# Uses the between-plot, within-plot, and within-core figs generated above.
g = ggarrange(disp.trt.Arl, disp.plot.Arl, disp.core.Arl,
              disp.trt.Lan, disp.plot.Lan, disp.core.Lan,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2,
          common.legend = TRUE, legend = "bottom")
g
#ggsave("Figures/DistToCentroid_ARLandLAN_Till.tiff", width=7.75, height=6, units = "in", device='tiff', dpi=400)



### Use ANOVA to find significant differences in distance to centroid, by site
### For between-plot-scale ('trt scale') data:
dat.Arl = subset(data.bytrt, site =="Arl")
dat.Lan = subset(data.bytrt, site =="Lan")

ano = aov(distances ~ trt*fraction, data = dat.Lan)
summary(ano)
model.tables(ano, "means")
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))
tuk = TukeyHSD(ano, conf.level = 0.95)
tuk
cld <- multcompLetters4(ano, tuk)
cld

### For within-plot-scale data:
dat.Arl = subset(data.byplot, site =="Arl")
dat.Lan = subset(data.byplot, site =="Lan")

ano = aov(distances ~ trt*fraction, data = dat.Arl)
summary(ano)
model.tables(ano, "means")
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))
tuk = TukeyHSD(ano, conf.level = 0.95)
tuk
cld <- multcompLetters4(ano, tuk)
cld

### For within-core-scale data:
dat.Arl = subset(data.bycore, site =="Arl" & fraction!="fresh")
dat.Lan = subset(data.bycore, site =="Lan" & fraction!="fresh")

ano = aov(distances ~ trt, data = dat.Lan)
summary(ano)
model.tables(ano, "means")
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))
tuk = TukeyHSD(ano, conf.level = 0.95)
tuk
cld <- multcompLetters4(ano, tuk)
cld