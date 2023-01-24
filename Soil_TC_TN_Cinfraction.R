library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(rcartocolor)
library(ggpattern)
library(multcompView) # for compact letter display


# Read in data
TCTN = read.csv(file ="TC_TN_forR.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA")
# Subset for site(s)
TCTN = subset(TCTN, site=="Arl" | site == "Lan")
head(TCTN)

# Find mean
TCTNmean = TCTN%>%
  group_by(site, trt, fraction) %>% 
  dplyr::summarize(n=n(),
                   "Total carbon, percent" = mean(C_pct, na.rm=TRUE),
                   "Total nitrogen, percent" = mean(N_pct, na.rm=TRUE),
                   "C:N ratio" = mean(C_N_ratio, na.rm=TRUE)) %>%
  gather("element", "mean", - c(site, trt, fraction, n), factor_key=TRUE)  

# Find SE
TCTNSE = TCTN%>%
  group_by(site, trt, fraction) %>% 
  dplyr::summarize(n=n(),
                   "Total carbon, percent" = sd(C_pct, na.rm=TRUE)/sqrt(n()),
                   "Total nitrogen, percent" = sd(N_pct, na.rm=TRUE)/sqrt(n()),
                   "C:N ratio" = sd(C_N_ratio, na.rm=TRUE)/sqrt(n())) %>%
  gather("element", "SE", - c(site, trt, fraction, n), factor_key=TRUE)  

# Join mean proportion and SE data together
TCTNmeanSE = merge(TCTNmean,TCTNSE, by=c("site","trt", "fraction", "n", "element"))      
head(TCTNmeanSE)


# Order fractions for graphing
TCTNmeanSE$fraction = ordered(TCTNmeanSE$fraction, levels=c("fresh", "mac", "freem", "occm"))
TCTNmeanSE$trt = ordered(TCTNmeanSE$trt, levels=c("NoTill", "ConvTill"))
trtpalette = c("#CC9933", "#117733")
TCTNmeanSE$site = ordered(TCTNmeanSE$site, levels=c("Arl", "Lan"))
site.labs=c("Arl" = "Arlington", "Lan" = "Lancaster")

# Add a column designating occluded fractions. This will be used to hash (or pattern) the occluded fracitons.
TCTNmeanSE = TCTNmeanSE %>%
  mutate(occluded = case_when(fraction == "occm" ~ "occluded",
                              fraction == "fresh" | fraction == "freem" | fraction == "mac" ~ "not"))


### Graph: Total C, BAR graph--SI FIGURE.

p = ggplot(subset(TCTNmeanSE, element == "Total carbon, percent"), aes(x = fraction, y = mean, fill=trt, pattern = occluded)) +
  geom_bar_pattern(stat = "identity", color="black", position=position_dodge(),
                   pattern_fill = "lightgray",
                   pattern_angle = 45,
                   pattern_density = 0.01,
                   pattern_spacing = 0.04
  ) +
  guides(fill = guide_legend(override.aes=list(pattern="none"))) +
  scale_pattern_manual(values=c(occluded="stripe", not="none"), guide = "none") +
  geom_errorbar(aes(ymin=mean-1.96*SE,ymax=mean+1.96*SE), position=position_dodge(0.9),width=0.15)
p = p + facet_grid(~site, labeller=labeller(site=site.labs))
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab("Total carbon (%)") # + ylim(1250,2500)
p = p + xlab("Soil fraction")
p = p + scale_x_discrete(labels=c("fresh"="Bulk \nsoil","mac"= "Macroaggregate", "freem"="\nFree \nmicroaggregate", "occm"="Occluded \nmicroaggregate"))
p = p + scale_fill_manual(values=trtpalette, labels = c("No-tillage", "Tillage"))
p = p + theme(legend.title = element_blank()) + theme(legend.position = c(0.1, 0.86), legend.box = "horizontal")
p = p + theme(legend.key = element_rect(color = "black"))
p # 750 x 350 for SI


### Graph: Total N, BAR graph--not published
p = ggplot(subset(TCTNmeanSE, element == "Total nitrogen, percent"), aes(x = fraction, y = mean, fill=trt)) +
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-1.96*SE,ymax=mean+1.96*SE), position=position_dodge(0.9),width=0.15)
p = p + facet_grid(~site, labeller=labeller(site=site.labs))
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab("Total nitrogen (%)") # + ylim(1250,2500)
p = p + xlab("Soil fraction")
p = p + scale_x_discrete(labels=c("fresh"="Bulk \nsoil","mac"= "Macro", "freem"="Free \nmicro", "occm"="Occluded \nmicro"))
p = p + scale_fill_manual(values=trtpalette, labels = c(#"No worm", "Low worm", "High worm",
  "No-tillage", "Tillage"))
p = p + theme(legend.title = element_blank()) + theme(legend.position = c(0.1, 0.86), legend.box = "horizontal")
p

### Graph: boxplots, C:N ratio

# Order fractions for graphing
TCTN$fraction = recode_factor(TCTN$fraction, fresh = "Bulk \nsoil", mac = "Macro",
                              freem  ="Free \nmicro.", occm = "Occluded \nmicro.")
TCTN$fraction = ordered(TCTN$fraction, levels=c("Bulk \nsoil", "Macro", "Free \nmicro.", "Occluded \nmicro."))
TCTN$trt = recode_factor(TCTN$trt, NoTill = "No-tillage", ConvTill = "Tillage")
TCTN$trt = ordered(TCTN$trt, levels=c("No-tillage", "Tillage"))
TCTN$site = recode_factor(TCTN$site, Arl = "Arlington", Lan = "Lancaster")

TCTN
p = ggplot(TCTN, aes(x = fraction, y = C_N_ratio, color=trt)) + 
  geom_boxplot()
p = p + facet_grid(~site)
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab("C:N ratio")
p = p + xlab("Soil fraction")
p = p + scale_color_manual(values=trtpalette)
p = p + theme(legend.title = element_blank()) + theme(legend.position = c(0.1, 0.86), legend.box = "horizontal")
#p = p + theme(axis.text.x = element_text(angle = 45))
p # 750 x 350 for SI

#### ANOVA, by site ####
# Arlington, TC
dat = subset(TCTN, site == "Arl")
ano = aov(C_pct ~ trt*fraction, data = dat)
summary(ano)
model.tables(ano, "means", se=TRUE)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))
tuk = TukeyHSD(ano, conf.level = 0.95)
tuk
cld = multcompLetters4(ano, tuk)
cld

# Arlington, TN
dat = subset(TCTN, site == "Arl")
ano = aov(N_pct ~ trt*fraction, data = dat)
summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))
tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

# Lancaster, TC
dat = subset(TCTN, site == "Lan")
ano = aov(C_pct ~ trt*fraction, data = dat)
summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))
tuk = TukeyHSD(ano, conf.level = 0.95)
tuk

# Lancaster, TN
dat = subset(TCTN, site == "Lan")
ano = aov(N_pct ~ trt*fraction, data = dat)
summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))
tuk = TukeyHSD(ano, conf.level = 0.95)
tuk


#### Graph of total C per g of soil by fraction--For manuscript figrue
# Need to first create agg mean

#Refresh data, and clean up dataframes prior to merging
TCTN = read.csv(file ="TC_TN_forR.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA")
head(TCTN)
TCTN2 = TCTN[, -c(1,8,9,11)]
TCTN2 = spread(TCTN2, fraction, C_pct)
head(TCTN2)

agg = read.csv("aggregate_fractions_proportions.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA", row.names=1)
head(agg)
agg = agg[, -c(13,14)]

TCTNagg = merge(TCTN2,agg, by=c("sample", "site","trt"))      
head(TCTNagg)
TCTNagg = subset(TCTNagg, site == "Arl" | site == "Lan")

# multiply out C percent and proportion of soil in each fraction.
# The column headings here are not super clena, but it'll be ok.
# Double check multiplier--depending on how original data is presented (proportion vs percent etc)
TCTNagg$freem_C = TCTNagg$freem.x * TCTNagg$freem.y *10
TCTNagg$occm_C = TCTNagg$occm * TCTNagg$occm.Mac* TCTNagg$Mac *10
TCTNagg$mac_C = TCTNagg$mac *  TCTNagg$Mac *10

# Option to add in estimates for silt + clay fractions
TCTNagg$SC_C = (TCTNagg$fresh*10-TCTNagg$mac_C-TCTNagg$freem_C)
TCTNagg$occSC_C = (TCTNagg$mac_C-TCTNagg$occm_C)


levels(TCTNmeanSE$site)
# Find proportion of each aggregate
TCTNaggmean = TCTNagg %>%
  group_by(site, trt) %>% 
  dplyr::summarize(n=n(),
                   "freem" = mean(freem_C, na.rm=TRUE),
                   "mac" = mean(mac_C, na.rm = TRUE),
                   "occm" = mean(occm_C, na.rm = TRUE),
                   "fresh" = mean((fresh*10), na.rm = TRUE),
                   "SC" = mean(SC_C, na.rm=TRUE),
                   "occSC" = mean(occSC_C, na.rm=TRUE)) %>%
  gather("fraction", "agg_fraction_C", - c(site, trt, n), factor_key=TRUE)  
TCTNaggmean
str(TCTNaggmean)

# Find SE
TCTNaggSE = subset(TCTNagg) %>%
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
TCTNaggSE

# Join mean proportion and SE data together
aggmeanSE = merge(TCTNaggmean,TCTNaggSE, by=c("site","trt", "n", "fraction"))      
str(aggmeanSE)

# Find the proportion of total C (in bulk soil/"fresh") in each fraction
aggmeanSE = aggmeanSE %>%
  group_by(site, trt) %>%
  mutate(prop = agg_fraction_C / agg_fraction_C[match('fresh', fraction)])
aggmeanSE$prop.SE = aggmeanSE$SE*aggmeanSE$prop/aggmeanSE$agg_fraction_C
View(aggmeanSE)

# Add column designating occluded fractions for shading/patterned bar in figure
aggmeanSE = aggmeanSE %>%
  mutate(occluded = case_when(fraction == "occm" | fraction == "occSC" ~ "occluded",
                              fraction == "fresh" | fraction == "freem" | fraction == "mac" | fraction == "SC"~ "not"))

# Option to save/read in saved data
#write.csv(aggmeanSE,"Carbon_in_aggregate_fractions.csv")
#aggmeanSE = read.csv("Carbon_in_aggregate_fractions.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA", row.names=1)


### Graph: amount of C in each soil fraction for MS
# Order fractions for graphing
aggmeanSE$fraction = ordered(aggmeanSE$fraction, levels=c("fresh", "mac", "freem", "SC", "occm", "occSC"))
aggmeanSE$trt = ordered(aggmeanSE$trt, levels=c("NoTill", "ConvTill"))
trtpalette = c("#CC9933", "#117733")
site.labs=c("Arl" = "Arlington", "Lan" = "Lancaster")
aggmeanSE$trt = recode_factor(aggmeanSE$trt, NoTill = "No-tillage", ConvTill = "Tillage")

p = ggplot(aggmeanSE, aes(x = fraction, y = agg_fraction_C, fill=trt, pattern = occluded)) +
  geom_bar_pattern(stat = "identity", color="black", position=position_dodge(),
                   pattern_fill = "black",
                   pattern_colour = "#444444",
                   pattern_angle = 45,
                   pattern_density = 0.01,
                   pattern_spacing = 0.04) +
  guides(fill = guide_legend(override.aes=list(pattern="none"))) +
  scale_pattern_manual(values=c(occluded="stripe", not="none"), guide = "none")
p = p + geom_errorbar(data=subset(aggmeanSE, fraction!="SC" & fraction!="occSC"), #no error bars for the estimates of SC fraction
                aes(ymin=agg_fraction_C-1.96*SE,ymax=agg_fraction_C+1.96*SE), position=position_dodge(0.9),width=0.15)
p = p + facet_grid(~site, labeller=labeller(site=site.labs))
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab("mg C (in fraction) per g bulk soil") + ylim(-0.5,38)
p = p + xlab("Soil fraction")
p = p + scale_x_discrete(labels=c("fresh"="Bulk \nsoil","mac"= "Macroagg.", "freem"="Free \nmicroagg.", "SC"="Silt + clay \n(estimate)",
                                  "occm"="Occluded \nmicroagg.", "occSC"="Occluded \nsilt + clay \n& POM \n(estimate)"))
p = p + scale_fill_manual(values=trtpalette)
p = p + theme(legend.title = element_blank()) + theme(legend.position = c(0.9, 0.83), legend.box = "horizontal")
p = p + theme(legend.key = element_rect(color = "black"))
p # 750 x 350

ggsave("CinFraction.till.tiff", width=7.7, height=3.7, units = "in", device='tiff', dpi=400)

head(TCTNagg)
# Statistics, treatment differences within fraction
# Must update site for different sites
# and must update the fraction in the aov statement
dat = subset(TCTNagg, site == "Lan")
ano = aov(occm_C ~ trt, data = dat)
summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

tuk = TukeyHSD(ano, conf.level = 0.95)
tuk


### Graph: proportion of C in fraction
# Order fractions for graphing
aggmeanSE$fraction = ordered(aggmeanSE$fraction, levels=c("fresh", "mac", "freem", "SC", "occm", "occSC"))
aggmeanSE$trt = ordered(aggmeanSE$trt, levels=c("NoTill", "ConvTill")) #"NoWorm", "LowWorm", "HighWorm",
trtpalette = c("#CC9933", "#117733") # "#FF6666","#CC0000","#661100", 
site.labs=c("Arl" = "Arlington", "Lan" = "Lancaster")
aggmeanSE$trt = recode_factor(aggmeanSE$trt, NoWorm = "No Amynthas",
                              LowWorm = "No Amynthas",
                              HighWorm  ="High Amynthas",
                              NoTill = "No-tillage", 
                              ConvTill = "Tillage")

p = ggplot(aggmeanSE, aes(x = fraction, y = prop, fill=trt, pattern = occluded)) +
  geom_bar_pattern(stat = "identity", color="black", position=position_dodge(),
           pattern_fill = "lightgray",
           pattern_angle = 45,
           pattern_density = 0.01,
           pattern_spacing = 0.04) +
  guides(fill = guide_legend(override.aes=list(pattern="none"))) +
  scale_pattern_manual(values=c(occluded="stripe", not="none"), guide = "none")
p = p + geom_errorbar(data=subset(aggmeanSE, fraction!="SC" & fraction!="occSC"), #no error bars for the estimates of SC fraction
                aes(ymin=prop-1.96*prop.SE,ymax=prop+1.96*prop.SE), position=position_dodge(0.9),width=0.15)
p = p + facet_grid(~site, labeller=labeller(site=site.labs))
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab("Proportion of total C in soil fraction") #+ ylim(-1,40)
p = p + xlab("Soil fraction")
p = p + scale_x_discrete(labels=c("fresh"="Bulk \nsoil","mac"= "Macroagg.", "freem"="Free \nmicroagg.", "SC"="Silt + clay \n(estimate)",
                                  "occm"="Occluded \nmicroagg.", "occSC"="Occluded \nsilt + clay \n& POM \n(estimate)"))
p = p + scale_fill_manual(values=trtpalette)
p = p + theme(legend.title = element_blank()) + theme(legend.position = c(0.9, 0.83), legend.box = "horizontal")
p # 770 x 370 for SI


ggsave("CinFraction.Proportion of whole.till.tiff", width=7.5, height=3.5, units = "in", device='tiff', dpi=400)