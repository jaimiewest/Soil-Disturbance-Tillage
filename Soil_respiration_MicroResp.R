library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)


CO2 = read.csv(file ="MicroResp_CO2_forR.csv", header=TRUE, stringsAsFactors=TRUE, na.strings = "NA")
row.names(CO2) = CO2$sample

# Choose tillage sites
CO2 = subset(CO2, site == "Arl" | site == "La")

head(CO2)

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


# Order fractions for graphing
CO2meanSE$trt = recode_factor(CO2meanSE$trt,
                      NoTill = "No-tillage", 
                      ConvTill = "Tillage")
CO2meanSE$trt = ordered(CO2meanSE$trt, levels = c("No-tillage", "Tillage"))
CO2meanSE$site = ordered(CO2meanSE$site, levels=c("Arl", "Lan"))

trtpalette = c("#CC9933", "#117733")
site.labs=c("Arl" = "Arlington", "Lan" = "Lancaster")

# ### Graph: CO2 respiration rate--point with error bars
# dodge = position_dodge(0.5)
# p = ggplot(subset(CO2meanSE, element == "CO2 rate, per g soil C"), aes(x=factor(trt), y=mean, color=trt)) +
#   geom_point(size=2, position = dodge) +
#   geom_errorbar(aes(ymin=mean-1.96*SE,ymax=mean+1.96*SE), position = dodge, width=0.2)
# p = p + facet_grid(~site, labeller=labeller(site=site.labs))
# p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
# p = p + ylab(expression(CO[2]*", \u03BCg C/g soil C/hour")) # + ylim(1250,2500)
# p = p + xlab("Treatment")
# p = p + scale_color_manual(values=trtpalette))
# p = p + theme(legend.title = element_blank()) #+ theme(legend.position = c(0.15, 0.15), legend.box = "horizontal")
# p = p + theme(axis.text.x = element_text(angle = 45))
# p = p + theme(legend.title= element_blank())
# p # 450x300


# Order fractions for graphing
CO2$trt = recode_factor(CO2$trt,
                        NoTill = "No-tillage", 
                        ConvTill = "Tillage")
CO2$trt = ordered(CO2$trt, levels = c("No-tillage", "Tillage"))
CO2$site = ordered(CO2$site, levels=c("Arl", "Lan"))
trtpalette = c("#CC9933", "#117733")
site.labs=c("Arl" = "Arlington", "Lan" = "Lancaster")
head(CO2)
### Graph: CO2 respiration rate--BOXPLOTS
#dodge = position_dodge(0.5)
#p = ggplot(CO2, aes(x=trt, y=mean_CO2rate_ugC_perg_perh, color=trt)) +
p = ggplot(CO2, aes(x=trt, y=mean_CO2rate_ugC_per_gC_perh, color=trt)) +
  geom_boxplot() +
  #geom_jitter(alpha=0.3) +
  expand_limits(y=c(0,1))
p = p + facet_grid(~site, labeller=labeller(site=site.labs))
p = p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p = p + ylab(expression(CO[2]*", \u03BCg C/g soil C/hour")) 
#p = p + ylab(expression(CO[2]*", \u03BCg C/g soil/hour"))
p = p + xlab(element_blank())
p = p + scale_color_manual(values=trtpalette)
p = p + theme(legend.title = element_blank()) + theme(legend.position = "none")#c(0.15, 0.15), legend.box = "horizontal")
p = p + scale_x_discrete(labels=c("NoTill" = "No-tillage", "ConvTill" = "Tillage"))
p = p + theme(legend.title= element_blank())
p # 450x300

#ggsave("CO2.pergsoil.till.tiff", width=3.5, height=3.2, units = "in", device='tiff', dpi=400)
#ggsave("CO2.pergsoil_C.till.tiff", width=3.5, height=3.2, units = "in", device='tiff', dpi=400)




#### ANOVA, repeat for each site and analysis ####
head(CO2)
dat = subset(CO2, site == "Arl")
#dat = subset(CO2, site == "Lan")

ano = aov(mean_CO2rate_ugC_per_gC_perh ~ trt, data = dat)
#ano = aov(mean_CO2rate_ugC_perg_perh ~ trt, data = dat)

summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

# tuk = TukeyHSD(ano, conf.level = 0.95)
# tuk
