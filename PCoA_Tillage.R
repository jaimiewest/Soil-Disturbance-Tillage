library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)
library(ggpubr)


# Read in phyloseq object
ps = readRDS("ps.Exp3")

trtpalette = c("#CC9933", "#117733") #Tillage sites

######################## #
#### PCoA, Arlington ####
######################## #
ps.A = subset_samples(ps,sample_data(ps)$site == "Arl")

# Prune away zero's
ps.A = prune_taxa(taxa_sums(ps.A) > 0, ps.A)

# Hellinger transformation
ps.A = transform_sample_counts(ps.A, function(x) (x / sum(x))^0.5 )

ps.PCoA = ordinate(ps.A, method="PCoA", distance="bray")
p.A = plot_ordination(ps.A, ps.PCoA, color = "trt", shape = "fraction")
p.A = p.A + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.A = p.A + scale_color_manual(values = trtpalette, labels = c("No-tillage", "Tillage"),
                               name = "Treatment")
p.A = p.A + scale_shape_manual(values = c(1, 17, 14), labels = c("Bulk soil", "Free micro", "Occluded micro"),
                               name = "Fraction")

######################## #
#### PCoA, Lancaster ####
######################## #

ps.L = subset_samples(ps,sample_data(ps)$site == "Lan")

# Prune away zero's
ps.L = prune_taxa(taxa_sums(ps.L) > 0, ps.L)

# Hellinger transformation
ps.L = transform_sample_counts(ps.L, function(x) (x / sum(x))^0.5 )

ps.PCoA = ordinate(ps.L, method="PCoA", distance="bray")
p.L = plot_ordination(ps.L, ps.PCoA, color = "trt", shape = "fraction")
p.L = p.L + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.L = p.L + scale_color_manual(values = trtpalette, labels = c("No-tillage", "Tillage"),
                               name = "Treatment")
p.L = p.L + scale_shape_manual(values = c(1, 17, 14), labels = c("Bulk soil", "Free micro", "Occluded micro"),
                               name = "Fraction")
p.L 

######################### #
#### Arrange PCoA  ####
######################### #

ggarrange(p.A, p.L,# bp + rremove("x.text"), 
          labels = c("A", "B"),
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "right")

ggsave("Figures/PCoA.till.tiff", width=7.5, height=4, units = "in", device='tiff', dpi=400)

######################## #
#### Statistics, Arlington ####
######################## #
# First, run PERMANOVA

# # Create veganotu function
# veganotu = function(physeq) {
#   require("vegan")
#   OTU = otu_table(physeq)
#   if (taxa_are_rows(OTU)) {
#     OTU = t(OTU)
#   }
#   return(as(OTU, "matrix"))
# }
ps.A.veg = veganotu(ps.A)
ps.A.df = data.frame(sample_data(ps.A))
DistVar.A = vegdist(ps.A.veg, method = "bray")
adonis2(DistVar.A ~ trt*fraction, data = ps.A.df, method = "bray")

# # Create the pairwise adonis function (JUST BELOW)
## Run the pairwise adonis function for pairwise PERMANOVA
pairwise.adonis(DistVar.A, ps.A.df$fraction)

pairwise.adonis(DistVar.A, ps.A.df$trt)

# # Create the pairwise adonis function
# pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='BH',reduce=NULL,perm=999)
# {
#   co <- combn(unique(as.character(factors)),2)
#   pairs <- c()
#   Df <- c()
#   SumsOfSqs <- c()
#   F.Model <- c()
#   R2 <- c()
#   p.value <- c()
# 
#   for(elem in 1:ncol(co)){
#     if(inherits(x, 'dist')){
#       x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
#                       factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
#     }
# 
#     else  (
#       if (sim.function == 'daisy'){
#         x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
#       }
#       else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
#     )
# 
#     ad <- adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])],
#                  permutations = perm);
#     pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
#     Df <- c(Df,ad$aov.tab[1,1])
#     SumsOfSqs <- c(SumsOfSqs, ad$aov.tab[1,2])
#     F.Model <- c(F.Model,ad$aov.tab[1,4]);
#     R2 <- c(R2,ad$aov.tab[1,5]);
#     p.value <- c(p.value,ad$aov.tab[1,6])
#   }
#   p.adjusted <- p.adjust(p.value,method=p.adjust.m)
# 
#   sig = c(rep('',length(p.adjusted)))
#   sig[p.adjusted <= 0.05] <-'.'
#   sig[p.adjusted <= 0.01] <-'*'
#   sig[p.adjusted <= 0.001] <-'**'
#   sig[p.adjusted <= 0.0001] <-'***'
#   pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
# 
#   if(!is.null(reduce)){
#     pairw.res <- subset (pairw.res, grepl(reduce,pairs))
#     pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
# 
#     sig = c(rep('',length(pairw.res$p.adjusted)))
#     sig[pairw.res$p.adjusted <= 0.1] <-'.'
#     sig[pairw.res$p.adjusted <= 0.05] <-'*'
#     sig[pairw.res$p.adjusted <= 0.01] <-'**'
#     sig[pairw.res$p.adjusted <= 0.001] <-'***'
#     pairw.res <- data.frame(pairw.res[,1:7],sig)
#   }
#   class(pairw.res) <- c("pwadonis", "data.frame")
#   return(pairw.res)
# }

## PERMDISP--Homogeneity of multivariate dispersions. 
#  A test to confirm assumption of homogeneity of variance prior to PERMANOVA.
# Must be done for factors separately.
# For fraction:
betadisp.A = betadisper(DistVar.A, ps.A.df$fraction)
# OR: For tillage trt:
betadisp.A = betadisper(DistVar.A, ps.A.df$trt)

betadisp.A
plot(betadisp.A)
anova(betadisp.A)

#TukeyHSD(betadisp.A)
#boxplot(betadisp.A, xlab = "", las = 2, cex.axis = 0.8)


######################## #
#### Statistics, Lancaster ####
######################## #
# First, run PERMANOVA
ps.L.veg = veganotu(ps.L)
ps.L.df = data.frame(sample_data(ps.L))
DistVar.L = vegdist(ps.L.veg, method = "bray")
adonis2(DistVar.L ~ trt*fraction, data = ps.L.df, method = "bray")

## Run the pairwise adonis function for pairwise PERMANOVA
pairwise.adonis(DistVar.L, ps.L.df$fraction)

pairwise.adonis(DistVar.L, ps.L.df$trt)

## PERMDISP--Homogeneity of multivariate dispersions. 
#  A test to confirm assumption of homogeneity of variance prior to PERMANOVA.
# Must be done for factors separately
# For fraction:
betadisp.L = betadisper(DistVar.L, ps.L.df$fraction)
# For tillage trt:
betadisp.L = betadisper(DistVar.L, ps.L.df$trt)

betadisp.L
plot(betadisp.L)
anova(betadisp.L)

TukeyHSD(betadisp.L)
#boxplot(betadisp.L, xlab = "", las = 2, cex.axis = 0.8)