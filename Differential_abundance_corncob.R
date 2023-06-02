#### Get set up ####
#install.packages("corncob")
library(corncob)
library(phyloseq)
library(dplyr)
library(ggplot2)



# Read in phyloseq object
ps = readRDS("ps.Exp3")
ps.site = subset_samples(ps,sample_data(ps)$site == "Arl" |
                           sample_data(ps)$site == "Lan")
ps.site

ps.site = prune_taxa(taxa_sums(ps.site) > 0, ps.site)
ps.site
ps = ps.site

#### Differential abundance, Tillagee treatments ####

# FIRST! Create empty data frame to hold all results
df = data.frame()

# Subset for one comparison: Within site and fraction
# Loop to run through each site x fraction combination
# Comparison if ConvTill vs the NoTill 'baseline'
for (i in c("Arl", "Lan")){
  ps.subset.site = prune_samples(sample_data(ps)$site == i, ps)
  for (l in c("fresh", "freem", "occm")){
    ps.subset = prune_samples(sample_data(ps.subset.site)$fraction == l, ps.subset.site)
    ps.subset = prune_taxa(taxa_sums(ps.subset)>0, ps.subset)   # Cut out any global zeros within this set

        ## Drop taxa with super low abundance.
    # First, make a list of taxa that are relatively low abundance
    ps.subset.relabund = transform_sample_counts(ps.subset, function(x) x / sum(x))
    AbundTaxa = taxa_names(filter_taxa(ps.subset.relabund, function(x) mean(x) > 0.00001, TRUE))
    # Keep only the abundant taxa
    ps.subset = prune_taxa(AbundTaxa,ps.subset)
    
    # Create post-hoc filters from the relative abundance data
    # We ultimately won't be interested in enriched taxa that are rare even after enrichment
    ps.sub = prune_samples(sample_data(ps.subset.relabund)$site == i, ps.subset.relabund)
    ps.sub = prune_samples(sample_data(ps.sub)$trt == "ConvTill", ps.sub)
    ps.sub = prune_samples(sample_data(ps.sub)$fraction == l, ps.sub)
    RareTrtTaxa = taxa_names(filter_taxa(ps.sub, function(x) mean(x) < 0.001, TRUE))
    
    # We give all parameters of interest (control and variable) to formula and phi.formula,
    # And then drop the parameter we want to test from the _null versions
    # (leaving 1 if there are no control variables, 
    # and the same parameters if we don't want to test for anything (as in phi.formula_null))
    # formula is the differential abundance
    # phi.formula is the differential variance
    # We may just need a very simple model for this dataset, testing for trt
    dT.ps.subset = differentialTest(formula = ~ trt, 
                                 phi.formula = ~ trt,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ trt,
                                 test = "Wald", boot = FALSE,
                                 data = ps.subset,
                                 fdr_cutoff = 0.05)
  
    # Making an empty dataframe to hold the full results
    df.ps  = data.frame()
    
    # Loop to pull out coefficients for each taxon
    for (j in 1:length(dT.ps.subset$significant_taxa)){
      # Get the significant model for that taxon
      sig_models = dT.ps.subset$significant_models[[j]]
      # Pull out the coefficients as above
      mu = data.frame(t(as.matrix(sig_models$coefficients[2,])))
      # Also grab the p_fdr estimate for that taxon's model
      p_fdr = dT.ps.subset$p_fdr[dT.ps.subset$significant_taxa][j]
      # Add that estimate onto our coefficient data frame
      mu$p_fdr = p_fdr
      # Create a column with the OTU ID
      mu$OTU= paste(row.names(data.frame(p_fdr)))
      # Add this row onto the df dataframe, which will collect the results
      # for all taxa as it iterates through this loop.
      df.ps = rbind(df.ps,mu)
    }

    # Make the column names nicer
    colnames(df.ps) = c("Estimate","SE","t","p","p_fdr","OTU")
    
    # Let's bring back in our taxonomy from the tax table
    SigOTUs = levels(as.factor(df.ps$OTU))
    pruned = prune_taxa(SigOTUs,ps.subset)
    taxtab = data.frame(tax_table(pruned))
    taxtab$OTU = c(taxa_names(pruned))
    joined.ps.subset = merge(df.ps,taxtab,by=c("OTU"))
    
    # Make column to designate if OTU is rare, for filtering later
    joined.ps.subset$RareTrtTaxa = ifelse(joined.ps.subset$OTU %in% RareTrtTaxa, "rare", "not rare")
    
    # Prep for merging with other treatments/fractions
    joined.ps.subset$site = paste(i)
    joined.ps.subset$fraction = paste(l)
    
    
    # Make final dataframe by joining together each differential test set
    df = rbind(df, joined.ps.subset)
  }
}

dim(df)
head(df)


# Fix up some taxon naming issues
ignoreList = c("","uncultured","uncultured bacterium","uncultured soil bacterium","uncultured forest soil bacterium", "Ambiguous_taxa",
               "uncultured actinobacterium","uncultured planctomycete","uncultured Chloroflexi bacterium","uncultured Syntrophobacteraceae bacterium",
               "uncultured crenarchaeote","uncultured Gemmatimonadetes bacterium", "uncultured Acidobacteria bacterium",
               "uncultured Planctomyces sp.", "metagenome", "Subgroup 6", "Blastocatellia (Subgroup 4)", "uncultured Holophaga sp.",
               "uncultured Hyphomicrobiaceae bacterium", "uncultured proteobacterium", "1921-2", "RB41", "A21b",
               "AD3", "Subgroup_7", "Subgroup_2", "WD260", "Subgroup_17", "Subgroup_10", "A4b", "SBR1031",
               "MB-A2-108", "KD4-96","Tychonema_CCAP_1459-11B","Nodosilinea_PCC-7104","Microcoleus_PCC-7113","SBR1031",
               "BIrii41", "JG30-KF-AS9", "67-14", "SC-I-84", "GOUTA6", "IS-44", "966-1", "B1-7BS", "MND1", "JG30-KF-CM45")

df = df %>%
  mutate(Name = ifelse(Genus %in% ignoreList | is.na(Genus),
                       ifelse(Family %in% ignoreList | is.na(Family),
                              ifelse(Class %in% ignoreList | is.na(Class),
                                     ifelse(Order %in% ignoreList | is.na(Order),
                                     ifelse(Phylum %in% ignoreList | is.na(Phylum), paste(OTU),
                                            paste(Phylum)),paste(Class)),paste(Order)),paste(Family)),paste(Genus))) %>%
  mutate(Name = ifelse(Name == "uncultured thaumarchaeote","Thaumarchaeota",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Xiphinematobacter","Xiphinematobacter",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus Solibacter","Solibacter",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Alysiosphaera","Alysiosphaera",Name)) %>%
  mutate(Name = ifelse(Name == "NB1-j","Myxococcales",Name)) %>%
  mutate(Name = ifelse(Name == "IMCC26256","Acidimicrobiia",Name)) %>%
  mutate(Name = ifelse(Name == "SBR1031","Anaerolineae",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Udaeobacter","Udaeobacter",Name))


# May want to save at this point
#write.csv(df,"Derived_data/Differential_Abundance_Tillage_trt_differences_0.00001_0.001.csv", row.names = FALSE)
# 
df <- read.csv(file = 'Derived_data/Differential_Abundance_Tillage_trt_differences_0.00001_0.001.csv')
head(df)

# How many total responders (> 0), or big responders (>1), that aren't overly rare?
total.OTUs=df%>%
  #filter(site == "Arl")%>%
  filter(Estimate > 1)%>%
  filter(RareTrtTaxa=="not rare")#%>%
  #count(fraction)%>% # for number of OTUs, by trt
  #summarise_each(funs(n_distinct)) # for total number of OTUs
total.OTUs

#### Make a table for SI of all the enriched/depleted OTUs
# First, need to grab all the sequences from the .fasta file
library("Biostrings")
Exp3Seq = readDNAStringSet("dna-sequences_Exp3.fasta")
OTU = names(Exp3Seq)
sequence = paste(Exp3Seq)
Exp3Seq.df <- data.frame(OTU, sequence)

# Then, add these seqs to the dataframe of taxa.
df2 = merge(total.OTUs, Exp3Seq.df, by = "OTU")
head(df2)

# Write to .csv
write.csv(df2,"Derived_data/Tillage_enriched_taxa_with_sequences.csv", row.names = FALSE)

#write.csv(df2,"Derived_data/Tillage_ALL_enriched_taxa_>0_withRare_with_sequences.csv", row.names = FALSE)


# Tons of taxa, will probably want to filter somewhat.
# Looking at biggest positive responders
Big.pos = df %>% filter(Estimate > 1)
dim(Big.pos)
# Remove the responders, who even after enrichment, are still 'rare'
Big.pos = Big.pos %>% filter(RareTrtTaxa=="not rare")
dim(Big.pos) # includes duplicates e.g. across fractions or sites
head(Big.pos)

# Re-name variables and put them in order
Big.pos$fraction = recode_factor(Big.pos$fraction,
                                      "fresh" = "Bulk soil",
                                      "freem" = "Free microaggregate",
                                      "occm" = "Occluded microaggregate")
Big.pos$site = recode_factor(Big.pos$site,
                                  Arl = "Arlington",
                                  Lan = "Lancaster")

p = ggplot(Big.pos,aes(x=Estimate, color=Phylum, y=reorder(OTU, Estimate)))
p = p + theme_bw() + theme(#panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  strip.background = element_rect(fill="white"))
p = p + geom_point(size=2)
p = p + geom_errorbar(aes(xmax=Estimate+1.96*SE,xmin=Estimate-1.96*SE)) #plot negative estimate for y-value
#p = p + theme(axis.text.x=element_text(angle=45,face="italic", size=6, hjust=1, vjust=1))
p = p + xlab("Coefficient of differential abundance for taxa enriched under tillage") + ylab("")
p = p + scale_y_discrete(breaks=Big.pos$OTU,labels=Big.pos$Name)
p = p + facet_grid(vars(site), vars(fraction), scales="free_y", space="free_y")
p = p + scale_shape_manual(values = c(1, 17, 14), name = "Soil fraction")
#p = p + geom_hline(xintercept=0)
p = p + theme(legend.position = "bottom")
p #720 x 650

#ggsave("Figures/DiffAbund_pos_Tillage.tiff", width=7.5, height=7.5, units = "in", device='tiff', dpi=400)


#### Negative responders/Depleted taxa i.e., the taxa that are enriched under NO-TILLAGE####
# Retain only biggest negative responders
Big.neg = df %>% filter(Estimate < -1)
dim(Big.neg)
# Remove the responders that are 'rare'
Big.neg = Big.neg %>% filter(RareTrtTaxa=="not rare")
dim(Big.neg)
Big.neg

# Re-name variables and put them in order
Big.neg$fraction = recode_factor(Big.neg$fraction,
                                 "fresh" = "Bulk soil",
                                 "freem" = "Free microaggregate",
                                 "occm" = "Occluded microaggregate")
Big.pos$site = recode_factor(Big.pos$site,
                             Arl = "Arlington",
                             Lan = "Lancaster")

# Here we will use the negative of the estimate in order to call this enriched under no-tillage (the control)
# as opposed to calling it taxa depleted under tillage

p = ggplot(Big.neg, aes(x=-Estimate, color=Phylum, y= reorder(OTU, -Estimate)))
p = p + theme_bw() + theme(#panel.grid.major = element_blank(), 
                                              panel.grid.minor = element_blank(), 
                                              strip.background = element_rect(fill="white"))
p = p + geom_point(size=2)
p = p + geom_errorbar(aes(xmax=-Estimate+1.96*SE,xmin=-Estimate-1.96*SE)) #plot negative estimate for y-value
#p = p + theme(axis.text.x=element_text(angle=45,face="italic", size=6, hjust=1, vjust=1))
p = p + xlab("Coefficient of differential abundance for taxa enriched under no-tillage") + ylab("")
p = p + scale_y_discrete(breaks=Big.neg$OTU,labels=Big.neg$Name)
p = p + facet_grid(vars(site), vars(fraction), scales="free_y", space="free_y")
#p = p + geom_hline(yintercept=0)
p = p + theme(legend.position = "bottom")
p # 675 x 500

#ggsave("Figures/DiffAbund_neg.Tillage.tiff", width=6.75, height=5, units = "in", device='tiff', dpi=400)


################ So, we have a list of the significant responders
# Note - corncob gives us an estimate in a log-likelihood model that estimates relative abundance.
# "Estimate" is not "log2-fold change". If we wanted to plot fold-change, we could derive it.