library(doMC)
library(BiocManager)
library(phyloseq)
library(picante)

### This code is likely best used on a high-throughput computing cluster.
### This code is based on Stegen et al., 2013
### found at https://github.com/stegen/Stegen_etal_ISME_2013

# Bring in ps object
ps.TRT = readRDS("ps.Exp3")
ps.TRT

# Site Subset (run separately for each site)
ps.TRT.site = subset_samples(ps.TRT, site=="Arl")
#ps.TRT.site = subset_samples(ps.TRT, site=="Lan")

# Hellinger transformation
ps.TRT.site = transform_sample_counts(ps.TRT.site, function(x) (x / sum(x))^0.5 )

# Prune away zero's
ps.TRT.site = prune_taxa(taxa_sums(ps.TRT.site) > 0, ps.TRT.site)
otu=otu_table(ps.TRT.site)

#head(otu)
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and 5 columns


## read in the phylogeny
phylo = read.tree("tree_Exp3.nwk");
phylo; # a summary of the phylogeny

## make sure the names on the phylogeny are ordered the same as the names in otu table
## and drop unused tips
match.phylo.otu = match.phylo.data(phylo, t(otu));

## calculate empirical betaMNTD
beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];

write.csv(beta.mntd.weighted,'betaMNTD_weighted.csv',quote=F);

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD
beta.reps = 999; # number of randomizations. Recommendation 999-9999

rand.weighted.bMNTD.comp = array(c(-100),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

# Set number of computing threads or cores
registerDoMC(cores = 16)

rand.weighted.bMNTD.comp = foreach (rep = 1:beta.reps) %dopar% {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
}

rand.weighted.bMNTD.comp.array = array(as.numeric(unlist(rand.weighted.bMNTD.comp)),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps))
rownames(rand.weighted.bMNTD.comp.array)=rownames(rand.weighted.bMNTD.comp[[1]])
colnames(rand.weighted.bMNTD.comp.array)=colnames(rand.weighted.bMNTD.comp[[1]])
head(rand.weighted.bMNTD.comp.array)

rand.weighted.bMNTD.comp=rand.weighted.bMNTD.comp.array

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
head(weighted.bNTI);
write.csv(weighted.bNTI,"weighted_bNTI.csv",quote=F);

pdf("weighted_bNTI_Histogram.pdf")
  hist(weighted.bNTI)
dev.off()

