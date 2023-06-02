
library(doMC)
library(BiocManager)
library(phyloseq)
library(picante)
library(vegan)

### This code is likely best used on a high-throughput computing cluster.
### This code is based on Stegen et al., 2013
### found at https://github.com/stegen/Stegen_etal_ISME_2013

# Bring in ps object
ps.TRT = readRDS("ps.Exp3")

# Site Subset (run separately for each site)
ps.TRT.site=subset_samples(ps.TRT, site=="Arl")
#ps.TRT.site=subset_samples(ps.TRT, site=="Lan")

# Hellinger transformation
ps.TRT.norm = transform_sample_counts(ps.TRT.site, function(x) (x / sum(x))^0.5 )

# Prune away zero's
ps.TRT.norm = prune_taxa(taxa_sums(ps.TRT.norm) > 0, ps.TRT.norm)

# Turn into integers
ComSize = max(sample_sums(ps.TRT.site)) #used in for integer-izing
ps.TRT.norm.int = transform_sample_counts(ps.TRT.norm, function(x) round((x*ComSize), digits=10))
sample_sums(ps.TRT.norm.int) #should all be equal

# read in OTU table
otu=as.matrix(otu_table(ps.TRT.norm.int))
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and 5 columns

# Need this to fix the fact that you can't exactly have half of a read created from my averaging step
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

rc.reps = 9999 #recommendation is 9999

# Set number of computing threads or cores.
registerDoMC(cores = 32)

raup_crick_abundance = function(spXsite, plot_names_in_col1=FALSE, classic_metric=FALSE,
                                split_ties=TRUE, reps=rc.reps, set_all_species_equal=FALSE, 
                                as.distance.matrix=TRUE, report_similarity=FALSE){
	
	##expects a species by site matrix for spXsite, with row names for plots, 
  #or optionally plots named in column 1.  By default calculates a modification of the 
  #Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). 
  #Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges 
  #from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null 
  #observations that are equal to the observed number of shared species to the calculation- 
  #this is highly recommended.  The argument report_similarity defaults to FALSE so the function 
  #reports a dissimilarity (which is appropriate as a measure of beta diversity).  
  #Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  
  #If ties are split (as we recommend) the dissimilarity (default) and similarity 
  #(set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, 
  #which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges 
  #from 0 to 1). If ties are not split (and there are ties between the observed and expected shared 
  #number of species) this conversion will not work. The argument reps specifies the number of 
  #randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights 
  #all species equally in the null model instead of weighting species by frequency of occupancy.  
	
	
	##Note that the choice of how many plots (rows) to include has a real impact on the metric, 
  #as species and their occurrence frequencies across the set of plots is used to determine gamma 
  #and the frequency with which each species is drawn from the null model	
  
  
	##this section moves plot names in column 1 (if specified as being present) into the row names 
  #of the matrix and drops the column of names
  
	if(plot_names_in_col1){
		row.names(spXsite)<-spXsite[,1]
		spXsite<-spXsite[,-1]
		}
	
	
	## count number of sites and total species richness across all plots (gamma)
	n_sites<-nrow(spXsite)
	gamma<-ncol(spXsite)

	##build a site by site matrix for the results, with the names of the sites in the row and col names:
	results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
	
	##make the spXsite matrix into a new, pres/abs. matrix:
	ceiling(spXsite/max(spXsite))->spXsite.inc
	
	##create an occurrence vector- used to give more weight to widely distributed species in the null model:
	occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
	
	##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
	abundance<-apply(spXsite, MARGIN=2, FUN=sum)
	
	##make_null:
	
	##looping over each pairwise community combination, generating null results:
	
	for(null.one in 1:(nrow(spXsite)-1)){
	  for(null.two in (null.one+1):nrow(spXsite)){
	    
	    null_bray_curtis<-NULL
	    null_bray_curtis = foreach(i=1:rc.reps, .packages = c("vegan","picante")) %dopar% {
				
				##two empty null communities of size gamma:
				com1<-rep(0,gamma)
				com2<-rep(0,gamma)
				
				##add observed number of species to com1, weighting by species occurrence frequencies:
				com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
				com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
				com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
				com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
				com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
				com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
				#sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
				rm('com1.samp.sp','com1.sp.counts');			
				
				##same for com2:
				com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
				com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
				com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
				com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
				com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
				com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
				# sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
				rm('com2.samp.sp','com2.sp.counts');

				null.spXsite = rbind(com1,com2); # null.spXsite;
				
				##calculate null bray curtis
				null_bray_curtis[i] = distance(otu_table(null.spXsite,taxa_are_rows = FALSE),
				                               method='bray');
				
				### OR----Calculates the null Bray-Curtis as it is on Danczakre github:
				#null_bray_curtis[i] = vegdist(null.spXsite, method='bray');
				
			}; # end reps loop

			## empirically observed bray curtis
			obs.bray = distance(otu_table(spXsite[c(null.one,null.two),], taxa_are_rows = FALSE), 
			                    method='bray');

			##how many null observations is the observed value tied with?
			num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
			
			##how many null values are smaller than the observed *dissimilarity*?
			num_less_than_in_null = sum(null_bray_curtis<obs.bray);
			rc = (num_less_than_in_null )/reps; # rc;
			if(split_ties){
				rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
			};
			
			
			if(!classic_metric){
					##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
					rc = (rc-.5)*2
			};

			results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
			print(c(null.one,null.two,date()));
		}; ## end null.two loop
		
	}; ## end null.one loop
	
	if(as.distance.matrix){ ## return as distance matrix if so desired
		results<-as.dist(results)
	}	
	
return(results)

} ## end function

RC.results = raup_crick_abundance(spXsite=otu)



write.csv(data.frame(as.matrix(RC.results)),"RC_results.csv",quote=F)



