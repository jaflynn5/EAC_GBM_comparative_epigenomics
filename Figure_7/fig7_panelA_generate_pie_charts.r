# Jennifer Flynn jaflynn@wustl.edu
# June 20, 2020
# This script reads DMR sets and TE locations and determines their overlap.

# Load in the appropriate libraries:
library(bedr)

# Define functions
determine_overlapping_nucleotides <- function(group1, group2)
{
  # Make sure the column names are correct
  colnames(group1)=c("chr", "start", "end")
  group1$chr=as.character(group1$chr)
  group1$start=as.numeric(group1$start)
  group1$end=as.numeric(group1$end)
  colnames(group2)=c("chr", "start", "end")
  group2$chr=as.character(group2$chr)
  group2$start=as.numeric(group2$start)
  group2$end=as.numeric(group2$end)
  
  # Overlap the regions
  overlapping_regions=bedr(input = list(a = bedr.sort.region(group1[is.valid.region(group1),]), b = bedr.sort.region(group2[is.valid.region(group2),])), method="intersect", params="-loj")
  overlapping_regions=overlapping_regions[overlapping_regions$chr.b!=".",]
  overlapping_regions$start.b=as.numeric(overlapping_regions$start.b)
  overlapping_regions$end.b=as.numeric(overlapping_regions$end.b)
  overlapping_regions=data.frame(chr=overlapping_regions$chr,
                                 start=apply(overlapping_regions, 1, function(x) if (x[5]<x[2]) x[2] else x[5]),
                                 end=apply(overlapping_regions, 1, function(x) if (x[3]>x[6]) x[6] else x[3]))
  overlapping_regions$chr=as.character(overlapping_regions$chr)
  overlapping_regions$start=as.numeric(as.character(overlapping_regions$start))
  overlapping_regions$end=as.numeric(as.character(overlapping_regions$end))
  
  # Sort the output
  overlapping_regions_sorted=bedr.sort.region(overlapping_regions)
  
  # Merge the output
  overlapping_regions_sorted_merged=bedr.merge.region(overlapping_regions_sorted)
  
  # Return the result
  return(overlapping_regions_sorted_merged)
}

# Load in the data
EAC_hyperDMRs=read.table("EAC_hyperDMRs.bed")
colnames(EAC_hyperDMRs)=c("chr", "start", "end")
EAC_hyperDMRs$chr=as.character(EAC_hyperDMRs$chr)
EAC_hypoDMRs=read.table("EAC_hypoDMRs.bed")
colnames(EAC_hypoDMRs)=c("chr", "start", "end")
EAC_hypoDMRs$chr=as.character(EAC_hypoDMRs$chr)
GBM_hyperDMRs=read.table("GBM_hyperDMRs.bed")
colnames(GBM_hyperDMRs)=c("chr", "start", "end")
GBM_hyperDMRs$chr=as.character(GBM_hyperDMRs$chr)
GBM_hypoDMRs=read.table("GBM_hypoDMRs.bed")
colnames(GBM_hypoDMRs)=c("chr", "start", "end")
GBM_hypoDMRs$chr=as.character(GBM_hypoDMRs$chr)
all_TEs=read.table("hg19_rmsk.bed") #Download from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
all_TEs=all_TEs[,1:3]
colnames(all_TEs)=c("chr", "start", "end")
all_TEs$chr=as.character(all_TEs$chr)
simple_repeats=read.table("simpleRepeat.bed") #Download from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
simple_repeats=simple_repeats[,1:3]
colnames(simple_repeats)=c("chr", "start", "end")
simple_repeats$chr=as.character(simple_repeats$chr)
genomic_promoters=read.table("refGene_Promoters2500bp_sorted_merged.bed")
colnames(genomic_promoters)=c("chr", "start", "end")
genomic_promoters$chr=as.character(genomic_promoters$chr)
epigenomic_and_or_genomic_promoter_sites=read.table("union_of_genomic_and_epigenomic_promoters.bed")# Output file from fig2_panelA_determine_epi_defined_promoters_and_enhancers.sh
colnames(epigenomic_and_or_genomic_promoter_sites)=c("chr", "start", "end")
epigenomic_and_or_genomic_promoter_sites$chr=as.character(epigenomic_and_or_genomic_promoter_sites$chr)
epigenomic_enhancer_sites_sorted_merged=read.table("active_enhancer_sites_sorted_merged.bed")# Output file from fig2_panelA_determine_epi_defined_promoters_and_enhancers.sh
colnames(epigenomic_enhancer_sites_sorted_merged)=c("chr", "start", "end")
epigenomic_enhancer_sites_sorted_merged$chr=as.character(epigenomic_enhancer_sites_sorted_merged$chr)

# Remove simple repeats from the analysis
tes_excluding_simple_repeats=bedr(input = list(a = bedr.sort.region(all_TEs[is.valid.region(all_TEs),]), b = bedr.sort.region(simple_repeats[is.valid.region(simple_repeats),])), method="subtract", params="")


#### EAC HyperDMRs ####
# Determine which DMR nucleotides overlap TEs (excluding simple repeats)
eac_hyperDMRs_overlapping_tes_sorted_merged=determine_overlapping_nucleotides(EAC_hyperDMRs, tes_excluding_simple_repeats)

# Determine which DMR-TE overlap nucleotides overlap genic promoter definitions
eac_hyperDMRs_overlapping_tes_and_genic_promoters_sorted_merged=determine_overlapping_nucleotides(eac_hyperDMRs_overlapping_tes_sorted_merged, genomic_promoters)

# Determine which DMR-TE overlap nucleotides do not overlap genic promoter definitions
eac_hyperDMRs_overlapping_tes_not_genic_promoters = bedr(input = list(a = bedr.sort.region(eac_hyperDMRs_overlapping_tes_sorted_merged), b = bedr.sort.region(genomic_promoters[is.valid.region(genomic_promoters),])), method="subtract", params="")

# Determine which DMR-TE overlap nucleotides that do not overlap to genic promoters overlap to epigenetic promoters
eac_hyperDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged=determine_overlapping_nucleotides(eac_hyperDMRs_overlapping_tes_not_genic_promoters, epigenomic_and_or_genomic_promoter_sites)
  
# Determine which DMR-TE overlap nucleotides do not overlap to genic promoters or epigenetic promoters
eac_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters = bedr(input = list(a = bedr.sort.region(eac_hyperDMRs_overlapping_tes_not_genic_promoters), b = bedr.sort.region(epigenomic_and_or_genomic_promoter_sites[is.valid.region(epigenomic_and_or_genomic_promoter_sites),])), method="subtract", params="")

# Determine which DMR-TE overlap nucleotides that do not overlap genomic or epigenomic promoters overlap active enhancers
eac_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers=determine_overlapping_nucleotides(eac_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters, epigenomic_enhancer_sites_sorted_merged)

# Determine which DMR-TE overlap nucleotides do not overlap genomic or epigenomic promoters or active enhancers
eac_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers = bedr(input = list(a = bedr.sort.region(eac_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters), b = bedr.sort.region(epigenomic_enhancer_sites_sorted_merged[is.valid.region(epigenomic_enhancer_sites_sorted_merged),])), method="subtract", params="")

# Create plot for EAC hyperDMRs
slices <- c(sum(eac_hyperDMRs_overlapping_tes_and_genic_promoters_sorted_merged$end-eac_hyperDMRs_overlapping_tes_and_genic_promoters_sorted_merged$start)/sum(eac_hyperDMRs_overlapping_tes_sorted_merged$end-eac_hyperDMRs_overlapping_tes_sorted_merged$start),
            sum(eac_hyperDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged$end-eac_hyperDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged$start)/sum(eac_hyperDMRs_overlapping_tes_sorted_merged$end-eac_hyperDMRs_overlapping_tes_sorted_merged$start),
            sum(eac_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers$end-eac_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers$start)/sum(eac_hyperDMRs_overlapping_tes_sorted_merged$end-eac_hyperDMRs_overlapping_tes_sorted_merged$start),
            sum(eac_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers$end-eac_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers$start)/sum(eac_hyperDMRs_overlapping_tes_sorted_merged$end-eac_hyperDMRs_overlapping_tes_sorted_merged$start))
pct <- round(slices*100, digits = 2)
lbls <- c("Genomic Promoters", "Epigenomic Promoters", "Active Enhancers", "Other")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices, labels = lbls, main="EAC HyperDMRs", cex=.75)


#### EAC HypoDMRs ####
# Determine which DMR nucleotides overlap TEs (excluding simple repeats)
eac_hypoDMRs_overlapping_tes_sorted_merged=determine_overlapping_nucleotides(EAC_hypoDMRs, tes_excluding_simple_repeats)

# Determine which DMR-TE overlap nucleotides overlap genic promoter definitions
eac_hypoDMRs_overlapping_tes_and_genic_promoters_sorted_merged=determine_overlapping_nucleotides(eac_hypoDMRs_overlapping_tes_sorted_merged, genomic_promoters)

# Determine which DMR-TE overlap nucleotides do not overlap genic promoter definitions
eac_hypoDMRs_overlapping_tes_not_genic_promoters = bedr(input = list(a = bedr.sort.region(eac_hypoDMRs_overlapping_tes_sorted_merged), b = bedr.sort.region(genomic_promoters[is.valid.region(genomic_promoters),])), method="subtract", params="")

# Determine which DMR-TE overlap nucleotides that do not overlap to genic promoters overlap to epigenetic promoters
eac_hypoDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged=determine_overlapping_nucleotides(eac_hypoDMRs_overlapping_tes_not_genic_promoters, epigenomic_and_or_genomic_promoter_sites)

# Determine which DMR-TE overlap nucleotides do not overlap to genic promoters or epigenetic promoters
eac_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters = bedr(input = list(a = bedr.sort.region(eac_hypoDMRs_overlapping_tes_not_genic_promoters), b = bedr.sort.region(epigenomic_and_or_genomic_promoter_sites[is.valid.region(epigenomic_and_or_genomic_promoter_sites),])), method="subtract", params="")

# Determine which DMR-TE overlap nucleotides that do not overlap genomic or epigenomic promoters overlap active enhancers
eac_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers=determine_overlapping_nucleotides(eac_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters, epigenomic_enhancer_sites_sorted_merged)

# Determine which DMR-TE overlap nucleotides do not overlap genomic or epigenomic promoters or active enhancers
eac_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers = bedr(input = list(a = bedr.sort.region(eac_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters), b = bedr.sort.region(epigenomic_enhancer_sites_sorted_merged[is.valid.region(epigenomic_enhancer_sites_sorted_merged),])), method="subtract", params="")

# Create plot for EAC hypoDMRs
slices <- c(sum(eac_hypoDMRs_overlapping_tes_and_genic_promoters_sorted_merged$end-eac_hypoDMRs_overlapping_tes_and_genic_promoters_sorted_merged$start)/sum(eac_hypoDMRs_overlapping_tes_sorted_merged$end-eac_hypoDMRs_overlapping_tes_sorted_merged$start),
            sum(eac_hypoDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged$end-eac_hypoDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged$start)/sum(eac_hypoDMRs_overlapping_tes_sorted_merged$end-eac_hypoDMRs_overlapping_tes_sorted_merged$start),
            sum(eac_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers$end-eac_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers$start)/sum(eac_hypoDMRs_overlapping_tes_sorted_merged$end-eac_hypoDMRs_overlapping_tes_sorted_merged$start),
            sum(eac_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers$end-eac_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers$start)/sum(eac_hypoDMRs_overlapping_tes_sorted_merged$end-eac_hypoDMRs_overlapping_tes_sorted_merged$start))
pct <- round(slices*100, digits = 2)
lbls <- c("Genomic Promoters", "Epigenomic Promoters", "Active Enhancers", "Other")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices, labels = lbls, main="EAC HypoDMRs", cex=.75)


#### GBM HyperDMRs ####
# Determine which DMR nucleotides overlap TEs (excluding simple repeats)
gbm_hyperDMRs_overlapping_tes_sorted_merged=determine_overlapping_nucleotides(GBM_hyperDMRs, tes_excluding_simple_repeats)

# Determine which DMR-TE overlap nucleotides overlap genic promoter definitions
gbm_hyperDMRs_overlapping_tes_and_genic_promoters_sorted_merged=determine_overlapping_nucleotides(gbm_hyperDMRs_overlapping_tes_sorted_merged, genomic_promoters)

# Determine which DMR-TE overlap nucleotides do not overlap genic promoter definitions
gbm_hyperDMRs_overlapping_tes_not_genic_promoters = bedr(input = list(a = bedr.sort.region(gbm_hyperDMRs_overlapping_tes_sorted_merged), b = bedr.sort.region(genomic_promoters[is.valid.region(genomic_promoters),])), method="subtract", params="")

# Determine which DMR-TE overlap nucleotides that do not overlap to genic promoters overlap to epigenetic promoters
gbm_hyperDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged=determine_overlapping_nucleotides(gbm_hyperDMRs_overlapping_tes_not_genic_promoters, epigenomic_and_or_genomic_promoter_sites)

# Determine which DMR-TE overlap nucleotides do not overlap to genic promoters or epigenetic promoters
gbm_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters = bedr(input = list(a = bedr.sort.region(gbm_hyperDMRs_overlapping_tes_not_genic_promoters), b = bedr.sort.region(epigenomic_and_or_genomic_promoter_sites[is.valid.region(epigenomic_and_or_genomic_promoter_sites),])), method="subtract", params="")

# Determine which DMR-TE overlap nucleotides that do not overlap genomic or epigenomic promoters overlap active enhancers
gbm_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers=determine_overlapping_nucleotides(gbm_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters, epigenomic_enhancer_sites_sorted_merged)

# Determine which DMR-TE overlap nucleotides do not overlap genomic or epigenomic promoters or active enhancers
gbm_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers = bedr(input = list(a = bedr.sort.region(gbm_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters), b = bedr.sort.region(epigenomic_enhancer_sites_sorted_merged[is.valid.region(epigenomic_enhancer_sites_sorted_merged),])), method="subtract", params="")

# Create plot for GBM hyperDMRs
slices <- c(sum(gbm_hyperDMRs_overlapping_tes_and_genic_promoters_sorted_merged$end-gbm_hyperDMRs_overlapping_tes_and_genic_promoters_sorted_merged$start)/sum(gbm_hyperDMRs_overlapping_tes_sorted_merged$end-gbm_hyperDMRs_overlapping_tes_sorted_merged$start),
            sum(gbm_hyperDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged$end-gbm_hyperDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged$start)/sum(gbm_hyperDMRs_overlapping_tes_sorted_merged$end-gbm_hyperDMRs_overlapping_tes_sorted_merged$start),
            sum(gbm_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers$end-gbm_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers$start)/sum(gbm_hyperDMRs_overlapping_tes_sorted_merged$end-gbm_hyperDMRs_overlapping_tes_sorted_merged$start),
            sum(gbm_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers$end-gbm_hyperDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers$start)/sum(gbm_hyperDMRs_overlapping_tes_sorted_merged$end-gbm_hyperDMRs_overlapping_tes_sorted_merged$start))
pct <- round(slices*100, digits = 2)
lbls <- c("Genomic Promoters", "Epigenomic Promoters", "Active Enhancers", "Other")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices, labels = lbls, main="GBM HyperDMRs", cex=.75)


#### GBM HypoDMRs ####
# Determine which DMR nucleotides overlap TEs (excluding simple repeats)
gbm_hypoDMRs_overlapping_tes_sorted_merged=determine_overlapping_nucleotides(GBM_hypoDMRs, tes_excluding_simple_repeats)

# Determine which DMR-TE overlap nucleotides overlap genic promoter definitions
gbm_hypoDMRs_overlapping_tes_and_genic_promoters_sorted_merged=determine_overlapping_nucleotides(gbm_hypoDMRs_overlapping_tes_sorted_merged, genomic_promoters)

# Determine which DMR-TE overlap nucleotides do not overlap genic promoter definitions
gbm_hypoDMRs_overlapping_tes_not_genic_promoters = bedr(input = list(a = bedr.sort.region(gbm_hypoDMRs_overlapping_tes_sorted_merged), b = bedr.sort.region(genomic_promoters[is.valid.region(genomic_promoters),])), method="subtract", params="")

# Determine which DMR-TE overlap nucleotides that do not overlap to genic promoters overlap to epigenetic promoters
gbm_hypoDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged=determine_overlapping_nucleotides(gbm_hypoDMRs_overlapping_tes_not_genic_promoters, epigenomic_and_or_genomic_promoter_sites)

# Determine which DMR-TE overlap nucleotides do not overlap to genic promoters or epigenetic promoters
gbm_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters = bedr(input = list(a = bedr.sort.region(gbm_hypoDMRs_overlapping_tes_not_genic_promoters), b = bedr.sort.region(epigenomic_and_or_genomic_promoter_sites[is.valid.region(epigenomic_and_or_genomic_promoter_sites),])), method="subtract", params="")

# Determine which DMR-TE overlap nucleotides that do not overlap genomic or epigenomic promoters overlap active enhancers
gbm_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers=determine_overlapping_nucleotides(gbm_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters, epigenomic_enhancer_sites_sorted_merged)

# Determine which DMR-TE overlap nucleotides do not overlap genomic or epigenomic promoters or active enhancers
gbm_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers = bedr(input = list(a = bedr.sort.region(gbm_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters), b = bedr.sort.region(epigenomic_enhancer_sites_sorted_merged[is.valid.region(epigenomic_enhancer_sites_sorted_merged),])), method="subtract", params="")

# Create plot for GBM hypoDMRs
slices <- c(sum(gbm_hypoDMRs_overlapping_tes_and_genic_promoters_sorted_merged$end-gbm_hypoDMRs_overlapping_tes_and_genic_promoters_sorted_merged$start)/sum(gbm_hypoDMRs_overlapping_tes_sorted_merged$end-gbm_hypoDMRs_overlapping_tes_sorted_merged$start),
            sum(gbm_hypoDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged$end-gbm_hypoDMRs_overlapping_tes_not_genic_promoters_overlapping_epigenetic_promoters_sorted_merged$start)/sum(gbm_hypoDMRs_overlapping_tes_sorted_merged$end-gbm_hypoDMRs_overlapping_tes_sorted_merged$start),
            sum(gbm_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers$end-gbm_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_overlapping_enhancers$start)/sum(gbm_hypoDMRs_overlapping_tes_sorted_merged$end-gbm_hypoDMRs_overlapping_tes_sorted_merged$start),
            sum(gbm_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers$end-gbm_hypoDMRs_overlapping_tes_not_overlapping_genic_or_epigenetic_promoters_or_enhancers$start)/sum(gbm_hypoDMRs_overlapping_tes_sorted_merged$end-gbm_hypoDMRs_overlapping_tes_sorted_merged$start))
pct <- round(slices*100, digits = 2)
lbls <- c("Genomic Promoters", "Epigenomic Promoters", "Active Enhancers", "Other")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices, labels = lbls, main="GBM HypoDMRs", cex=.75)
