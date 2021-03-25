# Jennifer Flynn jaflynn@wustl.edu
# June 27, 2016
# This script determines the percentage of DMR nucleotides that overlap certain genomic annotations

# Load libraries
library(bedr)
library(stringr)
library(ggplot2)

# Define functions

# input: list of DMR locations, locations of genomic regions of interest, and whether or not the background to be considered is EAC or GBM.
# output: x: list of length 2 - first element is the percentage of DMR nucleotides that overlap the genomic regions, second element is the enrichment of the DMR group within the genomic regions
calculate_percentage_overlap_and_enrichment <- function(DMR_group, genomic_regions, background)
{
  # Read in the DMRs and label columns
  colnames(DMR_group)=c("chr", "start", "end")
  DMR_group$chr=as.character(DMR_group$chr)
  DMR_group$start=as.numeric(DMR_group$start)
  DMR_group$end=as.numeric(DMR_group$end)
  
  # Determine number of nucleotides in DMR group
  DMR_nucleotides=sum(DMR_group$end-DMR_group$start)
  
  # Read in the genomic regions of interest and label columns
  colnames(genomic_regions)=c("chr", "start", "end")
  genomic_regions$chr=as.character(genomic_regions$chr)
  genomic_regions=genomic_regions[genomic_regions$chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"), ] 
  genomic_regions$start=as.numeric(genomic_regions$start)
  genomic_regions$end=as.numeric(genomic_regions$end)
  
  # Intersect regions of interest and DMR file
  DMR_group_list=paste0(DMR_group$chr, ":", DMR_group$start, "-", DMR_group$end)
  genomic_regions_list=paste0(genomic_regions$chr, ":", genomic_regions$start, "-", genomic_regions$end)
  overlap_regions=bedr(input = list(a = DMR_group_list, b = genomic_regions_list), method="intersect", params="")
  
  # Calculate number of nucleotides in overlap
  overlap_regions_df=data.frame(overlap_regions)
  overlap_regions_df$chr=apply(overlap_regions_df, 1, function(x) strsplit(x, split=":")[[1]][1])
  overlap_regions_df$region=apply(overlap_regions_df, 1, function(x) strsplit(x[1], split=":")[[1]][2])
  overlap_regions_df$start=as.numeric(apply(overlap_regions_df, 1, function(x) strsplit(x[3], split="-")[[1]][1]))
  overlap_regions_df$end=as.numeric(apply(overlap_regions_df, 1, function(x) strsplit(x[3], split="-")[[1]][2]))
  overlap_nucleotides=sum(overlap_regions_df$end-overlap_regions_df$start)
  
  # Calcualte percentage of DMR nucleotides in overlap
  DMR_percentage=overlap_nucleotides/DMR_nucleotides
  
  # Determine number of nucleotides in genomic regions for background (EAC or GBM depending on DMR group)
  if (background=="EAC")
  {
    background_file=EAC_background_regions_list
    background_nucleotides=EAC_background_nucleotides
  } else {
    background_file=GBM_background_regions_list
    background_nucleotides=GBM_background_nucleotides
  }
  background_overlap_regions=bedr(input = list(a = background_file, b = genomic_regions_list), method="intersect", params="")
  background_overlap_regions_df=data.frame(background_overlap_regions)
  background_overlap_regions_df$chr=apply(background_overlap_regions_df, 1, function(x) strsplit(x, split=":")[[1]][1])
  background_overlap_regions_df$region=apply(background_overlap_regions_df, 1, function(x) strsplit(x[1], split=":")[[1]][2])
  background_overlap_regions_df$start=as.numeric(apply(background_overlap_regions_df, 1, function(x) strsplit(x[3], split="-")[[1]][1]))
  background_overlap_regions_df$end=as.numeric(apply(background_overlap_regions_df, 1, function(x) strsplit(x[3], split="-")[[1]][2]))
  background_overlap_nucleotides=sum(background_overlap_regions_df$end-background_overlap_regions_df$start)
  
  # Calcualte the percentage of background within genomic region
  background_percentage=background_overlap_nucleotides/background_nucleotides
  
  #Determine enrichment of DMR group in genomic region (%DMR nucleotides overlapping region / %background nucleotides overlapping region)
  enrichment=log2(DMR_percentage/background_percentage)
  return(c(DMR_percentage, enrichment))
}

# Load in DMR data
EAC_unique_hyperDMRs=read.table("EAC_unique_hyperDMRs.bed")
GBM_unique_hyperDMRs=read.table("GBM_unique_hyperDMRs.bed")
EAC_unique_hypoDMRs=read.table("EAC_unique_hypoDMRs.bed")
GBM_unique_hypoDMRs=read.table("GBM_unique_hypoDMRs.bed")
shared_hyperDMRs=read.table("shared_EAC_GBM_hyperDMRs.bed")
shared_hypoDMRs=read.table("shared_EAC_GBM_hypoDMRs.bed")

# Genomic category files
cpg_islands_unmasked=read.table("cpgIslandExtUnmasked.bed")
promoters_2500bp=read.table("refGene_Promoters2500bp_sorted_merged.bed")
core_promoters_1000bp=read.table("refGene_corePromoters1000bp_sorted_merged.bed")
gene_bodies=read.table("refGene_transcription_sorted_merged.bed")
fantom_enhancers=read.table("human_permissive_enhancers_phase_1_and_2.bed")
vista_enhancers=read.table("vista_human_enhancers_sorted_merged.bed")
super_enhancers=read.table("All.Super.Enhancer.hg19_sorted_merged.bed")

# Background regions for EAC and GBM
EAC_background_regions=read.table("fivehundred_bp_genomic_bins_excluding_those_that_overlap_to_blacklisted_cpgs_and_those_without_a_cpg_and_those_without_MeDIP_andor_MRE_for_EAC.bed") #Please note this file must be downloaded from https://drive.google.com/drive/folders/1ft61CuNAG7_CRXoRwkBsB1vGUgxxbLX5?usp=sharing and unzipped prior to using here.
colnames(EAC_background_regions)=c("chr", "start", "end")
EAC_background_regions$start=as.numeric(EAC_background_regions$start)
EAC_background_regions$end=as.numeric(EAC_background_regions$end)
EAC_background_nucleotides=sum(EAC_background_regions$end-EAC_background_regions$start)
EAC_background_regions_list=paste0(EAC_background_regions$chr, ":", EAC_background_regions$start, "-", EAC_background_regions$end)
GBM_background_regions=read.table("fivehundred_bp_genomic_bins_excluding_those_that_overlap_to_blacklisted_cpgs_and_those_without_a_cpg_and_those_without_MeDIP_andor_MRE_for_GBM.bed") #Please note this file must be downloaded from https://drive.google.com/drive/folders/1ft61CuNAG7_CRXoRwkBsB1vGUgxxbLX5?usp=sharing and unzipped prior to using here.
colnames(GBM_background_regions)=c("chr", "start", "end")
GBM_background_regions$start=as.numeric(GBM_background_regions$start)
GBM_background_regions$end=as.numeric(GBM_background_regions$end)
GBM_background_nucleotides=sum(GBM_background_regions$end-GBM_background_regions$start)
GBM_background_regions_list=paste0(GBM_background_regions$chr, ":", GBM_background_regions$start, "-", GBM_background_regions$end)

# Calculate the percentage overlap and enrichment for each DMR group to each genomic category
DMR_groups=c("EAC_unique_hyperDMRs", "EAC_unique_hypoDMRs", "GBM_unique_hyperDMRs", "GBM_unique_hypoDMRs", "shared_hyperDMRs", "shared_hypoDMRs")
genomic_categories=c("cpg_islands_unmasked", "promoters_2500bp", "core_promoters_1000bp", "gene_bodies", "fantom_enhancers", "vista_enhancers", "super_enhancers")
for (group in DMR_groups)
{
  for (category in genomic_categories)
  {
    if (group=="EAC_unique_hyperDMRs" | group=="EAC_unique_hypoDMRs")
    {
      background_file="EAC"
    } else
    {
      background_file="GBM"
    }
    percentage_overlap_and_enrichment=calculate_percentage_overlap_and_enrichment(eval(parse(text = eval(parse(text = "group")))), eval(parse(text = eval(parse(text = "category")))),  background_file)
    
    # Gather data for plotting
    if (exists("plotting_data"))
    {
      data=data.frame(DMR_group=group, Genomic_Annotation=category, Percentage=percentage_overlap_and_enrichment[1], Log2Enrichment=percentage_overlap_and_enrichment[2])
      rbind(plotting_data, data)
    } else
    {
      plotting_data=data.frame(DMR_group=group, Genomic_Annotation=category, Percentage=percentage_overlap_and_enrichment[1], Log2Enrichment=percentage_overlap_and_enrichment[2])
    }
  }
}

# Plot the data
ggplot(plotting_data, aes(factor(Genomic_Annotation), Percentage, fill = DMR_group)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  labs(title="Genomic Annotation Distributions", x = "", y = "Fraction of DMR Nucleotides within Genomic Category") +
  ylim(0,1) +
  theme(axis.text.y = element_text(size=10), axis.title = element_text(size=7)) + 
  theme(axis.text.x = element_text(size=5, angle=90))

ggplot(plotting_data, aes(factor(Genomic_Annotation), Log2Enrichment, fill = DMR_group)) +
  geom_bar(stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  labs(y = "Enrichment of Genomic Category within DMR Group (log2)", x = "Genomic Annotation Categories") +
  theme(axis.text.y = element_text(size=10), axis.title = element_text(size=7)) +
  theme(axis.text.x = element_text(size=5, angle=90))
