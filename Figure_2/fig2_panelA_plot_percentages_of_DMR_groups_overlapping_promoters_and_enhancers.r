# Jennifer Flynn jaflynn@wustl.edu
# June 17, 2020

# Load the packages
library(bedr)

# Define functions
determine_dmr_overlap <- function(dmr_list, overlap_list)
{
  # Determine valid regions
  is.overlap_list.valid  <- is.valid.region(overlap_list);
  overlap_list <- overlap_list[is.overlap_list.valid]
  
  # Intersect the files
  dmr_overlap=bedr(input = list(a = bedr.sort.region(dmr_list), b = bedr.sort.region(overlap_list)), method="intersect", params="")
  dmr_overlap_df=data.frame(dmr_overlap)
  dmr_overlap_df$chr=apply(dmr_overlap_df, 1, function(x) strsplit(x, split=":")[[1]][1])
  dmr_overlap_df$region=apply(dmr_overlap_df, 1, function(x) strsplit(x[1], split=":")[[1]][2])
  dmr_overlap_df$start=as.numeric(apply(dmr_overlap_df, 1, function(x) strsplit(x[3], split="-")[[1]][1]))
  dmr_overlap_df$end=as.numeric(apply(dmr_overlap_df, 1, function(x) strsplit(x[3], split="-")[[1]][2]))
  dmr_overlap_nucleotides=sum(as.numeric(dmr_overlap_df$end)-as.numeric(dmr_overlap_df$start))
  return(dmr_overlap_nucleotides)
}

determine_dmrs_not_in_promoters_or_enhancers <- function(dmr_list, promoters_list, enhancers_list, promoters_and_enhancers_list)
{
  # Determine valid regions
  is.promoters_list.valid  <- is.valid.region(promoters_list);
  promoters_list <- promoters_list[is.promoters_list.valid]
  is.enhancers_list.valid  <- is.valid.region(enhancers_list);
  enhancers_list <- enhancers_list[is.enhancers_list.valid]
  is.promoters_and_enhancers_list.valid  <- is.valid.region(promoters_and_enhancers_list);
  promoters_and_enhancers_list <- promoters_and_enhancers_list[is.promoters_and_enhancers_list.valid]
  
  # Determine DMRs not within promoters or enhancers
  dmrs_not_in_promoters <- bedr.subtract.region(bedr.sort.region(dmr_list), bedr.sort.region(promoters_list), remove.whole.feature = FALSE)
  dmrs_not_in_promoters_or_enhancers <- bedr.subtract.region(bedr.sort.region(dmrs_not_in_promoters), bedr.sort.region(enhancers_list), remove.whole.feature = FALSE)
  dmrs_not_in_promoters_or_enhancers <- bedr.subtract.region(bedr.sort.region(dmrs_not_in_promoters_or_enhancers), bedr.sort.region(promoters_and_enhancers_list), remove.whole.feature = FALSE)
  dmrs_not_in_promoters_or_enhancers_df=data.frame(dmrs_not_in_promoters_or_enhancers)
  dmrs_not_in_promoters_or_enhancers_df$chr=apply(dmrs_not_in_promoters_or_enhancers_df, 1, function(x) strsplit(x, split=":")[[1]][1])
  dmrs_not_in_promoters_or_enhancers_df$region=apply(dmrs_not_in_promoters_or_enhancers_df, 1, function(x) strsplit(x[1], split=":")[[1]][2])
  dmrs_not_in_promoters_or_enhancers_df$start=as.numeric(apply(dmrs_not_in_promoters_or_enhancers_df, 1, function(x) strsplit(x[3], split="-")[[1]][1]))
  dmrs_not_in_promoters_or_enhancers_df$end=as.numeric(apply(dmrs_not_in_promoters_or_enhancers_df, 1, function(x) strsplit(x[3], split="-")[[1]][2]))
  dmrs_not_in_promoters_or_enhancers_nucleotides=sum(as.numeric(dmrs_not_in_promoters_or_enhancers_df$end)-as.numeric(dmrs_not_in_promoters_or_enhancers_df$start))
  return(dmrs_not_in_promoters_or_enhancers_nucleotides)
}

# Read in the data
# Read in the MERGED DMR groups:
EAC_hyperDMRs=read.table("EAC_merged_hyperDMRs.bed")
EAC_hyperDMRs_list=paste0(EAC_hyperDMRs$V1, ":", EAC_hyperDMRs$V2, "-", EAC_hyperDMRs$V3)
EAC_hypoDMRs=read.table("EAC_merged_hypoDMRs.bed")
EAC_hypoDMRs_list=paste0(EAC_hypoDMRs$V1, ":", EAC_hypoDMRs$V2, "-", EAC_hypoDMRs$V3)
GBM_hyperDMRs=read.table("GBM_merged_hyperDMRs.bed")
GBM_hyperDMRs_list=paste0(GBM_hyperDMRs$V1, ":", GBM_hyperDMRs$V2, "-", GBM_hyperDMRs$V3)
GBM_hypoDMRs=read.table("GBM_merged_hypoDMRs.bed")
GBM_hypoDMRs_list=paste0(GBM_hypoDMRs$V1, ":", GBM_hypoDMRs$V2, "-", GBM_hypoDMRs$V3)
  
# Read in regions annotated as active enhancers and not promoters (output from script: "fig2_panelA_determine_epi_defined_promoters_and_enhancers.sh")
active_enhancers_not_promoters=read.table("active_enhancers_not_overlapping_promoters.bed")
active_enhancers_not_promoters_list=paste0(active_enhancers_not_promoters$V1, ":", active_enhancers_not_promoters$V2, "-", active_enhancers_not_promoters$V3)
# Read in regions annotated as promoters and not active enhancers (output from script: "fig2_panelA_determine_epi_defined_promoters_and_enhancers.sh")
promoters_not_active_enhancers=read.table("promoters_not_overlapping_active_enhancers.bed")
promoters_not_active_enhancers_list=paste0(promoters_not_active_enhancers$V1, ":", promoters_not_active_enhancers$V2, "-", promoters_not_active_enhancers$V3)
# Read in regions annotated as promoters and active enhancers (output from script: "fig2_panelA_determine_epi_defined_promoters_and_enhancers.sh")
promoters_and_active_enhancers=read.table("promoters_and_active_enhancers.bed")
promoters_and_active_enhancers_list=paste0(promoters_and_active_enhancers$V1, ":", promoters_and_active_enhancers$V2, "-", promoters_and_active_enhancers$V3)

# Determine which DMRs overlap promoters only
data_for_plotting=data.frame(nucleotides=c(determine_dmr_overlap(EAC_hyperDMRs_list, promoters_not_active_enhancers_list),
                                           determine_dmr_overlap(EAC_hypoDMRs_list, promoters_not_active_enhancers_list),
                                           determine_dmr_overlap(GBM_hyperDMRs_list, promoters_not_active_enhancers_list),
                                           determine_dmr_overlap(GBM_hypoDMRs_list, promoters_not_active_enhancers_list),
                                           determine_dmr_overlap(EAC_hyperDMRs_list, active_enhancers_not_promoters_list),
                                           determine_dmr_overlap(EAC_hypoDMRs_list, active_enhancers_not_promoters_list),
                                           determine_dmr_overlap(GBM_hyperDMRs_list, active_enhancers_not_promoters_list),
                                           determine_dmr_overlap(GBM_hypoDMRs_list, active_enhancers_not_promoters_list),
                                           determine_dmr_overlap(EAC_hyperDMRs_list, promoters_and_active_enhancers_list),
                                           determine_dmr_overlap(EAC_hypoDMRs_list, promoters_and_active_enhancers_list),
                                           determine_dmr_overlap(GBM_hyperDMRs_list, promoters_and_active_enhancers_list),
                                           determine_dmr_overlap(GBM_hypoDMRs_list, promoters_and_active_enhancers_list),
                                           determine_dmrs_not_in_promoters_or_enhancers(EAC_hyperDMRs_list, promoters_not_active_enhancers_list, active_enhancers_not_promoters_list, promoters_and_active_enhancers_list),
                                           determine_dmrs_not_in_promoters_or_enhancers(EAC_hypoDMRs_list, promoters_not_active_enhancers_list, active_enhancers_not_promoters_list, promoters_and_active_enhancers_list),
                                           determine_dmrs_not_in_promoters_or_enhancers(GBM_hyperDMRs_list, promoters_not_active_enhancers_list, active_enhancers_not_promoters_list, promoters_and_active_enhancers_list),
                                           determine_dmrs_not_in_promoters_or_enhancers(GBM_hypoDMRs_list, promoters_not_active_enhancers_list, active_enhancers_not_promoters_list, promoters_and_active_enhancers_list)),
                                 dmr_group=rep(c("EAC_hyperDMRs", "EAC_hypoDMRs", "GBM_hyperDMRs", "GBM_hypoDMRs"),4),
                                 overlap_group=c(rep("promoters_only", 4), rep("enhancers_only", 4), rep("promoters_and_enhancers", 4), rep("neither_promoters_nor_enhancers", 4)))

# Plot the data
# EAC HyperDMRs
total_nucleotides=sum(c(data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hyperDMRs" & data_for_plotting$overlap_group=="promoters_only"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hyperDMRs" & data_for_plotting$overlap_group=="enhancers_only"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hyperDMRs" & data_for_plotting$overlap_group=="promoters_and_enhancers"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hyperDMRs" & data_for_plotting$overlap_group=="neither_promoters_nor_enhancers"]))
slices <- c(data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hyperDMRs" & data_for_plotting$overlap_group=="promoters_only"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hyperDMRs" & data_for_plotting$overlap_group=="enhancers_only"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hyperDMRs" & data_for_plotting$overlap_group=="promoters_and_enhancers"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hyperDMRs" & data_for_plotting$overlap_group=="neither_promoters_nor_enhancers"])
pct <- round((slices/sum(EAC_hyperDMRs$V3-EAC_hyperDMRs$V2))*100, digits = 2)
lbls <- c("Promoters Only", "Active Enhancers Only", "Promoters and Active Enhancers", "Other")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices, labels = lbls, main="EAC HyperDMRs", cex=.75)
dev.off()

# EAC HypoDMRs
total_nucleotides=sum(c(data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hypoDMRs" & data_for_plotting$overlap_group=="promoters_only"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hypoDMRs" & data_for_plotting$overlap_group=="enhancers_only"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hypoDMRs" & data_for_plotting$overlap_group=="promoters_and_enhancers"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hypoDMRs" & data_for_plotting$overlap_group=="neither_promoters_nor_enhancers"]))
slices <- c(data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hypoDMRs" & data_for_plotting$overlap_group=="promoters_only"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hypoDMRs" & data_for_plotting$overlap_group=="enhancers_only"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hypoDMRs" & data_for_plotting$overlap_group=="promoters_and_enhancers"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="EAC_hypoDMRs" & data_for_plotting$overlap_group=="neither_promoters_nor_enhancers"])
pct <- round((slices/sum(EAC_hypoDMRs$V3-EAC_hypoDMRs$V2))*100, digits = 2)
lbls <- c("Promoters Only", "Active Enhancers Only", "Promoters and Active Enhancers", "Other")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices, labels = lbls, main="EAC hypoDMRs", cex=.75)
dev.off()

# GBM HyperDMRs
total_nucleotides=sum(c(data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hyperDMRs" & data_for_plotting$overlap_group=="promoters_only"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hyperDMRs" & data_for_plotting$overlap_group=="enhancers_only"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hyperDMRs" & data_for_plotting$overlap_group=="promoters_and_enhancers"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hyperDMRs" & data_for_plotting$overlap_group=="neither_promoters_nor_enhancers"]))
slices <- c(data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hyperDMRs" & data_for_plotting$overlap_group=="promoters_only"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hyperDMRs" & data_for_plotting$overlap_group=="enhancers_only"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hyperDMRs" & data_for_plotting$overlap_group=="promoters_and_enhancers"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hyperDMRs" & data_for_plotting$overlap_group=="neither_promoters_nor_enhancers"])
pct <- round((slices/sum(GBM_hyperDMRs$V3-GBM_hyperDMRs$V2))*100, digits = 2)
lbls <- c("Promoters Only", "Active Enhancers Only", "Promoters and Active Enhancers", "Other")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices, labels = lbls, main="GBM HyperDMRs", cex=.75)
dev.off()

# GBM HypoDMRs
total_nucleotides=sum(c(data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hypoDMRs" & data_for_plotting$overlap_group=="promoters_only"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hypoDMRs" & data_for_plotting$overlap_group=="enhancers_only"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hypoDMRs" & data_for_plotting$overlap_group=="promoters_and_enhancers"], data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hypoDMRs" & data_for_plotting$overlap_group=="neither_promoters_nor_enhancers"]))
slices <- c(data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hypoDMRs" & data_for_plotting$overlap_group=="promoters_only"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hypoDMRs" & data_for_plotting$overlap_group=="enhancers_only"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hypoDMRs" & data_for_plotting$overlap_group=="promoters_and_enhancers"],
            data_for_plotting$nucleotides[data_for_plotting$dmr_group=="GBM_hypoDMRs" & data_for_plotting$overlap_group=="neither_promoters_nor_enhancers"])
pct <- round((slices/sum(GBM_hypoDMRs$V3-GBM_hypoDMRs$V2))*100, digits = 2)
lbls <- c("Promoters Only", "Active Enhancers Only", "Promoters and Active Enhancers", "Other")
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices, labels = lbls, main="GBM hypoDMRs", cex=.75)
dev.off()
