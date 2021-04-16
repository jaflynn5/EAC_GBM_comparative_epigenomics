# Jennifer Flynn jaflynn@wustl.edu
# 2/25/2016
# This script takes in a list of tumor suppressor genes, lists of DMRs, and a file containing the locations of genes, and creates venn diagrams depicting the overlap between TSGs with EAC and GBM DMRs within their core promoters.

# Load in the appropriate library(s):
library(bedr)
library(venneuler)
library(VennDiagram)

# Read in the data
full_tumor_suppressor_gene_file=read.table("Human_TSGs.txt", sep="\t", header=TRUE) # File downloaded from from "http://bioinfo.mc.vanderbilt.edu/TSGene/Human_TSGs.txt" on 02/25/2016
core_promoters_1000bp=read.table("refGene_corePromoters1000bp_withGeneNames_sorted.bed") # File adapted from downloaded file from "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/"

EAC_hyperDMRs=read.table("EAC_merged_hyperDMRs.bed")
colnames(EAC_hyperDMRs)=c("chr", "start", "end", "value")
EAC_hyperDMRs$chr=as.character(EAC_hyperDMRs$chr)
EAC_hypoDMRs=read.table("EAC_merged_hypoDMRs.bed")
colnames(EAC_hypoDMRs)=c("chr", "start", "end", "value")
EAC_hypoDMRs$chr=as.character(EAC_hypoDMRs$chr)
GBM_hyperDMRs=read.table("GBM_merged_hyperDMRs.bed")
colnames(GBM_hyperDMRs)=c("chr", "start", "end", "value")
GBM_hyperDMRs$chr=as.character(GBM_hyperDMRs$chr)
GBM_hypoDMRs=read.table("GBM_merged_hypoDMRs.bed")
colnames(GBM_hypoDMRs)=c("chr", "start", "end", "value")
GBM_hypoDMRs$chr=as.character(GBM_hypoDMRs$chr)

# Determine locations of TSGs
tsgs=unique(full_tumor_suppressor_gene_file$GeneSymbol)
tsgs_core_promoters=unique(core_promoters_1000bp[core_promoters_1000bp$V5%in%tsgs,])
colnames(tsgs_core_promoters)=c("chr", "start", "end", "strand", "gene_id")
tsgs_core_promoters$chr=as.character(tsgs_core_promoters$chr)
is.overlap_list.valid  <- is.valid.region(tsgs_core_promoters);
tsgs_core_promoters <- tsgs_core_promoters[is.overlap_list.valid,]

# Determine which TSGs overlap which DMR groups
tsgs_with_promoter_overlap_to_EAC_hyperDMRs=bedr(input = list(a = tsgs_core_promoters, b = EAC_hyperDMRs), method="intersect", params="-c")
tsgs_with_promoter_overlap_to_EAC_hyperDMRs=tsgs_with_promoter_overlap_to_EAC_hyperDMRs[tsgs_with_promoter_overlap_to_EAC_hyperDMRs$V6>0,]
tsgs_with_promoter_overlap_to_EAC_hyperDMRs=unique(tsgs_with_promoter_overlap_to_EAC_hyperDMRs$V5)
tsgs_with_promoter_overlap_to_EAC_hypoDMRs=bedr(input = list(a = tsgs_core_promoters, b = EAC_hypoDMRs), method="intersect", params="-c")
tsgs_with_promoter_overlap_to_EAC_hypoDMRs=tsgs_with_promoter_overlap_to_EAC_hypoDMRs[tsgs_with_promoter_overlap_to_EAC_hypoDMRs$V6>0,]
tsgs_with_promoter_overlap_to_EAC_hypoDMRs=unique(tsgs_with_promoter_overlap_to_EAC_hypoDMRs$V5)
tsgs_with_promoter_overlap_to_GBM_hyperDMRs=bedr(input = list(a = tsgs_core_promoters, b = GBM_hyperDMRs), method="intersect", params="-c")
tsgs_with_promoter_overlap_to_GBM_hyperDMRs=tsgs_with_promoter_overlap_to_GBM_hyperDMRs[tsgs_with_promoter_overlap_to_GBM_hyperDMRs$V6>0,]
tsgs_with_promoter_overlap_to_GBM_hyperDMRs=unique(tsgs_with_promoter_overlap_to_GBM_hyperDMRs$V5)
tsgs_with_promoter_overlap_to_GBM_hypoDMRs=bedr(input = list(a = tsgs_core_promoters, b = GBM_hypoDMRs), method="intersect", params="-c")
tsgs_with_promoter_overlap_to_GBM_hypoDMRs=tsgs_with_promoter_overlap_to_GBM_hypoDMRs[tsgs_with_promoter_overlap_to_GBM_hypoDMRs$V6>0,]
tsgs_with_promoter_overlap_to_GBM_hypoDMRs=unique(tsgs_with_promoter_overlap_to_GBM_hypoDMRs$V5)

# Draw Venn Diagrams
# HyperDMRs
EAC_hyperDMRs_unique=length(tsgs_with_promoter_overlap_to_EAC_hyperDMRs[!tsgs_with_promoter_overlap_to_EAC_hyperDMRs%in%tsgs_with_promoter_overlap_to_GBM_hyperDMRs])
GBM_hyperDMRs_unique=length(tsgs_with_promoter_overlap_to_GBM_hyperDMRs[!tsgs_with_promoter_overlap_to_GBM_hyperDMRs%in%tsgs_with_promoter_overlap_to_EAC_hyperDMRs])
hyperDMRs_shared=length(tsgs_with_promoter_overlap_to_EAC_hyperDMRs[tsgs_with_promoter_overlap_to_EAC_hyperDMRs %in% tsgs_with_promoter_overlap_to_GBM_hyperDMRs])
grid.newpage()
col1="#FF0000"
col2="#0000FF"
draw.pairwise.venn(area1 = EAC_hyperDMRs_unique+hyperDMRs_shared, area2 = GBM_hyperDMRs_unique+hyperDMRs_shared, cross.area = hyperDMRs_shared, category = c("EAC", "GBM"), euler.d = TRUE, alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), cex=c(3,3,3), cat.cex=c(2,2), cat.just=list(c(.5,.15), c(.5, -.1)))
dev.off()

# HypoDMRs
EAC_hypoDMRs_unique=length(tsgs_with_promoter_overlap_to_EAC_hypoDMRs[!tsgs_with_promoter_overlap_to_EAC_hypoDMRs%in%tsgs_with_promoter_overlap_to_GBM_hypoDMRs])
GBM_hypoDMRs_unique=length(tsgs_with_promoter_overlap_to_GBM_hypoDMRs[!tsgs_with_promoter_overlap_to_GBM_hypoDMRs%in%tsgs_with_promoter_overlap_to_EAC_hypoDMRs])
hypoDMRs_shared=length(tsgs_with_promoter_overlap_to_EAC_hypoDMRs[tsgs_with_promoter_overlap_to_EAC_hypoDMRs %in% tsgs_with_promoter_overlap_to_GBM_hypoDMRs])
grid.newpage()
col1="#FF0000"
col2="#0000FF"
draw.pairwise.venn(area1 = EAC_hypoDMRs_unique+hypoDMRs_shared, area2 = GBM_hypoDMRs_unique+hypoDMRs_shared, cross.area = hypoDMRs_shared, category = c("EAC", "GBM"), euler.d = TRUE, alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), cex=c(3,3,3), cat.cex=c(2,2), cat.just=list(c(.5,.15), c(.5, -.1)))
dev.off()
