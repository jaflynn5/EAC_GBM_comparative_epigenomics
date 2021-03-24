# Jennifer Flynn jaflynn@wustl.edu
# September 24, 2015
# This script reads in EAC and GBM DMRs and creates venn diagrams displaying their overlaps

# Load the packages
library(limma)
library(gplots)
library(venneuler)
library(VennDiagram)

# HyperDMRs
eac_hyperDMRs = read.table("EAC_hyperDMRs.bed")
gbm_hyperDMRs = read.table("GBM_hyperDMRs.bed")
shared_hyperDMRs = read.table("shared_EAC_GBM_hyperDMRs.bed")
num_eac_hyperDMRs=dim(eac_hyperDMRs)[1]
num_gbm_hyperDMRs=dim(gbm_hyperDMRs)[1]
num_shared_hyperDMRs=dim(shared_hyperDMRs)[1]
grid.newpage()
draw.pairwise.venn(area1 = num_eac_hyperDMRs, area2 = num_gbm_hyperDMRs, cross.area = num_shared_hyperDMRs, , fill = c("red", "seagreen"), euler.d = TRUE, alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), cex=c(3,3,3), cat.cex=c(2,2), cat.just=list(c(.5,.15), c(.5, -.1)))
phyper(num_shared_hyperDMRs-1, num_gbm_hyperDMRs, 5196471-num_gbm_hyperDMRs, num_eac_hyperDMRs, lower.tail = FALSE) #5196471: number of 500bp genomic regions with at least 1 CpG and MeDIP and/or MRE signal in at least one sample; please see script: fig1_panelA_determine_background_bins.sh

# HypoDMRs
eac_hypoDMRs = read.table("EAC_hypoDMRs.bed")
gbm_hypoDMRs = read.table("GBM_hypoDMRs.bed")
shared_hypoDMRs = read.table("shared_EAC_GBM_hypoDMRs.bed")
num_eac_hypoDMRs=dim(eac_hypoDMRs)[1]
num_gbm_hypoDMRs=dim(gbm_hypoDMRs)[1]
num_shared_hypoDMRs=dim(shared_hypoDMRs)[1]
grid.newpage()
draw.pairwise.venn(area1 = num_eac_hypoDMRs, area2 = num_gbm_hypoDMRs, cross.area = num_shared_hypoDMRs, , fill = c("lightblue", "plum"), euler.d = TRUE, alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), cex=c(3,3,3), cat.cex=c(2,2), cat.just=list(c(.5,.15), c(.5, -.1)))
phyper(num_shared_hypoDMRs-1, num_gbm_hypoDMRs, 5196471-num_gbm_hypoDMRs, num_eac_hypoDMRs, lower.tail = FALSE) #5196471: number of 500bp genomic regions with at least 1 CpG and MeDIP and/or MRE signal in at least one sample; please see script: fig1_panelA_determine_background_bins.sh
