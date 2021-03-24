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
shared_hyperDMRs = read.table("shared_EAC_GBM_hyperDMRs_one_line_one_DMR_sorted.bed")
num_eac_hyperDMRs=dim(eac_hyperDMRs)[1]
num_gbm_hyperDMRs=dim(gbm_hyperDMRs)[1]
num_shared_hyperDMRs=dim(shared_hyperDMRs)[1]
grid.newpage()
draw.pairwise.venn(area1 = num_eac_hyperDMRs, area2 = num_gbm_hyperDMRs, cross.area = num_shared_hyperDMRs, , fill = c("red", "seagreen"), euler.d = TRUE, alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), cex=c(3,3,3), cat.cex=c(2,2), cat.just=list(c(.5,.15), c(.5, -.1)))

# HypoDMRs
eac_hypoDMRs = read.table("EAC_hypoDMRs.bed")
gbm_hypoDMRs = read.table("GBM_hypoDMRs.bed")
shared_hypoDMRs = read.table("shared_EAC_GBM_hypoDMRs_one_line_one_DMR_sorted.bed")
num_eac_hypoDMRs=dim(eac_hypoDMRs)[1]
num_gbm_hypoDMRs=dim(gbm_hypoDMRs)[1]
num_shared_hypoDMRs=dim(shared_hypoDMRs)[1]
grid.newpage()
draw.pairwise.venn(area1 = num_eac_hypoDMRs, area2 = num_gbm_hypoDMRs, cross.area = num_shared_hypoDMRs, , fill = c("lightblue", "plum"), euler.d = TRUE, alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), cex=c(3,3,3), cat.cex=c(2,2), cat.just=list(c(.5,.15), c(.5, -.1)))
