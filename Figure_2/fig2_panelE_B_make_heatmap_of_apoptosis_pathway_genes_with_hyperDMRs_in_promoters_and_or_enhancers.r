# Jennifer Flynn jaflynn@wustl.edu
# January 20, 2020
# Determine which KEGG Apoptosis genes (hsa04210) have DMRs within their promoters and/or associated enhancers

# Load in the appropriate library(s):
library(bedr)
library(reshape2)
library(ggplot2)
library(gplots)

# Read in the data 
gene_list=read.table("Apoptosis_Gene_List") #apoptosis gene list (from KEGG pathway "Apoptosis - Homo sapiens", entry: hsa04210)
promoters=read.table("refGene_Promoters2500bp_withGeneNames_sorted.bed")
active_enhancers=read.table("active_enhancer_sites_sorted.bed") #Output from fig2_panelA_determine_epi_defined_promoters_and_enhancers.sh
tss=read.table("refGene_TSSs_withGeneNames_sorted.bed", sep="\t")

# Load in the merged DMR files
EAC_hyperDMRs=read.table("EAC_merged_hyperDMRs.bed")
EAC_hypoDMRs=read.table("EAC_merged_hypoDMRs.bed")
GBM_hyperDMRs=read.table("GBM_merged_hyperDMRs.bed")
GBM_hypoDMRs=read.table("GBM_merged_hypoDMRs.bed")

# Determine which Apoptosis pathway genes have a DMR in their genomic promoters
apoptosis_gene_promoters=unique(promoters[promoters$V5%in%gene_list$V1,])
colnames(apoptosis_gene_promoters)=c("chr", "start", "end", "strand", "gene_id")
apoptosis_gene_promoters$chr=as.character(apoptosis_gene_promoters$chr)
colnames(EAC_hyperDMRs)=c("chr", "start", "end", "value")
EAC_hyperDMRs$chr=as.character(EAC_hyperDMRs$chr)
apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs=bedr(input = list(a = apoptosis_gene_promoters[is.valid.region(apoptosis_gene_promoters),], b = EAC_hyperDMRs), method="intersect", params="-c")
apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs=apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs[apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs$V6>0,]
apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs=unique(apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs$V5)
colnames(GBM_hyperDMRs)=c("chr", "start", "end", "value")
GBM_hyperDMRs$chr=as.character(GBM_hyperDMRs$chr)
apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs=bedr(input = list(a = apoptosis_gene_promoters[is.valid.region(apoptosis_gene_promoters),], b = GBM_hyperDMRs), method="intersect", params="-c")
apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs=apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs[apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs$V6>0,]
apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs=unique(apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs$V5)

# Determine which Apoptosis pathway genes have a DMR in an associated active enhancer that is not within 2.5kb promoters
active_enhancers=unique(data.frame(chr=active_enhancers$V1, start=active_enhancers$V2, stop=active_enhancers$V3))
colnames(promoters)=c("chr", "start", "end")
active_enhancers$chr=as.character(active_enhancers$chr)
promoters$chr=as.character(promoters$chr)
active_enhancers_outside_2.5kb_promoters=bedr(input = list(a = active_enhancers[is.valid.region(active_enhancers),], b = promoters[is.valid.region(promoters),]), method="subtract", params="-A")
tss_bed=data.frame(chr=tss$V1, start=tss$V2-1, end=tss$V2, gene=tss$V4)
tss_bed$chr=as.character(tss_bed$chr)
active_enhancers_not_promoters_gene_assignment=bedr(input = list(a = active_enhancers_outside_2.5kb_promoters, b = tss_bed[is.valid.region(tss_bed),]), method="closest", params="-d")
active_enhancers_not_promoters_gene_assignment$V8=as.numeric(active_enhancers_not_promoters_gene_assignment$V8)
active_enhancers_not_promoters_gene_assignment_500k=active_enhancers_not_promoters_gene_assignment[active_enhancers_not_promoters_gene_assignment$V8<500000,]
active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap=bedr(input = list(a = active_enhancers_not_promoters_gene_assignment_500k, b = EAC_hyperDMRs), method="intersect", params="-c")
active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap$V9=as.numeric(active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap$V9)
active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap=active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap[active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap$V9>0,]
active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap=unique(active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap$V7)
active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap=bedr(input = list(a = active_enhancers_not_promoters_gene_assignment_500k, b = GBM_hyperDMRs), method="intersect", params="-c")
active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap$V9=as.numeric(active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap$V9)
active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap=active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap[active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap$V9>0,]
active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap=unique(active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap$V7)
apoptosis_genes_with_active_enhancer_overlap_to_EAC_hyperDMRs=active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap[active_enhancers_not_promoters_gene_assignment_500k_EAC_hyperDMR_overlap%in%gene_list$V1]
apoptosis_genes_with_active_enhancer_overlap_to_GBM_hyperDMRs=active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap[active_enhancers_not_promoters_gene_assignment_500k_GBM_hyperDMR_overlap%in%gene_list$V1]

# Plot the data
apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs_df=data.frame(gene=apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs, apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs=1)
apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs_df=data.frame(gene=apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs, apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs=1)
apoptosis_genes_with_active_enhancer_overlap_to_EAC_hyperDMRs_df=data.frame(gene=apoptosis_genes_with_active_enhancer_overlap_to_EAC_hyperDMRs, apoptosis_genes_with_active_enhancer_overlap_to_EAC_hyperDMRs=1)
apoptosis_genes_with_active_enhancer_overlap_to_GBM_hyperDMRs_df=data.frame(gene=apoptosis_genes_with_active_enhancer_overlap_to_GBM_hyperDMRs, apoptosis_genes_with_active_enhancer_overlap_to_GBM_hyperDMRs=1)
data_for_plotting=merge(apoptosis_genes_with_promoter_overlap_to_GBM_hyperDMRs_df, apoptosis_genes_with_promoter_overlap_to_EAC_hyperDMRs_df, by="gene", all=TRUE)
data_for_plotting=merge(data_for_plotting, apoptosis_genes_with_active_enhancer_overlap_to_EAC_hyperDMRs_df, by="gene", all=TRUE)
data_for_plotting=merge(data_for_plotting, apoptosis_genes_with_active_enhancer_overlap_to_GBM_hyperDMRs_df, by="gene", all=TRUE)
data_for_plotting[is.na(data_for_plotting)] <- 0
rownames(data_for_plotting)=data_for_plotting$gene
data_for_plotting$gene=NULL
data_for_plotting_matrix=as.matrix(data_for_plotting)
heatmap.2(data_for_plotting_matrix, trace="none", key = FALSE, dendrogram="none", Colv=FALSE, col=c("gray", "blue"))
