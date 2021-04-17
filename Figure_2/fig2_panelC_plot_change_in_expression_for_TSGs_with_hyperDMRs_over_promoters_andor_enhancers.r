# Jennifer Flynn jaflynn@wustl.edu
# June 22, 2020
# This script reads in a list of core promoter locations, lists of DMRs, a list of active enhancer locations, and x.

# Load in the appropriate library(s):
library(bedr)
library(ggplot2)

# Read in the data
EAC_hyperDMRs=read.table("EAC_hyperDMRs.bed")
GBM_hyperDMRs=read.table("GBM_hyperDMRs.bed")
full_tumor_suppressor_gene_file=read.table("Human_TSGs.txt", sep="\t", header=TRUE) # File downloaded from from "http://bioinfo.mc.vanderbilt.edu/TSGene/Human_TSGs.txt" on 02/25/2016
core_promoters_1000bp=read.table("refGene_corePromoters1000bp_withGeneNames_sorted.bed")
tss=read.table("refGene_TSSs_withGeneNames_sorted.bed", sep="\t")
active_enhancers=read.table("active_enhancer_sites_sorted_merged.bed")
eac_tcga_expression=read.table("RPKM_Exp_List_MSI-H.EndoGr3_second_SLC35E2_renamed", row.names = 1, header=T, sep="\t") #EAC expression (TCGA)
normal_endo_expression=read.table("RPKM_Exp_List_Normal_second_SLC35E2_renamed", row.names = 1, header=T, sep="\t") #normal endometrium expression (GTEx)
gbm_tcga_expression=read.table("brain_normal_and_cancer_seperated/GBM_TCGA.Normalized_SLC35E2_second_renamed", row.names = 1, header=T, sep="\t") #GBM expression data (TCGA)
normal_brain_expression=read.table("normal_cortex_GTEx.Normalized_SLC35E2_second_renamed", row.names = 1, header=T, sep="\t") #Normal brain expression (GTEx)

# Correct aliases for DEC1, TRP53COR, LOC401317, and PTPLAD2
full_tumor_suppressor_gene_file$GeneSymbol=as.character(full_tumor_suppressor_gene_file$GeneSymbol)
full_tumor_suppressor_gene_file$GeneSymbol[full_tumor_suppressor_gene_file$GeneSymbol=="Dec-01"]="DEC1"
full_tumor_suppressor_gene_file$GeneSymbol[full_tumor_suppressor_gene_file$GeneSymbol=="PTPLAD2"]="HACD4"
full_tumor_suppressor_gene_file$GeneSymbol[full_tumor_suppressor_gene_file$GeneSymbol=="TRP53COR1"]="TP53COR1"
full_tumor_suppressor_gene_file$GeneSymbol[full_tumor_suppressor_gene_file$GeneSymbol=="LOC401317"]="CREB5"

# Make TSS bed file
tss_bed=data.frame(chr=tss$V1, start=tss$V2-1, end=tss$V2, gene=tss$V4)

# Determine core promoters for TSGs
tsgs=unique(full_tumor_suppressor_gene_file$GeneSymbol)
tsgs_core_promoters=unique(core_promoters_1000bp[core_promoters_1000bp$V5%in%tsgs,])

# Determine which active enhancers do not overlap core promoters
colnames(active_enhancers)=c("chr", "start", "end")
colnames(core_promoters_1000bp)=c("chr", "start", "end", "strand", "gene_id")
active_enhancers$chr=as.character(active_enhancers$chr)
core_promoters_1000bp$chr=as.character(core_promoters_1000bp$chr)
active_enhancers_outside_core_promoters=bedr(input = list(a = active_enhancers[is.valid.region(active_enhancers),], b = core_promoters_1000bp[is.valid.region(core_promoters_1000bp),]), method="subtract", params="-A")

# Assign active enhancers that do not overlap promoters to the nearest TSS (within 500,000 bps)
tss_bed$chr=as.character(tss_bed$chr)
active_enhancers_not_promoters_gene_assignment=bedr(input = list(a = active_enhancers_outside_core_promoters, b = tss_bed[is.valid.region(tss_bed),]), method="closest", params="-d")
active_enhancers_not_promoters_gene_assignment$V8=as.numeric(active_enhancers_not_promoters_gene_assignment$V8)
active_enhancers_not_promoters_gene_assignment_500k=active_enhancers_not_promoters_gene_assignment[active_enhancers_not_promoters_gene_assignment$V8<500000,]

# Determine which TSGs have associated nearby active enhancers (outside of core promoters)
tsg_active_enhancers_not_promoters_gene_assignment=active_enhancers_not_promoters_gene_assignment_500k[active_enhancers_not_promoters_gene_assignment_500k$V7%in%tsgs,]

# Determine which TSGs have a hyperDMR overlapping an associated active enhancer
EAC_hyperDMRs$V1=as.character(EAC_hyperDMRs$V1)
tsgs_with_EAC_hyperDMRs_over_active_enhancer=bedr(input = list(a = tsg_active_enhancers_not_promoters_gene_assignment, b = EAC_hyperDMRs), method="intersect", params="-c")
tsgs_with_EAC_hyperDMRs_over_active_enhancer=tsgs_with_EAC_hyperDMRs_over_active_enhancer[tsgs_with_EAC_hyperDMRs_over_active_enhancer$V9>0,]
tsgs_with_EAC_hyperDMRs_over_active_enhancer=unique(tsgs_with_EAC_hyperDMRs_over_active_enhancer$V7)
GBM_hyperDMRs$V1=as.character(GBM_hyperDMRs$V1)
tsgs_with_GBM_hyperDMRs_over_active_enhancer=bedr(input = list(a = tsg_active_enhancers_not_promoters_gene_assignment, b = GBM_hyperDMRs), method="intersect", params="-c")
tsgs_with_GBM_hyperDMRs_over_active_enhancer=tsgs_with_GBM_hyperDMRs_over_active_enhancer[tsgs_with_GBM_hyperDMRs_over_active_enhancer$V9>0,]
tsgs_with_GBM_hyperDMRs_over_active_enhancer=unique(tsgs_with_GBM_hyperDMRs_over_active_enhancer$V7)

# Determine which TSGs have a hyperDMR overlapping their core promoter
colnames(tsgs_core_promoters)=c("chr", "start", "end", "strand", "gene_id")
tsgs_core_promoters$chr=as.character(tsgs_core_promoters$chr)
tsgs_with_EAC_hyperDMRs_over_core_promoter=bedr(input = list(a = tsgs_core_promoters[is.valid.region(tsgs_core_promoters),], b = EAC_hyperDMRs), method="intersect", params="-c")
tsgs_with_EAC_hyperDMRs_over_core_promoter=tsgs_with_EAC_hyperDMRs_over_core_promoter[tsgs_with_EAC_hyperDMRs_over_core_promoter$V6>0,]
tsgs_with_EAC_hyperDMRs_over_core_promoter=unique(tsgs_with_EAC_hyperDMRs_over_core_promoter$V5)
tsgs_with_GBM_hyperDMRs_over_core_promoter=bedr(input = list(a = tsgs_core_promoters[is.valid.region(tsgs_core_promoters),], b = GBM_hyperDMRs), method="intersect", params="-c")
tsgs_with_GBM_hyperDMRs_over_core_promoter=tsgs_with_GBM_hyperDMRs_over_core_promoter[tsgs_with_GBM_hyperDMRs_over_core_promoter$V6>0,]
tsgs_with_GBM_hyperDMRs_over_core_promoter=unique(tsgs_with_GBM_hyperDMRs_over_core_promoter$V5)

# Create dataframes listing which TSGs contain DMRs within promoters and/or enhancers for plotting for EAC
tsgs_with_EAC_hyperDMR_in_promoter=data.frame(gene=tsgs_with_EAC_hyperDMRs_over_core_promoter, tsgs_with_EAC_hyperDMR_in_promoter=1)
tsgs_with_EAC_hyperDMR_in_active_enhancer=data.frame(gene=tsgs_with_EAC_hyperDMRs_over_active_enhancer, tsgs_with_EAC_hyperDMRs_over_active_enhancer=1)
EAC_genes=merge(tsgs_with_EAC_hyperDMR_in_promoter, tsgs_with_EAC_hyperDMR_in_active_enhancer, by="gene", all=T)

# Determine the change in expression for each gene with an EAC hyperDMR in the TSG's promoter and/or enhancer
eac_tcga_expression$gene=rownames(eac_tcga_expression)
normal_endo_expression$gene=rownames(normal_endo_expression)
merged_eac_endo_data=merge(eac_tcga_expression, normal_endo_expression, by="gene")
rownames(merged_eac_endo_data)=merged_eac_endo_data$gene
merged_eac_endo_data$gene=NULL
merged_eac_endo_data_matrix=as.matrix(merged_eac_endo_data)
# Replace aliases
EAC_genes$gene=as.character(EAC_genes$gene)
EAC_genes$gene[EAC_genes$gene=="BRINP1"]="DBC1"
EAC_genes$gene[EAC_genes$gene=="LOC401317"]="CREB5"
EAC_genes$gene[EAC_genes$gene=="HCAR2"]="GPR109A"
EAC_genes$gene[EAC_genes$gene=="KDM8"]="JMJD5"
EAC_genes$gene[EAC_genes$gene=="NEURL1"]="NEURL"
EAC_genes$gene[EAC_genes$gene=="PLK5"]="PLK5P"
EAC_genes$gene[EAC_genes$gene=="RITA1"]="C12orf52"
EAC_genes$gene[EAC_genes$gene=="ZBTB18"]="ZNF238"
EAC_hyperDMR_TSGs=unique(EAC_genes$gene)
# Calculate log2 expression changes
EAC_hyperDMR_TSGs_expression_data=merged_eac_endo_data_matrix[rownames(merged_eac_endo_data_matrix)%in%EAC_hyperDMR_TSGs,]
EAC_hyperDMR_TSGs_expression_data=as.matrix(EAC_hyperDMR_TSGs_expression_data)
EAC_hyperDMR_TSGs_expression_changes=apply(EAC_hyperDMR_TSGs_expression_data, 1, function(x) if( is.na(log2(mean(x[1:57]) / mean(x[58-86])))) log2((mean(x[1:57])+.000001) / (mean(x[58-86])+.000001)) else (log2(mean(x[1:57]) / mean(x[58-86]))))
EAC_hyperDMR_TSGs_expression_changes_df=as.data.frame(EAC_hyperDMR_TSGs_expression_changes)
# Plot the data
ggplot(EAC_hyperDMR_TSGs_expression_changes_df, aes(x = reorder(EAC_hyperDMR_TSGs_expression_changes, EAC_hyperDMR_TSGs_expression_changes), y = EAC_hyperDMR_TSGs_expression_changes)) + 
  geom_bar(stat = "identity") + ylim(-5,1)
dev.off()
# Get the order of the genes plotted 
EAC_hyperDMR_TSGs_expression_changes_df$gene=rownames(EAC_hyperDMR_TSGs_expression_changes_df)
EAC_hyperDMR_TSGs_expression_changes_df_ordered <- EAC_hyperDMR_TSGs_expression_changes_df[order(EAC_hyperDMR_TSGs_expression_changes_df$EAC_hyperDMR_TSGs_expression_changes),]
EAC_genes_plotted=merge(EAC_genes, EAC_hyperDMR_TSGs_expression_changes_df_ordered, by="gene")
EAC_genes_plotted_ordered=EAC_genes_plotted[order(EAC_genes_plotted$EAC_hyperDMR_TSGs_expression_changes),]
EAC_genes_plotted_ordered[is.na(EAC_genes_plotted_ordered)]=0
EAC_genes_plotted_ordered_plotting=EAC_genes_plotted_ordered[,2:3]
EAC_genes_plotted_ordered_matrix=as.matrix(EAC_genes_plotted_ordered_plotting)
heatmap.2(EAC_genes_plotted_ordered_matrix, trace = "none", Colv = FALSE, Rowv = FALSE, col = c("gray", "orange"))
dev.off()

# Create dataframes listing which TSGs contain DMRs within promoters and/or enhancers for plotting for GBM
tsgs_with_GBM_hyperDMR_in_promoter=data.frame(gene=tsgs_with_GBM_hyperDMRs_over_core_promoter, tsgs_with_GBM_hyperDMR_in_promoter=1)
tsgs_with_GBM_hyperDMR_in_active_enhancer=data.frame(gene=tsgs_with_GBM_hyperDMRs_over_active_enhancer, tsgs_with_GBM_hyperDMRs_over_active_enhancer=1)
GBM_genes=merge(tsgs_with_GBM_hyperDMR_in_promoter, tsgs_with_GBM_hyperDMR_in_active_enhancer, by="gene", all=T)

# Determine the change in expression for each gene with a GBM hyperDMR in the TSG's promoter and/or enhancer
gbm_tcga_expression$gene=rownames(gbm_tcga_expression)
normal_brain_expression$gene=rownames(normal_brain_expression)
merged_gbm_brain_data=merge(gbm_tcga_expression, normal_brain_expression, by="gene")
rownames(merged_gbm_brain_data)=merged_gbm_brain_data$gene
merged_gbm_brain_data$gene=NULL
merged_gbm_brain_data_matrix=as.matrix(merged_gbm_brain_data)
# Calculate log2 expression changes
GBM_hyperDMR_TSGs_expression_data=merged_gbm_brain_data_matrix[rownames(merged_gbm_brain_data_matrix)%in%GBM_genes$gene,]
GBM_hyperDMR_TSGs_expression_data=as.matrix(GBM_hyperDMR_TSGs_expression_data)
GBM_hyperDMR_TSGs_expression_changes=apply(GBM_hyperDMR_TSGs_expression_data, 1, function(x) log2(mean(x[1:160]) / mean(x[161-188])))
GBM_hyperDMR_TSGs_expression_changes_df=as.data.frame(GBM_hyperDMR_TSGs_expression_changes)
# Plot the data
ggplot(GBM_hyperDMR_TSGs_expression_changes_df, aes(x = reorder(GBM_hyperDMR_TSGs_expression_changes, GBM_hyperDMR_TSGs_expression_changes), y = GBM_hyperDMR_TSGs_expression_changes)) + 
  geom_bar(stat = "identity") + ylim(-5,1)
dev.off()
# Get the order of the genes plotted 
GBM_hyperDMR_TSGs_expression_changes_df$gene=rownames(GBM_hyperDMR_TSGs_expression_changes_df)
GBM_hyperDMR_TSGs_expression_changes_df_ordered <- GBM_hyperDMR_TSGs_expression_changes_df[order(GBM_hyperDMR_TSGs_expression_changes_df$GBM_hyperDMR_TSGs_expression_changes),]
GBM_genes_plotted=merge(GBM_genes, GBM_hyperDMR_TSGs_expression_changes_df_ordered, by="gene")
GBM_genes_plotted_ordered=GBM_genes_plotted[order(GBM_genes_plotted$GBM_hyperDMR_TSGs_expression_changes),]
GBM_genes_plotted_ordered[is.na(GBM_genes_plotted_ordered)]=0
GBM_genes_plotted_ordered_plotting=GBM_genes_plotted_ordered[,2:3]
GBM_genes_plotted_ordered_matrix=as.matrix(GBM_genes_plotted_ordered_plotting)
heatmap.2(GBM_genes_plotted_ordered_matrix, trace = "none", Colv = FALSE, Rowv = FALSE, col = c("gray", "red"))
dev.off()
