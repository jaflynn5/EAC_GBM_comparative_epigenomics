# Jennifer Flynn jaflynn@wustl.edu
# March 19, 2020
# This script reads in expression data for EAC and normal endometruim and determines the expression z-scores for genes with enriched TFB motifs in ceDMRs (q<0.01 + present in at least 20% of ceDMRs) (motif enrichment was calculated using Homer v4.9 and opposite DMR groups as background)

##### Figure 4 Panels B and C #####
### Read in the Expression Data
# Load in the EAC expression data (TCGA)
eac_tcga_expression=read.table("RPKM_Exp_List_MSI-H.EndoGr3_second_SLC35E2_renamed", row.names = 1, header=T, sep="\t")
eac_tcga_expression$gene=rownames(eac_tcga_expression)

# Load in the normal endometrium expression (GTEx)
normal_endo_expression=read.table("RPKM_Exp_List_Normal_second_SLC35E2_renamed", row.names = 1, header=T, sep="\t")
normal_endo_expression$gene=rownames(normal_endo_expression)

# Merge the data
merged_eac_endo_data=merge(eac_tcga_expression, normal_endo_expression, by="gene")

# Z-score the data
# Step 1): Remove any genes with an average RPKM < 1:
# EAC, normal endo
rownames(merged_eac_endo_data)=merged_eac_endo_data$gene
merged_eac_endo_data$gene=NULL
merged_eac_endo_data_matrix=as.matrix(merged_eac_endo_data)
merged_eac_endo_data_matrix_expressed=merged_eac_endo_data_matrix[rowMeans(merged_eac_endo_data_matrix)>1,]

# Step 2): Add pseudocounts (I made the pseudocount 1/100th of the value of the smallest value)
# EAC, normal endo
smallest_value_eac_endo=min(merged_eac_endo_data_matrix_expressed[merged_eac_endo_data_matrix_expressed!=0])
pseudocount=smallest_value_eac_endo/100
merged_eac_endo_data_matrix_expressed[merged_eac_endo_data_matrix_expressed==0]=pseudocount

# Step 3): Log transform the data (see: https://www.researchgate.net/post/How_can_I_calculate_z-score_from_rpkm_or_counts_values)
# EAC, normal endo
merged_eac_endo_data_matrix_expressed_log2=log2(merged_eac_endo_data_matrix_expressed)

### EAC ce-hyperDMR enriched motifs ###
# Step 4): Extract the genes of interest (output from fig3_panelA_plot_percentage_of_ceDMRs_containing_significantly_enriched_motifs_only_showing_motifs_in_at_least_20_percent_of_ceDMRs.r)
EAC_ce_hyperDMR_enriched_motifs=c("EGR1", "TCF21", "ZFX", "ZNF711", "OLIG2", "TFAP4", "TBP", "ATOH1", "MYOD1", "MAZ", "NEUROG2", "ZNF263", "ZBTB7A", 
                                  "TAL1", "TCF12", "MYF5", "ZNF467", "RBPJ", "NEUROD1", "SMAD2", "SMAD3", "MAFA", "MAX", "GATA3", "RXRA", "PGR",
                                  "KLF14", "SMAD4", "FLI1", "MYCN", "MYBL1", "ETS1", "MYBL2", "NPAS1", "TEAD4", "ARNTL", "ETV1", "NR2F2", "MYOG", 
                                  "ETV2", "STAT4", "PPARG", "FOXO1", "BCL6")
EAC_ce_hyperDMR_enriched_motifs=unique(EAC_ce_hyperDMR_enriched_motifs)
EAC_ce_hyperDMR_enriched_motifs_expression_data=merged_eac_endo_data_matrix_expressed_log2[rownames(merged_eac_endo_data_matrix_expressed_log2)%in%EAC_ce_hyperDMR_enriched_motifs,]
EAC_ce_hyperDMR_enriched_motifs_expression_data=as.matrix(EAC_ce_hyperDMR_enriched_motifs_expression_data)
EAC_ce_hyperDMR_enriched_motifs_expression_data_z_scores=apply(EAC_ce_hyperDMR_enriched_motifs_expression_data, 1, function(x) ((x-mean(x))/sd(x)))
EAC_ce_hyperDMR_enriched_motifs_expression_data_z_scores_df=as.data.frame(EAC_ce_hyperDMR_enriched_motifs_expression_data_z_scores)
EAC_ce_hyperDMR_enriched_motifs_expression_data_z_scores_df$sample=rownames(EAC_ce_hyperDMR_enriched_motifs_expression_data_z_scores_df)
EAC_ce_hyperDMR_enriched_motifs_expression_data_z_scores_df$sample_type=c(rep("EAC", length(colnames(eac_tcga_expression))-1), rep("Normal_Endo", length(colnames(normal_endo_expression))-1))
EAC_ce_hyperDMR_enriched_motifs_expression_data_z_scores_df_melt=melt(EAC_ce_hyperDMR_enriched_motifs_expression_data_z_scores_df)

# Step 5): Plot the data
ggplot(EAC_ce_hyperDMR_enriched_motifs_expression_data_z_scores_df_melt, aes(x=reorder(sample, value, FUN = median), y=value)) + 
  geom_boxplot(aes(fill=reorder(sample_type, value, FUN = median)))  + scale_fill_manual(values=c("#C791CB", "#91CBA6")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

### EAC ce-hypoDMR enriched motifs ###
# Step 4): Extract the genes of interest (output from fig3_panelA_plot_percentage_of_ceDMRs_containing_significantly_enriched_motifs_only_showing_motifs_in_at_least_20_percent_of_ceDMRs.r)
EAC_ce_hypoDMR_enriched_motifs=c("ZEB1", "TCF3", "SNAI2", "TBX5", "FOXM1", "SOX2", "NKX3-2", "SOX15",
                                 "FOXA1", "FOXA2", "SOX10", "SOX3", "FOXA1", "TCF3", "SOX6", "NKX2-1",
                                 "NKX2-2", "GSC", "TCF12", "KLF5", "LHX2", "SIX2", "NKX2-5", "NFIX",
                                 "EOMES", "ZNF416", "TGIF2", "NKX3-1", "CRX", "TBX21", "SOX9", "LHX1")
EAC_ce_hypoDMR_enriched_motifs=unique(EAC_ce_hypoDMR_enriched_motifs)
EAC_ce_hypoDMR_enriched_motifs_expression_data=merged_eac_endo_data_matrix_expressed_log2[rownames(merged_eac_endo_data_matrix_expressed_log2)%in%EAC_ce_hypoDMR_enriched_motifs,]
EAC_ce_hypoDMR_enriched_motifs_expression_data=as.matrix(EAC_ce_hypoDMR_enriched_motifs_expression_data)
EAC_ce_hypoDMR_enriched_motifs_expression_data_z_scores=apply(EAC_ce_hypoDMR_enriched_motifs_expression_data, 1, function(x) ((x-mean(x))/sd(x)))
EAC_ce_hypoDMR_enriched_motifs_expression_data_z_scores_df=as.data.frame(EAC_ce_hypoDMR_enriched_motifs_expression_data_z_scores)
EAC_ce_hypoDMR_enriched_motifs_expression_data_z_scores_df$sample=rownames(EAC_ce_hypoDMR_enriched_motifs_expression_data_z_scores_df)
EAC_ce_hypoDMR_enriched_motifs_expression_data_z_scores_df$sample_type=c(rep("EAC", length(colnames(eac_tcga_expression))-1), rep("Normal_Endo", length(colnames(normal_endo_expression))-1))
EAC_ce_hypoDMR_enriched_motifs_expression_data_z_scores_df_melt=melt(EAC_ce_hypoDMR_enriched_motifs_expression_data_z_scores_df)

# Step 5): Plot the data
ggplot(EAC_ce_hypoDMR_enriched_motifs_expression_data_z_scores_df_melt, aes(x=reorder(sample, value, FUN = median), y=value)) + 
  geom_boxplot(aes(fill=reorder(sample_type, value, FUN = median)))  + scale_fill_manual(values=c("#91CBA6", "#C791CB")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
