# Jennifer Flynn jaflynn@wustl.edu
# April 22, 2020
# This script reads in list of genes with significantly enriched motifs in ce-DMRs, determines if their expression significantly changes, and plots the data. 

# Load in appropriate libraries #
library(gplots)
library(RColorBrewer)

# Define functions
calculate_expression_change_direction <- function(motif_list, merged_expression_log2, cancer_expression, normal_expression){
  motifs=unique(motif_list)
  motifs_data=merged_expression_log2[rownames(merged_expression_log2)%in%motifs,]
  direction=apply(motifs_data, MARGIN=1, FUN=function(x) if(mean(x[1:length(colnames(cancer_expression))-1])>mean(x[length(colnames(cancer_expression)):(length(colnames(normal_expression))+length(colnames(cancer_expression))-2)])) 1 else 0)
  return(direction)
}

calculate_expression_change_pval <- function(motif_list, merged_expression_log2, cancer_expression, normal_expression){
  motifs=unique(motif_list)
  motifs_data=merged_expression_log2[rownames(merged_expression_log2)%in%motifs,]
  p_values=apply(motifs_data, MARGIN=1, FUN=function(x) wilcox.test(x[1:length(colnames(cancer_expression))-1], x[length(colnames(cancer_expression)):(length(colnames(normal_expression))+length(colnames(cancer_expression))-2)])[[3]])
  p_adj=p.adjust(p_values, method = "BH", n = length(p_values))
  return(p_adj)
}

split_string <- function(string_to_split){
  unlist(strsplit(string_to_split, "[.]"))[[1]][1]
}

process_data <- function(my.data){
  colnames(my.data)=sapply(colnames(my.data), split_string)
  rownames_data=paste0(my.data$Chr, "_", my.data$Start, "_", my.data$End)
  my.data$PeakID=NULL
  my.data$Chr=NULL
  my.data$Start=NULL
  my.data$End=NULL
  columnnames=colnames(my.data)
  my.data <- data.frame(lapply(my.data, as.character), stringsAsFactors=FALSE)
  colnames(my.data)=columnnames
  my.data[my.data!=""]=1
  my.data[my.data==""]=0
  my.data=as.matrix(my.data)
  my.data=apply(my.data, 2, as.numeric)
  rownames(my.data)=rownames_data
  return(my.data)
}

plot_data_clusters <- function(data, clust_dist, color)
{
  # Plotting data, clusters defined
  par(mar = rep(0, 4))
  a=heatmap.2(t(data), trace="none", col = c("gray", color), dendrogram = "row", 
              distfun = function(x) dist(x, method="euclidean"),
              hclustfun = function(x) hclust(x, method="complete"))
  hc <- as.hclust( a$rowDendrogram )
  row_groups <- cutree(hc, h=clust_dist) 
  return(row_groups)
}


### Read in the Expression Data
# Load in the EAC expression data (TCGA)
eac_tcga_expression=read.table("RPKM_Exp_List_MSI-H.EndoGr3_second_SLC35E2_renamed", row.names = 1, header=T, sep="\t")
eac_tcga_expression$gene=rownames(eac_tcga_expression)
# Load in the normal endometrium expression (GTEx)
normal_endo_expression=read.table("RPKM_Exp_List_Normal_second_SLC35E2_renamed", row.names = 1, header=T, sep="\t")
normal_endo_expression$gene=rownames(normal_endo_expression)
# Load in the GBM expression data (TCGA)
gbm_tcga_expression=read.table("GBM_TCGA.Normalized_SLC35E2_second_renamed", row.names = 1, header=T, sep="\t")
gbm_tcga_expression$gene=rownames(gbm_tcga_expression)
# Load in the normal brain expression (GTEx)
normal_brain_expression=read.table("normal_cortex_GTEx.Normalized_SLC35E2_second_renamed", row.names = 1, header=T, sep="\t")
normal_brain_expression$gene=rownames(normal_brain_expression)

# Merge the data
merged_eac_endo_data=merge(eac_tcga_expression, normal_endo_expression, by="gene")
merged_gbm_brain_data=merge(gbm_tcga_expression, normal_brain_expression, by="gene")
rownames(merged_eac_endo_data)=merged_eac_endo_data$gene
merged_eac_endo_data$gene=NULL
merged_eac_endo_data_matrix=as.matrix(merged_eac_endo_data)
rownames(merged_gbm_brain_data)=merged_gbm_brain_data$gene
merged_gbm_brain_data$gene=NULL
merged_gbm_brain_data_matrix=as.matrix(merged_gbm_brain_data)

# Add pseudocounts (I made the pseudocount 1/100th of the value of the smallest value)
smallest_value_eac_endo=min(merged_eac_endo_data_matrix[merged_eac_endo_data_matrix!=0])
pseudocount=smallest_value_eac_endo/100
merged_eac_endo_data_matrix[merged_eac_endo_data_matrix==0]=pseudocount
smallest_value_gbm_brain=min(merged_gbm_brain_data_matrix[merged_gbm_brain_data_matrix!=0])
pseudocount=smallest_value_gbm_brain/100
merged_gbm_brain_data_matrix[merged_gbm_brain_data_matrix==0]=pseudocount

# Log transform the data (see: https://www.researchgate.net/post/How_can_I_calculate_z-score_from_rpkm_or_counts_values)
merged_eac_endo_data_matrix_log2=log2(merged_eac_endo_data_matrix)
merged_gbm_brain_data_matrix_log2=log2(merged_gbm_brain_data_matrix)


##### Figure 5A #####
EAC_ce_hyperDMR_enriched_motifs=c("PITX1", "CDX4", "HOXB13", "EGR1", "CDX2", "TCF21", "ZFX", "ZNF711", "ZNF238", "OLIG2", "HOXD13", "TFAP4", "TBP", "ATOH1", "E2F6",
                                  "ZNF281", "MYOD1", "MAZ", "NEUROG2", "ZNF263", "NRF1", "JUND", "BACH2", "ZBTB7A", "TAL1", "TCF12", "HAND2", "MEF2D", "MYF5", "GATA5",
                                  "JUN", "ZNF467", "RBPJ", "NFE2", "GATA4", "CTCFL", "RFX3", "CTCF", "JUNB", "BATF", "CREB1", "HOXA9", "ATF3", "THAP11", "GATA2",
                                  "NEUROD1", "EGR2", "SMAD2", "FOSL2", "FOSL1", "ZNF528", "E2F1", "TEAD1", "SRF", "HOXC9", "NFATC1", "SMAD3", "ELK1", "E2F3", "RFX2",
                                  "ATF2", "GATA1", "MAFA", "BACH1", "MAX", "GATA6", "USF1", "GATA3", "AR", "TEAD2", "PGR", "KLF14", "SMAD4", "STAT1", "FLI1",
                                  "NR2C2", "ELK4", "MYCN", "MYBL1", "EBF1", "RORC", "ETS1", "MYBL2", "ATF4", "E2F4", "CLOCK", "NPAS1", "TEAD4", "TWIST1", "ZNF182", 
                                  "MYC", "ARNTL", "ETV1", "NFKB1", "NFKB2", "NR2F2", "STAT5A", "NFE2L2", "MAFF", "MYOG", "RFX1", "NR5A2", "ETV2", "ZNF692", "MAFK",
                                  "CEBPB", "ATF7", "STAT4", "PPARG", "PPARA", "RXRA", "RXRB", "RXRG", "RELA", "DDIT3", "RARA", "FOXO1", "BHLHE40", "PRDM4", "BCL6",
                                  "E2F7", "FOXO3", "TCF7L2", "NR4A1", "FOXK1", "ELF3", "RUNX1", "STAT6", "NFATC2", "NFATC3", "NFATC4", "NFAT5", "FOXK2", "MAFB", "REST",
                                  "TFE3", "GABPA", "ESR1", "ESR2", "SP5", "HLF", "ZSCAN22", "VDR", "STAT3", "SPI1", "ELF1", "GLI3") # Extract the genes of interest #Using the homer key (http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults.html), I identified the gene names for the enriched motifs (and used alternate names as indicated in "EAC_cehyperDMR_homer_results_using_EAC_cehypoDMRs_as_background_filtering_at_q_0.05.txt").
EAC_ce_hyperDMR_enriched_motifs_qvals=calculate_expression_change_pval(EAC_ce_hyperDMR_enriched_motifs, merged_eac_endo_data_matrix_log2, eac_tcga_expression, normal_endo_expression)
EAC_ce_hyperDMR_enriched_motifs_direction=calculate_expression_change_direction(EAC_ce_hyperDMR_enriched_motifs, merged_eac_endo_data_matrix_log2, eac_tcga_expression, normal_endo_expression)
EAC_ce_hyperDMR_enriched_motifs_and_sig_loss_exp=data.frame(qvals=EAC_ce_hyperDMR_enriched_motifs_qvals, direction=EAC_ce_hyperDMR_enriched_motifs_direction)
EAC_ce_hyperDMR_enriched_motifs_and_sig_loss_exp=rownames(EAC_ce_hyperDMR_enriched_motifs_and_sig_loss_exp[EAC_ce_hyperDMR_enriched_motifs_and_sig_loss_exp$qvals<0.05 & EAC_ce_hyperDMR_enriched_motifs_and_sig_loss_exp$direction==0,])
data=read.table("homerMotifs.enriched_and_down_regulated.motifs.locations_EAC_hyper.txt", header=T, sep="\t") #See key (EAC_cehyperDMR_enriched_motifs_key.txt) for identifying columns to include from the homer output locations file.
data=process_data(data)
cluster_distance=42
row_groupings=plot_data_clusters(data, cluster_distance, "firebrick3") # From this output (row_groupings), I could then identify clusters of TFs. I called a cluster any cluster that had at least 2 TFs in it. The locations of these clusters were used as input to GREAT. See manuscript for more details.

##### Figure 5B #####
GBM_ce_hyperDMR_enriched_motifs=c("ARNT", "BHLHE40", "CLOCK", "CTCF", "CTCFL", "CUX2", "E2F1", "E2F4", "E2F6", "EGR1", "EGR2", "EPAS1", "FOXA1", "FOXF1", "FOXO3", 
                                  "GATA3", "GATA6", "GSC", "HNF1B", "HOXA11", "HOXB13", "IRF9", "KLF3", "KLF9", "MAFB", "MAFF", "MEF2A", "MEF2C", "MEF2D",
                                  "MITF", "MYB", "MYBL2", "MYC", "NFATC1", "NFKB1", "NR1H2", "NRF1", "POU3F1", "POU3F2", "SLC22A2", "PAX5", "PAX7", "PHOX2A", "POU1F1", "SP1",
                                  "POU3F3", "SP5", "STAT1", "STAT4", "TBP", "TBR1", "T", "TFE3", "USF1", "USF2", "YY1", "ZBTB33", "ZNF136", "ZNF165", "ZNF281", "ZNF317", "ZNF669", "FOXK1") # Extract the genes of interest #Using the homer key (http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults.html), I identified the gene names for the enriched motifs (and used alternate names as indicated in "GBM_cehyperDMR_homer_results_using_GBM_cehypoDMRs_as_background_filtering_at_q_0.05.txt").
GBM_ce_hyperDMR_enriched_motifs_qvals=calculate_expression_change_pval(GBM_ce_hyperDMR_enriched_motifs, merged_gbm_brain_data_matrix_log2, gbm_tcga_expression, normal_brain_expression)
GBM_ce_hyperDMR_enriched_motifs_direction=calculate_expression_change_direction(GBM_ce_hyperDMR_enriched_motifs, merged_gbm_brain_data_matrix_log2, gbm_tcga_expression, normal_brain_expression)
GBM_ce_hyperDMR_enriched_motifs_and_sig_loss_exp=data.frame(qvals=GBM_ce_hyperDMR_enriched_motifs_qvals, direction=GBM_ce_hyperDMR_enriched_motifs_direction)
GBM_ce_hyperDMR_enriched_motifs_and_sig_loss_exp=rownames(GBM_ce_hyperDMR_enriched_motifs_and_sig_loss_exp[GBM_ce_hyperDMR_enriched_motifs_and_sig_loss_exp$qvals<0.05 & GBM_ce_hyperDMR_enriched_motifs_and_sig_loss_exp$direction==0,])
data=read.table("homerMotifs.enriched_and_down_regulated.motifs.locations_GBM_hyper.txt", header=T, sep="\t") #See key (GBM_cehyperDMR_enriched_motifs_key.txt) for identifying columns to include from the homer output locations file.
data=process_data(data)
cluster_distance=25
row_groupings=plot_data_clusters(data, cluster_distance, "firebrick3") # From this output (row_groupings), I could then identify clusters of TFs. I called a cluster any cluster that had at least 2 TFs in it. The locations of these clusters were used as input to GREAT. See manuscript for more details.

##### Figure 6A #####
EAC_ce_hypoDMR_enriched_motifs=c("AR", "ASCL1", "ATF2", "ATF3", "ATF4", "ATOH1", "BACH2", "BARX1", "BATF", "CEBPA", "CREB1", "DDIT3", "EBF1", "ERG",
                                 "EHF", "ELF1", "ELF3", "ELF5", "ELK4", "ESR1", "ESR2", "ESRRA", "ETS1", "ETV1", "ETV2", "FLI1", "FOS", "FOSL2", "FOXA2", "FOXO1",
                                 "GABPA", "GATA1", "GATA2", "GATA3", "GATA4", "GATA5", "GLI3", "GRHL2", "HAND2", "HNF4A", "HOXA2", "HOXB4", "HOXC9", "HSF1", "IRF8", "ISL1",
                                 "JUN", "JUNB", "JUND", "KLF5", "KLF6", "LHX1", "LHX2", "LHX3", "MAFK", "MYF5", "MYOD1", "MYOG", "NANOG", "NEUROD1",
                                 "NEUROG2", "NFAT5", "NFATC2", "NFATC3", "NFATC4", "NFE2", "NFIA", "NFIX", "NKX6-1", "NR1H4", "NR2E1", "NR5A2", "OLIG2", "PAX7", "PBX1", "PDX1", "PITX1",
                                 "POU3F1", "POU3F3", "PPARA", "PPARG", "PRDM4", "RBPJ", "RELA", "RFX6", "RUNX1", "RUNX2", "RXRA", "RXRB", "RXRG", "SMAD2", "SMAD3", "SNAI2", "SOX10", "SOX15",
                                 "SOX2", "SOX3", "SPDEF", "SPI1", "TAL1","TBX5", "TCF12", "TCF21", "TCF3", "TCF4", "TCF7L1", "TEAD2", "TFAP2A", "TFAP2C", "TFAP4", "TFCP2L1", "TGIF2", "THRB", "TWIST1", "VDR",
                                 "ZBTB7A", "ZEB1", "ZFX", "ZNF16", "ZNF416", "ZNF467", "ZNF692", "ZNF711") # Extract the genes of interest #Using the homer key (http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults.html), I identified the gene names for the enriched motifs (and used alternate names as indicated in "EAC_cehypoDMR_homer_results_using_EAC_cehyperDMRs_as_background_filtering_at_q_0.05.txt").
EAC_ce_hypoDMR_enriched_motifs_qvals=calculate_expression_change_pval(EAC_ce_hypoDMR_enriched_motifs, merged_eac_endo_data_matrix_log2, eac_tcga_expression, normal_endo_expression)
EAC_ce_hypoDMR_enriched_motifs_direction=calculate_expression_change_direction(EAC_ce_hypoDMR_enriched_motifs, merged_eac_endo_data_matrix_log2, eac_tcga_expression, normal_endo_expression)
EAC_ce_hypoDMR_enriched_motifs_and_sig_loss_exp=data.frame(qvals=EAC_ce_hypoDMR_enriched_motifs_qvals, direction=EAC_ce_hypoDMR_enriched_motifs_direction)
EAC_ce_hypoDMR_enriched_motifs_and_sig_loss_exp=rownames(EAC_ce_hypoDMR_enriched_motifs_and_sig_loss_exp[EAC_ce_hypoDMR_enriched_motifs_and_sig_loss_exp$qvals<0.05 & EAC_ce_hypoDMR_enriched_motifs_and_sig_loss_exp$direction==0,])
data=read.table("homerMotifs.enriched_and_up_regulated.motifs.locations_EAC_hypo.txt", header=T, sep="\t") #See key (EAC_cehypoDMR_enriched_motifs_key.txt) for identifying columns to include from the homer output locations file.
data=process_data(data)
cluster_distance=30
row_groupings=plot_data_clusters(data, cluster_distance, "cornflowerblue") # From this output (row_groupings), I could then identify clusters of TFs. I called a cluster any cluster that had at least 2 TFs in it. The locations of these clusters were used as input to GREAT. See manuscript for more details.

##### Figure 6B #####
GBM_ce_hypoDMR_enriched_motifs=c("BARX1", "CRX", "CUX2", "DMRT1", "DMRTB1", "ELF3", "ELF5", "EOMES", "ESRRA", "FOXA1", "FOXA2", "FOXA3", "FOXM1", "GFI1B", "GRHL2",
                                 "GSC", "HNF1A", "HNF1B", "KLF1", "KLF3", "KLF4", "KLF5", "LHX1", "LHX2", "LHX3", "NANOG", "NFIX", "NKX2-1", "NKX2-2", "NKX2-5",
                                 "NKX3-1", "NKX3-2", "NR1H2", "ONECUT1", "PAX5", "PAX6", "PAX7", "PAX8", "PBX1", "PBX3", "PKNOX1", "POU3F1", "POU3F3", "POU5F1", "SIX1",
                                 "SIX2", "SLC22A2", "SNAI2", "SOX10", "SOX15", "SOX17", "SOX2", "SOX3", "SOX4", "SOX6", "SOX9", "TBX21", "TBX5", "TCF12", "TCF3",
                                 "TGIF2", "THRA", "TP63", "ZEB1", "ZNF264", "ZNF416") # Extract the genes of interest #Using the homer key (http://homer.ucsd.edu/homer/motif/HomerMotifDB/homerResults.html), I identified the gene names for the enriched motifs (and used alternate names as indicated in "GBM_cehypoDMR_homer_results_using_GBM_cehyperDMRs_as_background_filtering_at_q_0.05.txt").
GBM_ce_hypoDMR_enriched_motifs_qvals=calculate_expression_change_pval(GBM_ce_hypoDMR_enriched_motifs, merged_gbm_brain_data_matrix_log2, gbm_tcga_expression, normal_brain_expression)
GBM_ce_hypoDMR_enriched_motifs_direction=calculate_expression_change_direction(GBM_ce_hypoDMR_enriched_motifs, merged_gbm_brain_data_matrix_log2, gbm_tcga_expression, normal_brain_expression)
GBM_ce_hypoDMR_enriched_motifs_and_sig_loss_exp=data.frame(qvals=GBM_ce_hypoDMR_enriched_motifs_qvals, direction=GBM_ce_hypoDMR_enriched_motifs_direction)
GBM_ce_hypoDMR_enriched_motifs_and_sig_loss_exp=rownames(GBM_ce_hypoDMR_enriched_motifs_and_sig_loss_exp[GBM_ce_hypoDMR_enriched_motifs_and_sig_loss_exp$qvals<0.05 & GBM_ce_hypoDMR_enriched_motifs_and_sig_loss_exp$direction==0,])
data=read.table("homerMotifs.enriched_and_up_regulated.motifs.locations_GBM_hypo.txt", header=T, sep="\t") #See key (GBM_cehypoDMR_enriched_motifs_key.txt) for identifying columns to include from the homer output locations file.
data=process_data(data)
cluster_distance=22
row_groupings=plot_data_clusters(data, cluster_distance, "cornflowerblue") # From this output (row_groupings), I could then identify clusters of TFs. I called a cluster any cluster that had at least 2 TFs in it. The locations of these clusters were used as input to GREAT. See manuscript for more details.
