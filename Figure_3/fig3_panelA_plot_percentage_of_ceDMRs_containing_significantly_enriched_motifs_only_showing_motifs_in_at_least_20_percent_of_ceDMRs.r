# Jennifer Flynn jaflynn@wustl.edu
# March 19, 2020
# This script plots a heatmap showing the percentage of ceDMRs containing each enriched motif (only plotting enriched motifs contained within at least 20% of ceDMR groups). Motif enrichment was calculated using Homer v4.9 and opposite DMR groups as background. See manuscript for further details.

# Read in the enriched TF motif information
EAC_hyperDMR_homer_output=read.table("EAC_hyperDMR_homer_knownResults.txt", sep = "\t")
EAC_hypoDMR_homer_output=read.table("EAC_hypoDMR_homer_knownResults.txt", sep = "\t")
GBM_hyperDMR_homer_output=read.table("GBM_hyperDMR_homer_knownResults.txt", sep = "\t")
GBM_hypoDMR_homer_output=read.table("GBM_hypoDMR_homer_knownResults.txt", sep = "\t")

colnames(EAC_hyperDMR_homer_output)=c("Motif Name",	"Consensus",	"P-value", "Log P-value",	"q-value (Benjamini) EAC hyperDMRs",	"# of Target Sequences with Motif",	"% of Target Sequences with Motif EAC hyperDMRs",	"# of Background Sequences with Motif",	"% of Background Sequences with Motif")
colnames(EAC_hypoDMR_homer_output)=c("Motif Name",	"Consensus",	"P-value", "Log P-value",	"q-value (Benjamini) EAC hypoDMRs",	"# of Target Sequences with Motif",	"% of Target Sequences with Motif EAC hypoDMRs",	"# of Background Sequences with Motif",	"% of Background Sequences with Motif")
colnames(GBM_hyperDMR_homer_output)=c("Motif Name",	"Consensus",	"P-value", "Log P-value",	"q-value (Benjamini) GBM hyperDMRs",	"# of Target Sequences with Motif",	"% of Target Sequences with Motif GBM hyperDMRs",	"# of Background Sequences with Motif",	"% of Background Sequences with Motif")
colnames(GBM_hypoDMR_homer_output)=c("Motif Name",	"Consensus",	"P-value", "Log P-value",	"q-value (Benjamini) GBM hypoDMRs",	"# of Target Sequences with Motif",	"% of Target Sequences with Motif GBM hypoDMRs",	"# of Background Sequences with Motif",	"% of Background Sequences with Motif")

EAC_hyperDMR_homer_output$Motif_and_name=paste0(EAC_hyperDMR_homer_output$`Motif Name`,"_", EAC_hyperDMR_homer_output$Consensus)
EAC_hypoDMR_homer_output$Motif_and_name=paste0(EAC_hypoDMR_homer_output$`Motif Name`,"_", EAC_hypoDMR_homer_output$Consensus)
GBM_hyperDMR_homer_output$Motif_and_name=paste0(GBM_hyperDMR_homer_output$`Motif Name`,"_", GBM_hyperDMR_homer_output$Consensus)
GBM_hypoDMR_homer_output$Motif_and_name=paste0(GBM_hypoDMR_homer_output$`Motif Name`,"_", GBM_hypoDMR_homer_output$Consensus)


# Combine the data
merged_homer_output=merge(EAC_hyperDMR_homer_output, EAC_hypoDMR_homer_output, by="Motif_and_name")
merged_homer_output=merge(merged_homer_output, GBM_hyperDMR_homer_output, by="Motif_and_name")
merged_homer_output=merge(merged_homer_output, GBM_hypoDMR_homer_output, by = "Motif_and_name")
homer_percentages_and_qvals_full=data.frame(motif=merged_homer_output$`Motif Name.x`, 
                                            EAC_hyperDMRs_percent=merged_homer_output$`% of Target Sequences with Motif EAC hyperDMRs`, 
                                            EAC_hypoDMRs_percent=merged_homer_output$`% of Target Sequences with Motif EAC hypoDMRs`,
                                            GBM_hyperDMRs_percent=merged_homer_output$`% of Target Sequences with Motif GBM hyperDMRs`,
                                            GBM_hypoDMRs_percent=merged_homer_output$`% of Target Sequences with Motif GBM hypoDMRs`,
                                            EAC_hyperDMRs_qvals=merged_homer_output$`q-value (Benjamini) EAC hyperDMRs`,
                                            EAC_hypoDMRs_qvals=merged_homer_output$`q-value (Benjamini) EAC hypoDMRs`,
                                            GBM_hyperDMRs_qvals=merged_homer_output$`q-value (Benjamini) GBM hyperDMRs`,
                                            GBM_hypoDMRs_qvals=merged_homer_output$`q-value (Benjamini) GBM hypoDMRs`)
# If the qval is greater than 0.01 or if the percentage of target sites with the motif is < 20%, then set their percentage to NA for plotting.
homer_percentages_and_qvals_full$EAC_hyperDMRs_percent[homer_percentages_and_qvals_full$EAC_hyperDMRs_qvals>0.01]=0
homer_percentages_and_qvals_full$EAC_hypoDMRs_percent[homer_percentages_and_qvals_full$EAC_hypoDMRs_qvals>0.01]=0
homer_percentages_and_qvals_full$GBM_hyperDMRs_percent[homer_percentages_and_qvals_full$GBM_hyperDMRs_qvals>0.01]=0
homer_percentages_and_qvals_full$GBM_hypoDMRs_percent[homer_percentages_and_qvals_full$GBM_hypoDMRs_qvals>0.01]=0
homer_percentages_and_qvals_full$EAC_hyperDMRs_qvals=NULL
homer_percentages_and_qvals_full$EAC_hypoDMRs_qvals=NULL
homer_percentages_and_qvals_full$GBM_hyperDMRs_qvals=NULL
homer_percentages_and_qvals_full$GBM_hypoDMRs_qvals=NULL

homer_percentages_and_qvals_full$EAC_hyperDMRs_percent=as.numeric(gsub(pattern = "%",replacement = "", x = homer_percentages_and_qvals_full$EAC_hyperDMRs_percent))
homer_percentages_and_qvals_full$EAC_hypoDMRs_percent=as.numeric(gsub(pattern = "%",replacement = "", x = homer_percentages_and_qvals_full$EAC_hypoDMRs_percent))
homer_percentages_and_qvals_full$GBM_hyperDMRs_percent=as.numeric(gsub(pattern = "%",replacement = "", x = homer_percentages_and_qvals_full$GBM_hyperDMRs_percent))
homer_percentages_and_qvals_full$GBM_hypoDMRs_percent=as.numeric(gsub(pattern = "%",replacement = "", x = homer_percentages_and_qvals_full$GBM_hypoDMRs_percent))

homer_percentages_and_qvals_full$EAC_hyperDMRs_percent[homer_percentages_and_qvals_full$EAC_hyperDMRs_percent<20]=0
homer_percentages_and_qvals_full$EAC_hypoDMRs_percent[homer_percentages_and_qvals_full$EAC_hypoDMRs_percent<20]=0
homer_percentages_and_qvals_full$GBM_hyperDMRs_percent[homer_percentages_and_qvals_full$GBM_hyperDMRs_percent<20]=0
homer_percentages_and_qvals_full$GBM_hypoDMRs_percent[homer_percentages_and_qvals_full$GBM_hypoDMRs_percent<20]=0

homer_percentages=homer_percentages_and_qvals_full[rowSums(homer_percentages_and_qvals_full[, -1], na.rm=TRUE)>0, ]
homer_percentages[is.na(homer_percentages)]=0
row.names(homer_percentages)=homer_percentages$motif
homer_percentages$motif=NULL
homer_percentages_matrix=as.matrix(homer_percentages)

# Create a heatmap
my_palette <- colorRampPalette(c("#DBF6FF", "#57ABE8", "#FCF9AF", "#F5714A"))(n = 30)
heatmap.2(homer_percentages_matrix, trace="none",col=my_palette)
dev.off()
