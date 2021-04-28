# Jennifer Flynn jaflynn@wustl.edu
# January 15, 2020
# This script reads in the list of GBM hyperDMRs (500bp regions), determines if each hyperDMR overlaps the list of enhancer regions for each tissue type (adult brain and fetal tissue types: E072, E067, E073, E068, E069, E071, E074, E053, E054, E082, E070, E125, E081), and excludes GBM hyperDMRs that overlap to refGene promoters

# Load libraries
library(bedr)
library(ggplot2)

# Read in the list of GBM hyperDMRs
gbm_hyperDMRs=read.table("GBM_hyperDMRs.bed")
colnames(gbm_hyperDMRs)=c("chr", "start", "end")
gbm_hyperDMRs$chr=as.character(gbm_hyperDMRs$chr)

# Read in the promoter locations
promoters=read.table("refGene_Promoters2500bp_sorted_merged.bed")
colnames(promoters)=c("chr", "start", "end")
promoters$chr=as.character(promoters$chr)

# Read in the list of chromHMM 18-state files for fetal and adult brain tissues (files can be downloaded from: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/indivModels/default_init/)
fetal_and_adult_brain_chromHMM_18_state_files=c("E072_15_coreMarks_dense.bed",
                                                "E067_15_coreMarks_dense.bed",
                                                "E073_15_coreMarks_dense.bed",
                                                "E068_15_coreMarks_dense.bed",
                                                "E069_15_coreMarks_dense.bed",
                                                "E071_15_coreMarks_dense.bed",
                                                "E074_15_coreMarks_dense.bed",
                                                "E053_15_coreMarks_dense.bed",
                                                "E054_15_coreMarks_dense.bed",
                                                "E082_15_coreMarks_dense.bed",
                                                "E070_15_coreMarks_dense.bed",
                                                "E125_15_coreMarks_dense.bed",
                                                "E081_15_coreMarks_dense.bed")

# Determine which GBM hyperDMRs overlap regions of the genome that contain enhancers (only state 7_Enh) in the tissues listed above that are outside of promoters
for (tissue_file in fetal_and_adult_brain_chromHMM_18_state_files)
{
  sample=strsplit(basename(tissue_file), split = "_")[[1]][1]
  data=read.table(tissue_file, sep = "\t", skip=1, header = FALSE)
  data=data[data$V4=="7_Enh",1:3]
  colnames(data)=c("chr", "start", "end")
  data$chr=as.character(data$chr)
  gbm_hyperDMRs_within_enhancers=bedr(input = list(a = gbm_hyperDMRs, b = data), method="intersect", params="-c")
  gbm_hyperDMRs_within_enhancers=gbm_hyperDMRs_within_enhancers[gbm_hyperDMRs_within_enhancers$V4>0,]
  gbm_hyperDMRs_within_enhancers_outside_promoters=bedr(input = list(a = gbm_hyperDMRs_within_enhancers, b = promoters[is.valid.region(promoters),]), method="subtract", params="-A")
  gbm_hyperDMRs_within_enhancers_outside_promoters=gbm_hyperDMRs_within_enhancers_outside_promoters[gbm_hyperDMRs_within_enhancers_outside_promoters$V4>0,]
  if (exists("full_data"))
  {
    temp_data=data.frame(DMR=rownames(gbm_hyperDMRs_within_enhancers_outside_promoters), value=1)
    colnames(temp_data)=c("DMR", sample)
    full_data=merge(full_data, temp_data, by="DMR", all=TRUE)
  } else {
    full_data=data.frame(DMR=rownames(gbm_hyperDMRs_within_enhancers_outside_promoters), value=1)
    colnames(full_data)=c("DMR", sample)
  }
}
full_data[is.na(full_data)] <- 0


#### Figure 5A ####
# Create heatmap
dataset_for_plotting=full_data[,2:ncol(full_data)]
col_order <- c("E125", "E082", "E081", "E074", "E073", "E072", "E071", "E070", "E069", "E068", "E067", "E053", "E054")
dataset_for_plotting <- dataset_for_plotting[, col_order]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E054),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E053),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E067),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E068),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E069),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E070),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E071),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E072),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E073),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E074),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E081),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E082),]
dataset_for_plotting <- dataset_for_plotting[order(-dataset_for_plotting$E125),]
dataset_for_plotting_matrix=as.matrix(dataset_for_plotting)
my_palette <- colorRampPalette(c("gray", "blue"))(n = 2)
heatmap.2(dataset_for_plotting_matrix, trace="none", dendrogram="column", col=my_palette)
dev.off()

# Create a barplot showing the percentage of GBM hyperDMRs that overlapped the enhancer state (state 7_Enh) in each cell type, but not to RefGene 2.5kb promoters
total_GBM_hyperDMRs=dim(gbm_hyperDMRs)[1]
percentage_DMRs_overlapping_enh_7_not_promoters=data.frame(sample=c("E081", "E125", "E070", "E082", "E054", "E053", "E074", "E073", "E067", "E068", "E072", "E071", "E069"),
                                                           percentage=c(sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E081"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E125"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E070"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E082"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E054"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E053"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E074"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E073"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E067"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E068"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E072"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E071"])/total_GBM_hyperDMRs,
                                                                        sum(dataset_for_plotting_matrix[,colnames(dataset_for_plotting_matrix)=="E069"])/total_GBM_hyperDMRs))
percentage_DMRs_overlapping_enh_7_not_promoters$sample=as.factor(percentage_DMRs_overlapping_enh_7_not_promoters$sample)
percentage_DMRs_overlapping_enh_7_not_promoters$sample <- factor(percentage_DMRs_overlapping_enh_7_not_promoters$sample, levels = c("E081", "E125", "E070", "E082", "E054", "E053", "E074", "E073", "E067", "E068", "E072", "E071", "E069"))
ggplot(data=percentage_DMRs_overlapping_enh_7_not_promoters, aes(x=sample, y=percentage)) + geom_bar(stat="identity") + coord_flip()
dev.off()


#### Figure 5B ####
# Create boxplot
percentage_DMRs_overlapping_enh_7_not_promoters$type=c(rep("fetal", 6), rep("adult", 7))
ggplot(percentage_DMRs_overlapping_enh_7_not_promoters, aes(x=type, y=percentage)) + 
  geom_boxplot()
t.test(percentage_DMRs_overlapping_enh_7_not_promoters$percentage[percentage_DMRs_overlapping_enh_7_not_promoters$type=="adult"],
       percentage_DMRs_overlapping_enh_7_not_promoters$percentage[percentage_DMRs_overlapping_enh_7_not_promoters$type=="fetal"])


#### Figure 5C ####
# Determine which GBM ce-hyperDMRs overlap with H3K4me1 and H3K27ac peaks in different adult and fetal brain samples
adult_and_fetal_brain_enh_potential_GBM_hyperDMRs=full_data$DMR
adult_and_fetal_brain_enh_potential_GBM_hyperDMRs_list=unlist(as.list(adult_and_fetal_brain_enh_potential_GBM_hyperDMRs))

# Determine which adult and fetal brain enhancer potential GBM hyperDMRs overlap with adult and/or fetal H3K4me1 peaks (H3K4me1 peak files can be downloaded from: http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/)
h3k4me1_peak_files=c("E053-H3K4me1.narrowPeak_sorted.bed",
                     "E054-H3K4me1.narrowPeak_sorted.bed",
                     "E067-H3K4me1.narrowPeak_sorted.bed",
                     "E068-H3K4me1.narrowPeak_sorted.bed",
                     "E069-H3K4me1.narrowPeak_sorted.bed",
                     "E070-H3K4me1.narrowPeak_sorted.bed",
                     "E071-H3K4me1.narrowPeak_sorted.bed",
                     "E072-H3K4me1.narrowPeak_sorted.bed",
                     "E073-H3K4me1.narrowPeak_sorted.bed",
                     "E074-H3K4me1.narrowPeak_sorted.bed",
                     "E081-H3K4me1.narrowPeak_sorted.bed",
                     "E082-H3K4me1.narrowPeak_sorted.bed",
                     "E125-H3K4me1.narrowPeak_sorted.bed")

for (peak_file in h3k4me1_peak_files)
{
  sample=paste0("h3k4me1_",strsplit(basename(peak_file), split = "-")[[1]][1])
  data=read.table(peak_file, sep = "\t", header = FALSE)
  colnames(data)=c("chr", "start", "end")
  data$chr=as.character(data$chr)
  gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks=bedr(input = list(a = bedr.sort.region(as.character(adult_and_fetal_brain_enh_potential_GBM_hyperDMRs_list)), b = data), method="intersect", params="-c")
  gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks=gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks[gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks$V4>0,]
  
  if (exists("full_peak_data"))
  {
    temp_data=data.frame(DMR=gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks$index, value=1)
    colnames(temp_data)=c("DMR", sample)
    full_peak_data=merge(full_peak_data, temp_data, by="DMR", all=TRUE)
  } else {
    full_peak_data=data.frame(DMR=gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks$index, value=1)
    colnames(full_peak_data)=c("DMR", sample)
  }
}

# Read in the H3K27ac peak files for the fetal and adult brain tissues (for which it exists). H3K27ac peak files can be downloaded from: http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/.
h3k27ac_peak_files=c("E067-H3K27ac.narrowPeak_sorted",
                     "E068-H3K27ac.narrowPeak_sorted",
                     "E069-H3K27ac.narrowPeak_sorted",
                     "E071-H3K27ac.narrowPeak_orted",
                     "E072-H3K27ac.narrowPeak_sorted",
                     "E073-H3K27ac.narrowPeak_sorted",
                     "E074-H3K27ac.narrowPeak_sorted",
                     "E125-H3K27ac.narrowPeak_sorted")

# Determine which adult and fetal brain enhancer potential GBM hyperDMRs overlap with adult and/or fetal H3K4me1 and/or H3K27ac peaks
for (peak_file in h3k27ac_peak_files)
{
  sample=paste0("h3k27ac_",strsplit(basename(peak_file), split = "-")[[1]][1])
  data=read.table(peak_file, sep = "\t", header = FALSE)
  colnames(data)=c("chr", "start", "end")
  data$chr=as.character(data$chr)
  gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks_and_or_h3k27ac_peaks=bedr(input = list(a = bedr.sort.region(as.character(adult_and_fetal_brain_enh_potential_GBM_hyperDMRs_list)), b = data), method="intersect", params="-c")
  gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks_and_or_h3k27ac_peaks=gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks_and_or_h3k27ac_peaks[gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks_and_or_h3k27ac_peaks$V4>0,]
  temp_data=data.frame(DMR=gbm_hyperDMRs_within_enhancers_overlapping_h3k4me1_peaks_and_or_h3k27ac_peaks$index, value=1)
  colnames(temp_data)=c("DMR", sample)
  full_peak_data=merge(full_peak_data, temp_data, by="DMR", all=TRUE)
}

full_peak_data[is.na(full_peak_data)] <- 0
full_peak_data=as.matrix(full_peak_data[,2:dim(full_peak_data)[2]])

heatmap.2(full_peak_data, trace="none", col=c("gray", "orange"), Colv = FALSE)
dev.off()
heatmap.2(full_peak_data, trace="none", col=c("gray", "green"), Colv = FALSE)
dev.off()

# Create a barplot showing the percentage of GBM ce-hyperDMRs that overlapped the enhancer state (state 7_Enh) in each cell type, but not to RefGene 2.5kb promoters, and that also overlapped H3K27ac and H3K4me1 peaks
total_GBM_ce_hyperDMRs=dim(dataset_for_plotting_matrix)[1]
percentage_DMRs_overlapping_enh_7_and_histone_mods_not_promoters=data.frame(sample=c("E067_H3K27ac", "E068_H3K27ac", "E069_H3K27ac", "E071_H3K27ac", "E072_H3K27ac",
                                                                                     "E073_H3K27ac", "E074_H3K27ac", "E125_H3K27ac", "E067_H3K4me1", "E068_H3K4me1",
                                                                                     "E069_H3K4me1", "E071_H3K4me1", "E072_H3K4me1", "E073_H3K4me1", "E074_H3K4me1",
                                                                                     "E053_H3K4me1", "E054_H3K4me1", "E070_H3K4me1", "E081_H3K4me1", "E082_H3K4me1",
                                                                                     "E125_H3K4me1"),
                                                                            percentage=c(sum(full_peak_data[,colnames(full_peak_data)=="h3k27ac_E067"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k27ac_E068"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k27ac_E069"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k27ac_E071"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k27ac_E072"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k27ac_E073"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k27ac_E074"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k27ac_E125"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E067"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E068"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E069"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E071"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E072"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E073"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E074"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E053"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E054"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E070"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E081"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E082"])/total_GBM_ce_hyperDMRs,
                                                                                         sum(full_peak_data[,colnames(full_peak_data)=="h3k4me1_E125"])/total_GBM_ce_hyperDMRs))
percentage_DMRs_overlapping_enh_7_and_histone_mods_not_promoters$sample=as.factor(percentage_DMRs_overlapping_enh_7_and_histone_mods_not_promoters$sample)
percentage_DMRs_overlapping_enh_7_and_histone_mods_not_promoters$sample <- factor(percentage_DMRs_overlapping_enh_7_and_histone_mods_not_promoters$sample, 
                                                                                  levels = c("E067_H3K27ac", "E068_H3K27ac", "E069_H3K27ac", "E071_H3K27ac",
                                                                                             "E072_H3K27ac", "E073_H3K27ac", "E074_H3K27ac", "E125_H3K27ac",
                                                                                             "E067_H3K4me1", "E068_H3K4me1", "E069_H3K4me1", "E071_H3K4me1",
                                                                                             "E072_H3K4me1", "E073_H3K4me1", "E074_H3K4me1", "E053_H3K4me1",
                                                                                             "E054_H3K4me1", "E070_H3K4me1", "E081_H3K4me1", "E082_H3K4me1",
                                                                                             "E125_H3K4me1"))

ggplot(data=percentage_DMRs_overlapping_enh_7_and_histone_mods_not_promoters, aes(x=sample, y=percentage)) + geom_bar(stat="identity") + coord_flip()
dev.off()


#### Figure 5D ####
# Create boxplot H3K27ac
percentage_DMRs_overlapping_enh_7_and_h3k27ac_peaks_not_promoters=percentage_DMRs_overlapping_enh_7_and_histone_mods_not_promoters[1:8,]
percentage_DMRs_overlapping_enh_7_and_h3k27ac_peaks_not_promoters$type=c(rep("adult", 7), "fetal")
ggplot(percentage_DMRs_overlapping_enh_7_and_h3k27ac_peaks_not_promoters, aes(x=type, y=percentage)) + 
  geom_boxplot()


#### Figure 5E ####
# Create boxplot H3K4me1
percentage_DMRs_overlapping_enh_7_and_h3K4me1_peaks_not_promoters=percentage_DMRs_overlapping_enh_7_and_histone_mods_not_promoters[9:dim(percentage_DMRs_overlapping_enh_7_and_histone_mods_not_promoters)[1],]
percentage_DMRs_overlapping_enh_7_and_h3K4me1_peaks_not_promoters$type=c(rep("adult", 7), rep("fetal", 6))
ggplot(percentage_DMRs_overlapping_enh_7_and_h3K4me1_peaks_not_promoters, aes(x=type, y=percentage)) + 
  geom_boxplot()
t.test(percentage_DMRs_overlapping_enh_7_and_h3K4me1_peaks_not_promoters$percentage[percentage_DMRs_overlapping_enh_7_and_h3K4me1_peaks_not_promoters$type=="adult"],
       percentage_DMRs_overlapping_enh_7_and_h3K4me1_peaks_not_promoters$percentage[percentage_DMRs_overlapping_enh_7_and_h3K4me1_peaks_not_promoters$type=="fetal"])
