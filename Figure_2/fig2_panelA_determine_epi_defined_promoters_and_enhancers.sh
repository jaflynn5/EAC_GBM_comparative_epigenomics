# Jennifer Flynn jaflynn@wustl.edu
# June 16, 2020
# This script reads in all 80 chromHMM files and creates the following files:
# 1): a file that lists the locations of any regions annotated as promoter (states 1, 2, 3, and 4) in at least one of the chromHMM18 state files
# 2): a file that lists the locations of any regions annotated as active enhancers (state 9 and 10) in at least one of the chromHMM18 state files
# 3): a file that lists the locations of any regions called a promoter by either the genomic definions or chromHMM18 state files
# 4): a file that lists the locations of any regions annotated as promoter (epigenetically and/or genetically) that do not overlap active enhancer annotations in any chromHMM18 state file
# 5): a file that lists the locations of any regions annoated as an active enhancer in any chromHMM18 state file that does not overlap a promoter annoatation (epigenetically and/or genetically)
# 6): a file that lists the locations of any regions annotated as an active enhancer in any chromHMM18 state file that is also annotated as a promoter (epigenetically and/or genetically)

# Load in the appropriate module(s):
module load bedtools/2.27.1

# Read in the path to the chromHMM18 state files
chromHMM18_annotation_files="/chromHMM18_annotation_files/*" #Please see the README.txt for information on where to download these files (as they are too large to include here)

# Determine regions annotated as "promoters" in any tissue by the chromHMM18 state model
cat $chromHMM18_annotation_files | awk '{if ($4=="1_TssA" || $4=="2_TssFlnk" || $4=="3_TssFlnkU" || $4=="4_TssFlnkD") print $0}' | sort -k 1,1 -k2,2n | bedtools merge -i - > "epigenomic_promoter_sites_sorted_merged.bed"

# Determine regions annotated as "active enhancers" in any tissue by the chromHMM18 state model
cat $chromHMM18_annotation_files | awk '{if ($4=="9_EnhA1" || $4=="10_EnhA2") print $0}' | sort -k 1,1 -k2,2n  > "active_enhancer_sites_sorted.bed"
bedtools merge -i "active_enhancer_sites_sorted.bed" > "active_enhancer_sites_sorted_merged.bed"

# Determine regions annotated as "promoters" in any tissue by the chromHMM18 state models and/or by genomic annotation
genomic_promoters="refGene_Promoters2500bp_sorted_merged.bed"
cat $genomic_promoters "epigenomic_promoter_sites_sorted_merged.bed" | sort -k 1,1 -k2,2n | bedtools merge -i - > "union_of_genomic_and_epigenomic_promoters.bed"

# Determine which regions are annotated as promoters that do not overlap active enhancers
bedtools subtract -a "union_of_genomic_and_epigenomic_promoters.bed" -b "active_enhancer_sites_sorted_merged.bed" > "promoters_not_overlapping_active_enhancers.bed"

# Determine which regions are annotated as active enhancers that do not overlap promoters
bedtools subtract -a "active_enhancer_sites_sorted_merged.bed" -b "union_of_genomic_and_epigenomic_promoters.bed" > "active_enhancers_not_overlapping_promoters.bed"

# Determine which regions are annotated as promoters and active enhancers
bedtools intersect -a "union_of_genomic_and_epigenomic_promoters.bed" -b "active_enhancer_sites_sorted_merged.bed" > "promoters_and_active_enhancers.bed"
