# Jennifer Flynn jaflynn@wustl.edu
# June 16, 2020
# This script reads in all 80 chromHMM files and creates a file that lists the locations of any regions annotated as promoter (states 1, 2, 3, and 4) in at least one of the files.

# Load in the appropriate module(s):
module load bedtools/2.27.1

# Determine regions annotated as "promoters" in any tissue by the chromHMM18 state model
chromHMM18_annotation_files="/chromHMM18_state_files/*" # Please see README.txt for information on where to download these files. They are not stored in this repository due to their size.
cat $chromHMM18_annotation_files | awk '{if ($4=="1_TssA" || $4=="2_TssFlnk" || $4=="3_TssFlnkU" || $4=="4_TssFlnkD") print $0}' | sort -k 1,1 -k2,2n | bedtools merge -i - > "epigenomic_promoter_sites_sorted_merged.bed"
