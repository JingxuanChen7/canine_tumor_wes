#######################################################################################
###                                                                                 ###
###    Copyright (C) 2020,  Burair Alsaihati                                        ###
###                                                                                 ###
###    This program is free software: you can redistribute it and/or modify         ###
###    it under the terms of the GNU General Public License as published by         ###
###    the Free Software Foundation version 3                                       ###
###                                                                                 ###
###    This program is distributed in the hope that it will be useful,              ###
###    but WITHOUT ANY WARRANTY; without even the implied warranty of               ###
###    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                ###
###    GNU General Public License for more details.                                 ###
###                                                                                 ###
###    You should have received a copy of the GNU General Public License            ###
###    along with this program.  If not, see <https://www.gnu.org/licenses/>.       ###
###                                                                                 ###
###    Email: burair.alsaihati25@uga.edu, burair_99@yahoo.com, szhao@uga.edu        ###
###                                                                                 ###
#######################################################################################
library(data.table)
library(R.utils)

################# Parameters related to the VAF file structure ######################## 
# Please don't change any of these parameters
residue_column_count <- 7;
meta_row_count <- 1;
################# End of parameters related to the VAF file structure ######################## 

depth_cutoff <- 10; # Cutoff for minimum coverage. Any variant with coverage < 10 will be assigned a VAF value of NA.
seperator <- "/"

#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Burair_pan_scripts/breed_prediction_test"
base_dir <- "/scratch/jc33471/canine_tumor_test/breed_prediction"

################# Input files ######################## 
# make sure to modify the paths to the correct ones
VAF_file <- paste(base_dir, "PanCancer_disc_val_merged_germline_VAF_01_01_2021.txt.gz", sep=seperator);
depth_file <- paste(base_dir, "PanCancer_disc_val_merged_germline_depths_01_01_2021.txt.gz", sep=seperator);

################# Output files ########################
# make sure to modify the paths to the correct ones
VAF_output_file <- paste(base_dir, "germline_VAF_matrix.reset_low_coverage.txt.gz", sep=seperator);


################ Main code ############################

VAF_data <- fread(VAF_file, header=F, sep="\t")
                  #check.names=F, stringsAsFactors=F);
VAF_data <- setDF(VAF_data)

depth_data <- fread(depth_file, header=F, sep="\t")
                    #check.names=F, stringsAsFactors=F);
depth_data <- setDF(depth_data)


# filtering VAF based on BED depth
for(i in (residue_column_count+1):ncol(depth_data)) {
  depths_as_strings <- as.vector(unlist(depth_data[-c(1:meta_row_count), i]));
  bed_depths_as_numbers <- as.numeric(sapply(depths_as_strings, function(x) {strsplit(x, ",", fixed=T)[[1]][1]},  USE.NAMES=F));
  vcf_depths_as_numbers <- as.numeric(sapply(depths_as_strings, function(x) {strsplit(x, ",", fixed=T)[[1]][2]},  USE.NAMES=F));
  na_depth <- which(is.na(bed_depths_as_numbers));
  vcf_na_depth <- which(is.na(vcf_depths_as_numbers));
  low_depth <- which(bed_depths_as_numbers < depth_cutoff);
  reset_rows <- union(c(na_depth, low_depth), vcf_na_depth) + meta_row_count;
  VAF_data[reset_rows, i] <- "NA";
  print(paste("col", (i-residue_column_count), "of", (ncol(depth_data)-residue_column_count)));
  flush.console();
}

# removing duplicated variants if any
# now make variant names
variant_names <- apply(VAF_data[-c(1:meta_row_count), c(1:5)], MARGIN=1, function(x) {paste(as.vector(unlist(x)), collapse="_")});
names(variant_names) <- NULL;
duplicated_indices <- which(duplicated(variant_names) == TRUE) + meta_row_count;
length(duplicated_indices);
if(length(duplicated_indices) > 0) {
  VAF_data <- VAF_data[-duplicated_indices,];
}

gz_file <- gzfile(VAF_output_file, "w");


write.table(VAF_data, file=gz_file, sep="\t", quote=F, row.names=F, col.names=F);
close(gz_file)

#fwrite(VAF_data,file=VAF_output_file, sep="\t", quote=F, row.names=F, col.names=F,compress ="gzip")


