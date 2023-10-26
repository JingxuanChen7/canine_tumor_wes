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


build_meta_data <- function(meta_data_file, exclud_fail_samples = T, include_unpaired = T) {
	removed_columns <- c("Sample_id", "SampleName", "Status");

	## constructing meta_data data frame
	temp_meta_data <- read.table(meta_data_file, header=T, sep="\t", check.names=F, stringsAsFactors=F);
	# check paired samples
	table_SampleName <- table(temp_meta_data$SampleName)
	paired <- temp_meta_data[temp_meta_data$SampleName %in% names(table_SampleName)[table_SampleName>1],]
	
	# parse only paired samples
	tumor_rows <- which(paired[, "Status"] == "Tumor");
	normal_rows <- which(paired[, "Status"] == "Normal");
	sample_names <- paired[tumor_rows, "SampleName"];
	status_col <- which(colnames(paired) == "Status");
	removed_column_indices <- which(colnames(paired) %in% removed_columns);
	# if paired, keep normal QC results
	meta_data <- paired[normal_rows, -removed_column_indices];
	rownames(meta_data) <- sample_names;
	meta_data[, "TumorID"] <- paired[tumor_rows, "Sample_id"];
	meta_data[, "NormalID"] <- paired[normal_rows, "Sample_id"];
	
	if (include_unpaired) {
	  unpaired <- temp_meta_data[temp_meta_data$SampleName %in% names(table_SampleName)[table_SampleName==1],]
	  unpaired_rows <- which(unpaired[, "Status"] == "Tumor" | unpaired[, "Status"] == "Normal")
	  sample_names_unpaired <- unpaired[unpaired_rows, "SampleName"];
	  meta_data_unpaired <- unpaired[unpaired_rows, -removed_column_indices];
	  rownames(meta_data_unpaired) <- sample_names_unpaired;
	  # treat all unpaired as "normal" to fit afterwards codes
	  meta_data_unpaired[, "TumorID"] <- unpaired[unpaired_rows, "Sample_id"];
	  meta_data_unpaired[, "NormalID"] <- unpaired[unpaired_rows, "Sample_id"];
	  meta_data <- rbind(meta_data, meta_data_unpaired)
	}

	if (exclud_fail_samples){
		meta_data <- meta_data[meta_data$The_reason_to_exclude =="Pass QC",]
	}
	return(meta_data);
}

add_tumor_normal_columns <- function(input_meta_data, column_names) {
	output_meta_data <- input_meta_data;
	output_meta_data[, "NormalCol"] <- NA;
	output_meta_data[, "TumorCol"] <- NA;
	for(i in 1:nrow(output_meta_data)) {
		tumor_col <- which(column_names == output_meta_data[i, "TumorID"]);
		normal_col <- which(column_names == output_meta_data[i, "NormalID"]);
		if(length(normal_col) > 1) {
			normal_col <- normal_col[1];
		}
		output_meta_data[i, "NormalCol"] <- normal_col;
		output_meta_data[i, "TumorCol"] <- tumor_col;
	}
	return(output_meta_data);
}