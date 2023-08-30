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


build_meta_data <- function(meta_data_file, exclud_fail_samples = T) {
	removed_columns <- c("Sample_id", "SampleName", "Status");

	## constructing meta_data data frame
	temp_meta_data <- read.table(meta_data_file, header=T, sep="\t", check.names=F, stringsAsFactors=F);
	tumor_rows <- which(temp_meta_data[, "Status"] == "Tumor");
	normal_rows <- which(temp_meta_data[, "Status"] == "Normal");

	sample_names <- temp_meta_data[tumor_rows, "SampleName"];
	status_col <- which(colnames(temp_meta_data) == "Status");

	removed_column_indices <- which(colnames(temp_meta_data) %in% removed_columns);
	meta_data <- temp_meta_data[tumor_rows, -removed_column_indices];
	rownames(meta_data) <- sample_names;
	meta_data[, "TumorID"] <- temp_meta_data[tumor_rows, "Sample_id"];
	meta_data[, "NormalID"] <- temp_meta_data[normal_rows, "Sample_id"];
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