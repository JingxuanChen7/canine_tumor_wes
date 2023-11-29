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


library("ComplexHeatmap");
library(data.table)
library(dplyr)
library(tibble)
library("Polychrome")

# WES(WGS)
# In the germline VAF file, VAF values for low coverage bases are expected to have NA values
# Varaints with enough coverage that are not called by GATK HaplotypeCaller or fail VariantFilteration should have a VAF value of 0

# These parameters are related to the structure of the input file containing VAF values (don't change any of them)
residue_column_count <- 7; # number of columns describing each variant
meta_row_count <- 1; # number of rows dedicated to meta data (sample ids, and others if applicable) in the VAF input file
seperator <- "/"

# modify directory for I/O files and scripts
file_base_dir <- "/scratch/jc33471/canine_tumor/breed_prediction"
script_dir <- "/home/jc33471/canine_tumor_wes/scripts/breed_prediction"
#file_base_dir <- "/Users/jingxuan/Downloads/tmp/breed"
#script_dir <- "/Users/jingxuan/GitHub/canine_tumor_wes/scripts/breed_prediction"
############ Script customization parameters ########################
# You may modify these parameters as desired
non_na_percentage_cutoff <- 0.8; # all samples must have known VAF values in at least 80% of the breed-specific variants
global_sufficient_cov_cutoff <- 0.8; # all variants must have sufficient coverage in at least 80% of the samples in the discovery dataset

# only use discovery breed
examined_breeds <- c("Shih Tzu", "Schnauzer","Golden Retriever", "Rottweiler", "Greyhound", "Maltese","Yorkshire Terrier","Boxer","Poodle","Cocker Spaniel","Labrador Retriever", "Boston Terrier");
breed_pallete <- c("lightblue",'blue',"#009EFF", "purple", "gray", "yellow","red","#964B00", "orange","#E58FAC","#c2b280","green","black");
breed_order <- 1:13; # This will define the order of which heatmap breed color legends will be displayed

examined_breeds <- c("Shih Tzu", "Schnauzer","Golden Retriever", "Rottweiler", "Greyhound", "Maltese","Yorkshire Terrier","Boxer","Poodle","Cocker Spaniel","Labrador Retriever", "Boston Terrier",
                     'Dachshund','Appenzeller Sennenhund','Collie','Ibizan Hound',
                     'Saint Bernard','German Spitz','Japanese Spitz','Keeshond',
                     'Leonberger','Small Swiss Hound','Swiss Hound',
                     'Petit Basset Griffon Vendeen','Pyrenean Shepherd','Small Munsterlander',
                     'White Swiss Shepherd Dog','Bernese Mountain Dog','Bouvier des Flandres',
                     'English Toy Terrier','Greater Swiss Mountain Dog','Sealyham Terrier',
                     'Smooth Fox Terrier','Toy Fox Terrier','Vizsla');
breed_pallete_more <- glasbey.colors(23)
breed_pallete <- c("lightblue",'blue',"#009EFF", "purple", "gray", "yellow","red","#964B00", "orange","#E58FAC","#c2b280","green",breed_pallete_more,"black");
breed_order <- 1:36; # This will define the order of which heatmap breed color legends will be displayed

cancer_types <- c("MT", "OM", "HSA","BCL","TCL","UCL", "OSA", "GLM","UC","HS","BT","NonT");
disease_order <- c(1:12); # This will define the order of which heatmap disease color legends will be displayed
# See Glasbey palette "Polychrome: Creating and Assessing Qualitative Palettes With Many Colors" (https://www.biorxiv.org/content/10.1101/303883v1.full)
cancer_pallete <- c("blue", "#FF00B6", "#009EFF", "#6f4e37", "#9A4D43", "black", "#FFD300","red","orange","lightblue","green","white");
# See Kelly palette "Polychrome: Creating and Assessing Qualitative Palettes With Many Colors" (https://www.biorxiv.org/content/10.1101/303883v1.full)


############ code dependency paths ########################
# Code to build sample meta data
build_meta_data_code_path <- paste(script_dir, "build_sample_meta_data.R", sep=seperator);

############ Input and output paths ########################
# Please modify these file paths as needed

output_base <- paste(file_base_dir,"output_include_WGS", sep=seperator);

# Input file containing VAF values for all samples for each germline variant: samples as columns and variants as rows
VAF_input_file <- paste(file_base_dir,"germline_VAF_matrix.reset_low_coverage.txt.gz", sep=seperator);
# Input file containing all breed-specific variants
specific_variants_file <- paste(output_base, "all_breed_specific_variants.txt", sep=seperator);
# Input file containing all samples meta data
meta_data_file <- paste(file_base_dir,"breed_prediction_metadata.txt", sep=seperator);

output_png1 <- paste(output_base, "breeds_heatmap_main_beforeQC.png", sep=seperator); # this heatmap won't contain samples with unknown breeds
output_png2 <- paste(output_base, "breeds_heatmap_assignment_beforeQC.png", sep=seperator); # this heatmap will contain samples with unknown breeds
output_png <- c(output_png1, output_png2);
# output_clusters files will have the list of samples ordered as they appear in the heatmaps (used for supplementary tables and breed validation/prediction results)
output_clusters <- c(paste(output_base, "main_clusters.txt", sep=seperator), 
                     paste(output_base, "assignment_clusters.txt", sep=seperator));

################ Main code ############################

VAF_data <- fread(VAF_input_file, header=F, sep="\t")
VAF_data <- setDF(VAF_data)
#, check.names=F, stringsAsFactors=F);
specific_variants_data <- read.table(specific_variants_file, header=T, sep="\t", check.names=F, stringsAsFactors=F);
# find duplicated variants
overlaps <- specific_variants_data[duplicated(specific_variants_data[,c(2,3)],fromLast = T) | duplicated(specific_variants_data[,c(2,3)],fromLast = F),]
# get rid of duplicated variants
specific_variants_data <- specific_variants_data[!(rownames(specific_variants_data) %in% rownames(overlaps)),]

### building meta_data data frame
source(build_meta_data_code_path);
meta_data <- build_meta_data(meta_data_file, exclud_fail_samples = T, include_unpaired = T);
meta_data <- add_tumor_normal_columns(meta_data, unlist(VAF_data[meta_row_count, ]));
# clean breed column
meta_data <- mutate(meta_data, Breed = case_when(Breed == "No breed provided" ~ NA, TRUE ~ Breed))

#which(is.na(meta_data$Breed_Predicted_Results))
#is.na(meta_data$Breed_Predicted_Results)
#meta_data <- meta_data[!is.na(meta_data$Breed_Predicted_Results),]

##### add QC info to metadata for plotting
# meta_data <- dplyr::select(meta_data, -c("Breed","DiseaseAcronym"))
# more_meta_data <- read.table(paste(file_base_dir,"assignment_clusters_meta.txt", sep=seperator),header=T, sep="\t", check.names=F, stringsAsFactors=F) %>% select(c("SampleName","Breed","DiseaseAcronym"))
# meta_data <- dplyr::inner_join(tibble::rownames_to_column(meta_data), more_meta_data, by = c("rowname" = "SampleName"))
# meta_data <- tibble::column_to_rownames(meta_data, var = "rowname")

# now make variant names
variant_names <- apply(VAF_data[-c(1:meta_row_count), c(2:5)], MARGIN=1, function(x) {paste(as.vector(unlist(x)), collapse="_")});
names(variant_names) <- NULL;

# now reading VAF values for normal samples
#normal_VAF_data <- as.matrix(VAF_data[-c(1:meta_row_count), as.vector(meta_data$NormalCol)]); # this line has been modified because data.matrix convert factor to interval codes instead of VAF values (float)
normal_VAF_data <- apply(as.matrix(VAF_data[-c(1:meta_row_count), meta_data[, "NormalCol"]]), 2,as.numeric);
colnames(normal_VAF_data) <- rownames(meta_data);
rownames(normal_VAF_data) <- variant_names;


# now getting the heatmap samples
heatmap_breed_samples <- rownames(meta_data)[which(meta_data[, "Breed"] %in% examined_breeds)];
na_breed_samples <- rownames(meta_data)[which(is.na(meta_data[, "Breed"]) == TRUE)];

heatmap_samples_1 <- c(heatmap_breed_samples); # Samples for first heatmap (no unknown breeds)
heatmap_samples_2 <- c(heatmap_samples_1, na_breed_samples); # samples for second heatmap (with unknown breeds)
sample_list <- list(heatmap_samples_1, heatmap_samples_2);

# new samples
common_breeds <- c("Shih Tzu", "Schnauzer","Golden Retriever", "Rottweiler", "Greyhound", "Maltese","Yorkshire Terrier","Boxer","Poodle","Cocker Spaniel","Labrador Retriever", "Boston Terrier");
common_breed_samples <- rownames(meta_data)[which(meta_data[, "Breed"] %in% common_breeds)];
uncommon_samples <- meta_data[heatmap_breed_samples[!heatmap_breed_samples %in% common_breed_samples],]

# now converting breed-specific variants to variant names for the heatmaps
text_tokens_to_variant_name <- function(text_tokens) {
  gene <- text_tokens[1];
  locus_split <- strsplit(text_tokens[2], fixed=TRUE, split=":")[[1]];
  chromosome <- locus_split[1];
  position <- locus_split[2];
  allele_split <- strsplit(text_tokens[3], fixed=TRUE, split=">")[[1]];
  ref_allele <- allele_split[1];
  alt_allele <- allele_split[2];
  protein <- text_tokens[4];
  return(paste( chromosome, position, ref_allele, alt_allele, sep="_"));
}
heatmap_variants <- apply(specific_variants_data, 1, function(x) {text_tokens_to_variant_name(x)});
heatmap_variants <- intersect(heatmap_variants, rownames(normal_VAF_data));
tmp_heatmap_variants <- heatmap_variants
# important for merging dataset! check high qual variants
non_na_count_cutoff <- ncol(normal_VAF_data) * global_sufficient_cov_cutoff;
for(variant in heatmap_variants) {
  variant_VAF <- normal_VAF_data[variant,];
  non_na_samples <- names(which(is.na(variant_VAF) == FALSE));
  if(length(non_na_samples) < non_na_count_cutoff) {
    tmp_heatmap_variants <- tmp_heatmap_variants[tmp_heatmap_variants != variant]
  }
}
heatmap_variants <- tmp_heatmap_variants
  
  

set.seed(8888);
heatmap_data_list <- list();
for(heatmap_version in c(1,2)) {
  heatmap_samples <- sample_list[[heatmap_version]];
  # building the heatmap data and removing bad samples
  heatmap_data <- normal_VAF_data[heatmap_variants, heatmap_samples];
  
  # check high qual samples
  VAF_non_na_counts <- apply(heatmap_data, 2, function(x) {length(which(is.na(x) == FALSE))});
  sample_count_cutoff <- length(heatmap_variants) * non_na_percentage_cutoff;
  bad_sample_indices <- which(VAF_non_na_counts < sample_count_cutoff);
  if(length(bad_sample_indices) > 0) {
    heatmap_data <- heatmap_data[, -bad_sample_indices];
  }
  heatmap_data_list[[heatmap_version]] <- heatmap_data;
}

# rm(normal_VAF_data, VAF_data);

# This variable is for debugging only (it will store the sample order for each heatmap)
backup_samples <- list();


for(heatmap_version in c(1,2)) {
  heatmap_data <- heatmap_data_list[[heatmap_version]];
  heatmap_samples <- colnames(heatmap_data);
  
  # randomly assigning VAF values to low coverage samples (NA)
  for(variant in rownames(heatmap_data)) {
    na_columns <- which(is.na(heatmap_data[variant,]) == TRUE);
    if(length(na_columns) > 0) {
      known_variant_mut_rates <- heatmap_data[variant, -na_columns];
      random_assignments <- sample(known_variant_mut_rates, length(na_columns), replace=TRUE);
      heatmap_data[variant, na_columns] <- random_assignments;
    }
  }
  
  # defining heatmap annotations and legends
  heatmap_breeds <- c(examined_breeds, "Unknown/Missing")
  #"Unknown");
  breed_colors <- breed_pallete[c(1:length(examined_breeds), length(breed_pallete))]; # Unknown is always assigned last "black" color
  names(breed_colors) <- heatmap_breeds;
  breed_colors <- breed_colors[breed_order];
  breed_info <- meta_data[heatmap_samples, "Breed"];
  breed_info <- factor(breed_info,levels = heatmap_breeds)
  
  disease_colors <- cancer_pallete[disease_order];
  names(disease_colors) <- cancer_types;
  disease_info <- meta_data[heatmap_samples, "DiseaseAcronym"];
  disease_info <- factor(disease_info,levels = cancer_types)
  
  
  
  if(heatmap_version == 1) {
    # do nothing
  } else {
    breed_info[which(heatmap_samples %in% na_breed_samples)] <- "Unknown/Missing";
  }
  
  annotation_legend_param <- list(labels_gp=gpar(fontsize=20), title_gp=gpar(fontsize=20, fontface="plain"), 
                                  grid_height=unit(6, "mm"));
  disease_legend_param <- annotation_legend_param;
  breed_legend_param <- annotation_legend_param;
  disease_legend_param[["ncol"]] <- 2;
  breed_legend_param[["ncol"]] <- 2;
  
  top_annotation <- HeatmapAnnotation(simple_anno_size = unit(1.2,"cm"),
                                      annotation_name_gp = gpar(fontsize=20),
                                      col=list(`Provided breeds` = breed_colors, 
                                               Disease = disease_colors), 
                                      `Provided breeds`=breed_info,
                                      Disease=disease_info,
                                      annotation_legend_param=list(`Provided breeds`=breed_legend_param, 
                                                                   Disease=disease_legend_param));
  
  heatmap_object <- Heatmap(name="ht", heatmap_data, col=colorRampPalette(c("white","blue"))(256),
                            top_annotation=top_annotation, 
                            show_row_names=FALSE, 
                            show_column_names=FALSE, 
                            show_heatmap_legend=FALSE, 
                            show_row_dend=FALSE, 
                            column_dend_height=unit(22, "mm"),
                            use_raster = F); # previous 22 mm
  
  png(output_png[heatmap_version], height=12.5, width=11.5, units="in", res=500);
  heatmap_object <- draw(heatmap_object, annotation_legend_side = "bottom",
                         heatmap_legend_side = "bottom",
                         merge_legend = TRUE);
  decorate_heatmap_body("ht", {grid.rect(gp = gpar(fill = "transparent", col = "black"))})
  dev.off();
  
  ordered_samples <- heatmap_samples[column_order(heatmap_object)];
  backup_samples[[heatmap_version]] <- ordered_samples;
  temp_dataset <- data.frame("SampleName" = ordered_samples,
                             "Provided_Breeds"= meta_data[ordered_samples, "Breed"],
                             #"BreedCluster"= meta_data[ordered_samples, "BreedCluster"],
                             #"FinalBreed" = meta_data[ordered_samples, "FinalBreed"],
                             "BreedQC" = meta_data[ordered_samples, "The_reason_to_exclude"],
                             "DiseaseAcronym"= meta_data[ordered_samples,"DiseaseAcronym"]
                             #"Result"= meta_data[ordered_samples, "Breed_Predicted_Results"],
                             #"DataType"= meta_data[ordered_samples, "DataType"]
  );
  write.table(temp_dataset, file=output_clusters[heatmap_version], quote=FALSE, sep="\t", row.names = FALSE, col.names=TRUE);
  
}

