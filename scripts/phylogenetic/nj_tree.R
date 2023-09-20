#!/usr/bin/env Rscript
library(ape)
library(parallel)
library(ggtree)
library(dplyr)
library(ggnewscale)
options(ignore.negative.edge=TRUE)

# parse arguments from command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript --vanilla nj_tree.R <thread> <in_vcf> <prefix> <out_dir> <metadata>", call.=FALSE)
} else if (length(args) == 5) {
    thread <- args[1]
    in_seq <- args[2]
    prefix <- args[3]
    out_dir <- args[4]
    metadata <- args[5]
}
# thread <- 10
# prefix <- "DLA-12_breedSample"
# in_seq <- paste0("/scratch/jc33471/canine_tumor_test/breed_prediction/",prefix,".min4.fasta")
# out_dir <- "/scratch/jc33471/canine_tumor_test/breed_prediction/"
# metadata <- "/home/jc33471/canine_tumor_wes/metadata/data_collection_old.csv"

# output tree and figure files
out_tree <- paste0(out_dir,"/",prefix,".nwk")
out_fig <- paste0(out_dir,"/",prefix,".pdf")


if (! file.exists(out_tree)){
  # bootstrap nj tree
  all_align <- read.dna(in_seq, format = "fasta", as.character = F, as.matrix = T)
  all_tr <- nj(dist.dna(all_align, model="raw", pairwise.deletion = T)) # best tree
  #all_boot <- boot.phylo(all_tr, all_align, function(x) nj(dist.dna(x, model="raw", pairwise.deletion = T)), B = 100, mc.cores = as.numeric(thread), trees = T)
  #all_tr$node.label <- all_boot$BP
  write.tree(all_tr, file = out_tree, append = FALSE)
}

# color palette for breed
setcolors <- c(
  "Shih Tzu" = "lightblue",
  "Schnauzer" = "blue",
  "Golden Retriever" = "#009EFF",
  "Rottweiler" = "purple",
  "Greyhound" = "gray",
  "Maltese" = "yellow",
  "Yorkshire Terrier" = "red",
  "Boxer" = "#964B00",
  "Poodle" = "orange",
  "Cocker Spaniel" = "#E58FAC",
  "No breed provided" = "black"
)

disease_colors <- c(
  "MT" = "blue", 
  "OM" = "#FF00B6", 
  "HSA" = "#009EFF",
  "BCL" = "#6f4e37",
  "TCL" = "#9A4D43",
  "UCL" = "black", 
  "OSA" = "#FFD300", 
  "GLM" = "red"
)


# out_tree <- "/scratch/jc33471/canine_tumor_test/breed_prediction/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy_breed_sample.nwk"
# out_fig <- "/scratch/jc33471/canine_tumor_test/breed_prediction/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy_breed_sample.pdf"
# 
# out_tree <- "/scratch/jc33471/canine_tumor_test/breed_prediction/breed_specific_breed_sample.nwk"
# out_fig <- "/scratch/jc33471/canine_tumor_test/breed_prediction/breed_specific_breed_sample.pdf"

all_tr <- read.tree(out_tree) 

# parse metadata

meta <- read.csv(metadata, stringsAsFactors = F) %>% filter(Status=="Normal") %>%
  select(c(Sample_ID,Case_ID,DiseaseAcronym2,Breed_info)) %>%
  mutate(taxa = paste(Case_ID,DiseaseAcronym2,Breed_info,sep = "_"))

# tree
p <- ggtree(all_tr) %<+% meta +
  geom_tiplab( aes(label = " "), size=1.5, align = T ) +
  geom_tippoint(aes(fill = Breed_info), shape = 21, size = 1.5) +
  #geom_nodelab(aes(x = branch, label = label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 1), size = 2, color = "black") +
  scale_fill_manual(values = setcolors, name = "Breeds")

# data frame for heat map
dftree_heat_breed <- meta; rownames(dftree_heat_breed) <- dftree_heat_breed$Sample_ID; dftree_heat_breed <- select(dftree_heat_breed, Breed_info)
dftree_heat_dis <- meta; rownames(dftree_heat_dis) <- dftree_heat_dis$Sample_ID; dftree_heat_dis <- select(dftree_heat_dis, DiseaseAcronym2)

# add heatmap to tree
p2 <- gheatmap(p, dftree_heat_breed, width=0.1, colnames=F, colnames_position = "top", font.size = 3.5, colnames_offset_y = 0, color = NA) +
  scale_fill_manual(values = c(setcolors), name = "Breed") 

p3 <- p2 + new_scale_fill()

p_final <- gheatmap(p3, dftree_heat_dis, offset = 0.01, width=0.1, colnames=F, colnames_position = "top", font.size = 3.5, colnames_offset_y = 0, color = NA) +
  scale_fill_manual(values = c(disease_colors), name = "Disease") 

# print figure
ggsave(plot = p_final, out_fig, device = "pdf", width = 12, height = 16, units = "in" )
