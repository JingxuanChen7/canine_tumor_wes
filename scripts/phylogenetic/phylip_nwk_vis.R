#!/usr/bin/env Rscript
library(ape)
library(parallel)
library(ggtree)
library(dplyr)
library(ggnewscale)
library(treeio)
options(ignore.negative.edge=TRUE)

# parse arguments from command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript --vanilla phylip_nwk_vis.R <prefix> <in_tree> <in_id> <metadata> <out_fig>", call.=FALSE)
} else if (length(args) == 5) {
  prefix <- args[1]
  in_tree <- args[2]
  in_id <- args[3]
  metadata <- args[4]
  out_fig <- args[5]
}

# prefix <- "breedSample_breedSpecific"
# in_tree <- paste0("/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/",prefix,".nwk")
# in_id <- paste0("/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/",prefix,"_merged_germline_variants_SNP.mdist.id")
# metadata <- "/home/jc33471/canine_tumor_wes/metadata/data_collection_old.csv"
# out_fig <- paste0("/scratch/jc33471/canine_tumor/phylogenetics/merge_vcf/",prefix,".pdf")

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
  "Labrador Retriever" = "#c2b280",
  "Boston Terrier" = "green",
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

# read tree file from phylip
all_tr <- read.tree(in_tree)

# read id file
all_id <- read.table(in_id) %>% mutate(id = row_number())

# rename taxa
all_tr <- rename_taxa(all_tr, all_id, 3, 1)


if (prefix == "breedSample_all" | prefix == "breedPlusMissingSample_all"){
  shih_node <- getMRCA(all_tr, tip = c("SRR7781092","SRR7780790"))
  all_tr <- root(all_tr, node = shih_node)
}


# parse meta for coloring, exclude duplicated sample IDs
meta <- read.csv(metadata, stringsAsFactors = F) %>% filter(Status=="Normal") %>%
  select(c(Sample_ID,Case_ID,DiseaseAcronym2,Breed_info)) %>%
  mutate(taxa = paste(Case_ID,Sample_ID,sep = "_")) %>%
  distinct(Sample_ID, .keep_all = T)



# tree
p <- ggtree(all_tr) %<+% meta +
  geom_tiplab( aes(label = taxa), size=1.5 ) +
  geom_tippoint(aes(fill = Breed_info), shape = 21, size = 1.5) +
  theme_tree2() +
  geom_nodelab(aes(x = branch, label = label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 1), size = 2, color = "black") +
  scale_fill_manual(values = setcolors, name = "Breeds")

# data frame for heat map
dftree_heat_breed <- meta; rownames(dftree_heat_breed) <- dftree_heat_breed$Sample_ID; dftree_heat_breed <- select(dftree_heat_breed, Breed_info)
dftree_heat_dis <- meta; rownames(dftree_heat_dis) <- dftree_heat_dis$Sample_ID; dftree_heat_dis <- select(dftree_heat_dis, DiseaseAcronym2)

# add heatmap to tree
p2 <- gheatmap(p, dftree_heat_breed, width=0.1, colnames=F, colnames_position = "top", font.size = 3.5, colnames_offset_y = 0, color = NA) +
  scale_fill_manual(values = c(setcolors), name = "Breed") 

# p3 <- p2 + new_scale_fill()

# p_final <- gheatmap(p3, dftree_heat_dis, offset = 0.005, width=0.1, colnames=F, colnames_position = "top", font.size = 3.5, colnames_offset_y = 0, color = NA) +
#     scale_fill_manual(values = c(disease_colors), name = "Disease") 

# print figure
ggsave(plot = p2, out_fig, device = "pdf", width = 16, height = 16, units = "in" )
