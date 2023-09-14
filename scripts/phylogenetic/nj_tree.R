#!/usr/bin/env Rscript
library(ape)
library(parallel)
library(ggtree)
library(dplyr)

# parse arguments from command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript --vanilla nj_tree.R <thread> <in_vcf> <prefix> <tree_dir> <out_tree>", call.=FALSE)
} else if (length(args) == 5) {
    thread <- args[1]
    in_seq <- args[2]
    prefix <- args[3]
    tree_dir <- args[4]
    out_tree <- args[5]
}
thread <- 10
in_seq <- "/scratch/jc33471/canine_tumor_test/breed_prediction/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy_breed_sample.min4.fasta"
out_tree <- "/scratch/jc33471/canine_tumor_test/breed_prediction/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy_breed_sample.nwk"
out_fig <- "/scratch/jc33471/canine_tumor_test/breed_prediction/PanCancer_57WGS_disc_val_sep_germline_VAF_0119.reset_low_coverage_copy_breed_sample.pdf"
metadata <- "/home/jc33471/canine_tumor_wes/metadata/data_collection_old.csv"

# bootstrap nj tree
all_align <- read.dna(in_seq, format = "fasta", as.character = F, as.matrix = T)
all_tr <- nj(dist.dna(all_align, model="raw", pairwise.deletion = T)) # best tree
all_boot <- boot.phylo(all_tr, all_align, function(x) nj(dist.dna(x, model="raw", pairwise.deletion = T)), B = 100, mc.cores = as.numeric(thread), trees = T)
all_tr$node.label <- all_boot$BP
write.tree(all_tr, file = out_tree, append = FALSE)

# parse metadata
out_tree <- "/scratch/jc33471/canine_tumor_test/breed_prediction/breed_specific_breed_sample.nwk"
out_fig <- "/scratch/jc33471/canine_tumor_test/breed_prediction/breed_specific_breed_sample.pdf"

all_tr <- read.tree(out_tree)
meta <- read.csv(metadata, stringsAsFactors = F) %>% filter(Status=="Normal") %>%
  select(c(Sample_ID,Case_ID,DiseaseAcronym2,Breed_info)) %>%
  mutate(taxa = paste(Case_ID,DiseaseAcronym2,Breed_info,sep = "_"))
p <- ggtree(all_tr) %<+% meta +
  geom_tiplab( aes(label = taxa), size=1.5, align = F ) +
  geom_tippoint(aes(fill = Breed_info), shape = 21, size = 1.5)

ggsave(plot = p, out_fig, device = "pdf", width = 12, height = 23, units = "in" )
