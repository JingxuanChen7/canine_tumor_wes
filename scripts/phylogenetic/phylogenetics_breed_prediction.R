library(dplyr)
library(ggtree)

vaf_assign <- read.table("/scratch/jc33471/canine_tumor_test/breed_prediction/assignment_clusters.txt", header = T,sep = "\t",stringsAsFactors = F) %>%
  filter(DataType != "WGS")
meta <- read.csv("/home/jc33471/canine_tumor_wes/metadata/data_collection_old.csv", stringsAsFactors = F) %>% filter(Status=="Normal") %>%
  select(c(Sample_ID,Case_ID))
vaf_assign <- left_join(vaf_assign, meta, by = c("SampleName" = "Case_ID"))

prefix <- "all_breedPlusMissingSample_breedSpecific"
out_dir <- "/scratch/jc33471/canine_tumor_test/breed_prediction/"
out_tree <- paste0(out_dir,"/",prefix,".nwk")

all_tr <- read.tree(out_tree) 
p <- ggtree(all_tr)

# the node for pure breeds are inferred from all_breedPlusMissingSample_breedSpecific.pdf
ShihTzu_node <- getMRCA(all_tr, tip = c("SRR7780865","SRR7780861")); ShihTzu_node_taxa <- get_taxa_name(p, node = ShihTzu_node)

Schnauzer_node <- getMRCA(all_tr, tip = c("ERR1672168-ERR1672255","SRR7781005")); Schnauzer_node_taxa <- get_taxa_name(p, node = Schnauzer_node)

GoldenRetriever_node <- getMRCA(all_tr, tip = c("SAMN03436829","SAMN03436711")); GoldenRetriever_node_taxa <- get_taxa_name(p, node = GoldenRetriever_node)

Rottweiler_node <- getMRCA(all_tr, tip = c("SAMN03436933","SAMN03436267")); Rottweiler_node_taxa <- get_taxa_name(p, node = Rottweiler_node)

Greyhound_node <- getMRCA(all_tr, tip = c("SAMN03436724","SAMN03732750")); Greyhound_node_taxa <- get_taxa_name(p, node = Greyhound_node)

Maltese_node <- getMRCA(all_tr, tip = c("SRR7780806","SRR7780948")); Maltese_node_taxa <- get_taxa_name(p, node = Maltese_node)

YorkshireTerrier_node <- getMRCA(all_tr, tip = c("SRR7781069","SRR7781079")); YorkshireTerrier_node_taxa <- get_taxa_name(p, node = YorkshireTerrier_node)

Boxer_node <- getMRCA(all_tr, tip = c("SAMN03436531","ERR1672159-ERR1672246")); Boxer_node_taxa <- get_taxa_name(p, node = Boxer_node)

Poodle_node <- getMRCA(all_tr, tip = c("SAMN03436608","SAMN03436573")); Poodle_node_taxa <- get_taxa_name(p, node = Poodle_node)

CockerSpaniel_node <- getMRCA(all_tr, tip = c("SAMN03436280","SAMN03436476")); CockerSpaniel_node_taxa <- get_taxa_name(p, node = CockerSpaniel_node)



phylo_assign <- vaf_assign %>%
  mutate(BreedPhylo = case_when(Sample_ID %in% ShihTzu_node_taxa ~ "Shih Tzu",
                                Sample_ID %in% Schnauzer_node_taxa ~ "Schnauzer",
                                Sample_ID %in% GoldenRetriever_node_taxa ~ "Golden Retriever",
                                Sample_ID %in% Rottweiler_node_taxa ~ "Rottweiler",
                                Sample_ID %in% Greyhound_node_taxa ~ "Greyhound",
                                Sample_ID %in% Maltese_node_taxa ~ "Maltese",
                                Sample_ID %in% YorkshireTerrier_node_taxa ~ "Yorkshire Terrier",
                                Sample_ID %in% Boxer_node_taxa ~ "Boxer",
                                Sample_ID %in% Poodle_node_taxa ~ "Poodle",
                                Sample_ID %in% CockerSpaniel_node_taxa ~ "Cocker Spaniel",
                                TRUE ~ "Unknown/Missing"))
