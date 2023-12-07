library(dplyr)
library(ggtree)
library(ape)
library(treeio)
library(tidytree)
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")

project_dir <- "/home/jc33471/canine_tumor_wes"
## Sample accessino to case id
temp_meta_data <- read.csv("/home/jc33471/canine_tumor_wes/metadata/this_wunpaired_meta.csv", stringsAsFactors = F) #%>% select(c(Sample_ID,Case_ID))

# check paired samples
table_SampleName <- table(temp_meta_data$Case_ID)
paired <- temp_meta_data[temp_meta_data$Case_ID %in% names(table_SampleName)[table_SampleName>1],]
normal_rows <- which(paired[, "Status"] == "Normal");
meta_data <- paired[normal_rows, c("Sample_ID","Case_ID")];
# check un paired samples
unpaired <- temp_meta_data[temp_meta_data$Case_ID %in% names(table_SampleName)[table_SampleName==1],]
unpaired_rows <- which(unpaired[, "Status"] == "Tumor" | unpaired[, "Status"] == "Normal")
meta_data_unpaired <- unpaired[unpaired_rows, c("Sample_ID","Case_ID")];
meta_data <- rbind(meta_data, meta_data_unpaired)

#### breed predcition 
df_assignment <- read.table("/home/jc33471/canine_tumor_wes/results/breed_prediction/assignment_clusters.txt", header = T, sep = "\t",stringsAsFactors = F)
df_assignment <- left_join(df_assignment, meta_data, by = c("SampleName" = "Case_ID"))
# the breed cluster are assigned based on breeds_heatmap_assignment_beforeQC.png
df_validation <- df_assignment %>%
  mutate(BreedCluster = case_when(
    as.numeric(rownames(.)) >= 1 & as.numeric(rownames(.)) <= 33 ~ "Rottweiler",
    as.numeric(rownames(.)) >= 34 & as.numeric(rownames(.)) <= 215 ~ "Golden Retriever",
    as.numeric(rownames(.)) >= 216 & as.numeric(rownames(.)) <= 244 ~ "Shih Tzu",
    as.numeric(rownames(.)) >= 245 & as.numeric(rownames(.)) <= 271 ~ "Greyhound",
    as.numeric(rownames(.)) >= 272 & as.numeric(rownames(.)) <= 321 ~ "Labrador Retriever",
    as.numeric(rownames(.)) >= 322 & as.numeric(rownames(.)) <= 340 ~ "Schnauzer",
    as.numeric(rownames(.)) >= 341 & as.numeric(rownames(.)) <= 413 ~ "Maltese",
    as.numeric(rownames(.)) >= 414 & as.numeric(rownames(.)) <= 423 ~ "Boston Terrier",
    as.numeric(rownames(.)) >= 424 & as.numeric(rownames(.)) <= 439 ~ "Yorkshire Terrier",
    as.numeric(rownames(.)) >= 441 & as.numeric(rownames(.)) <= 493 ~ "Poodle",
    as.numeric(rownames(.)) >= 610 & as.numeric(rownames(.)) <= 650 ~ "Boxer",
    as.numeric(rownames(.)) >= 651 & as.numeric(rownames(.)) <= 680 ~ "Cocker Spaniel",
    TRUE ~ "Unknown/Missing"
  ))

#### phylo breed validation
out_tree <- "/home/jc33471/canine_tumor_wes/results/phylogenetics/breedPlusMissingSample_breedSpecific.nwk"
all_tr <- read.tree(out_tree) 
in_id <- "/home/jc33471/canine_tumor_wes/results/phylogenetics/breedPlusMissingSample_breedSpecific_merged_germline_variants_SNP.mdist.id"
# read id file because of phylip only allow sample name <10 chars
all_id <- read.table(in_id) %>% mutate(id = row_number())
# rename taxa
all_tr <- rename_taxa(all_tr, all_id, 3, 1)
#write.tree(all_tr,"/home/jc33471/canine_tumor_wes/results/phylogenetics/breedPlusMissingSample_breedSpecific_named.nwk")


p <- ggtree(all_tr)


# the node for pure breeds are inferred from all_breedPlusMissingSample_breedSpecific.pdf
ShihTzu_node <- getMRCA(all_tr, tip = c("SRR7780994","SRR7780861")); ShihTzu_node_taxa <- get_taxa_name(p, node = ShihTzu_node)

Schnauzer_node <- getMRCA(all_tr, tip = c("SRR9911360","SRR9911377")); Schnauzer_node_taxa <- get_taxa_name(p, node = Schnauzer_node)

GoldenRetriever_node <- getMRCA(all_tr, tip = c("SAMN03436940","SAMN03436734")); GoldenRetriever_node_taxa <- get_taxa_name(p, node = GoldenRetriever_node)

Rottweiler_node <- getMRCA(all_tr, tip = c("SAMN03436933","SAMN03436267")); Rottweiler_node_taxa <- get_taxa_name(p, node = Rottweiler_node)

Greyhound_node <- getMRCA(all_tr, tip = c("SAMN03732738","SAMN03732751")); Greyhound_node_taxa <- get_taxa_name(p, node = Greyhound_node)

Maltese_node <- getMRCA(all_tr, tip = c("SRR9911350","SRR7780771")); Maltese_node_taxa <- get_taxa_name(p, node = Maltese_node)

YorkshireTerrier_node <- getMRCA(all_tr, tip = c("SRR7781069","SRR13571597")); YorkshireTerrier_node_taxa <- get_taxa_name(p, node = YorkshireTerrier_node)

Boxer_node <- getMRCA(all_tr, tip = c("SRR10351705","SRR10351788")); Boxer_node_taxa <- get_taxa_name(p, node = Boxer_node)

Poodle_node <- getMRCA(all_tr, tip = c("SRR17139270","SRR7781108")); Poodle_node_taxa <- get_taxa_name(p, node = Poodle_node)

CockerSpaniel_node <- getMRCA(all_tr, tip = c("SAMN03436931","SAMN03436476")); CockerSpaniel_node_taxa <- get_taxa_name(p, node = CockerSpaniel_node)

LabradorRetriever_node <- getMRCA(all_tr, tip = c("SRR17139222","SRR17139233")); LabradorRetriever_node_taxa <- get_taxa_name(p, node = LabradorRetriever_node)

BostonTerrier_node <- getMRCA(all_tr, tip = c("SRR10351748","SRR10351800")); BostonTerrier_node_taxa <- get_taxa_name(p, node = BostonTerrier_node)


phylo_assign <- data.frame(Sample_ID = all_tr$tip.label) %>%
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
                                Sample_ID %in% LabradorRetriever_node_taxa ~ "Labrador Retriever",
                                Sample_ID %in% BostonTerrier_node_taxa ~ "Boston Terrier",
                                TRUE ~ "Unknown/Missing"))

# join VAF cluster and phylogenetic results
df_joined <- left_join(df_validation, phylo_assign, by = "Sample_ID") %>%
  mutate(BreedFinal = case_when(BreedCluster == BreedPhylo ~ BreedCluster,
                                is.na(BreedPhylo) ~ BreedCluster,
                                BreedCluster == "Unknown/Missing" & BreedPhylo != "Unknown/Missing" & !is.na(BreedPhylo) ~ BreedPhylo,
                                BreedCluster != "Unknown/Missing" & BreedPhylo == "Unknown/Missing"  ~ BreedCluster,
                                BreedCluster != "Unknown/Missing" & BreedPhylo != "Unknown/Missing" & !is.na(BreedPhylo) & BreedCluster != BreedPhylo ~ "Conflict",
                                TRUE ~ NA)) %>%
  mutate(Provided_Breeds = case_when(is.na(Provided_Breeds) ~ "Unknown/Missing",
                                     TRUE ~ Provided_Breeds)) %>%
  mutate(BreedPhylo = case_when(is.na(BreedPhylo) ~ "Tumor_only",
                                     TRUE ~ BreedPhylo))

# print out the inconsistency
inconsistency <- df_joined %>%
  filter(BreedCluster != BreedPhylo | is.na(BreedPhylo))

# number of missing breeds
table(df_joined$BreedCluster, useNA = "ifany")
table(df_joined$BreedPhylo, useNA = "ifany")
table(df_joined$BreedFinal, useNA = "ifany")
table(df_joined$Provided_Breeds, useNA = "ifany")

count(df_joined, BreedCluster, BreedPhylo, BreedFinal) 
# output table with assignment results and phylogenetic assignment results
write.csv(phylo_assign, file = "/scratch/jc33471/canine_tumor_test/breed_prediction/assignment_clusters_phylogenetics.csv", quote = T, row.names = F)

inconsistency <- phylo_assign %>%
  filter(BreedCluster != BreedPhylo)
# output table for inconsistent clusters
write.csv(inconsistency, file = "/scratch/jc33471/canine_tumor_test/breed_prediction/assignment_inconsistency.csv", quote = T, row.names = F)


