library(reshape2)
library(dplyr)
library(readxl)
# install.packages("readxl")

meta_dir <- "~/Github/canine_tumor_wes/metadata"

dog10K <- read_excel(paste0(meta_dir, "/WGS/Dog10K_supp_table.xlsx"), skip = 2) %>%
  filter(Category == "Breed_Dogs") %>%
  select(c("Sample Name","Effective Autosomal Mean Coverage","Breed/Type")) %>%
  rename("Sample_id" = "Sample Name", "Coverage" = "Effective Autosomal Mean Coverage", "Breed" = "Breed/Type") %>%
  mutate(SampleName = Sample_id, DiseaseAcronym = "NA", Status = "Normal", The_reason_to_exclude = "Pass QC")

write.csv(dog10K, file = paste0(meta_dir, "/WGS/Dog10K_breeds.csv"), quote = F, row.names = F)


nih <- read_excel(paste0(meta_dir, "/WGS/NIH_supp.xlsx"), skip = 1) %>%
  #filter(BioProject == "PRJNA448733") %>%
  select(c("Name_ID_SRA","CoverageAll","Phylo_Results")) %>%
  rename("Case_ID" = "Name_ID_SRA", "Coverage" = "CoverageAll", "Breed" = "Phylo_Results")

dbvdc <- read_excel(paste0(meta_dir, "/WGS/DBVDC_supp.xlsx"), skip = 1) %>%
  #filter(`Study accession` == "PRJEB16012") %>%
  select(c("Sample ID","Coverage depth","Breed")) %>%
  rename("Case_ID" = "Sample ID", "Coverage" = "Coverage depth")

all_wgs <- rbind(dog10K,nih,dbvdc)
hist(all_wgs$Coverage, breaks = 50)

# breed frequency
dog10K_freq <- table(dog10K$Breed)[order(table(dog10K$Breed), decreasing = T)]
nih_freq <- table(nih$Breed)[order(table(nih$Breed), decreasing = T)]
