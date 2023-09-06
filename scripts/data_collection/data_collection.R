library(reshape2)
library(dplyr)
library(R.utils)
library(tidyr)
library(stringr)

# remote path
#meta_dir <- "/home/jc33471/canine_tumor_wes/metadata"
# local path
meta_dir <- "~/Github/canine_tumor_wes/metadata"

# PRJDB16014 NIPS not tumor
Bioproject <- "PRJDB16014"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
# run interactively in linux
# system(paste0("esearch -db sra -query ", Bioproject," | efetch -format runinfo > ", srarunlist))
df_PRJDB16014 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample,Sample_Name,BREED,Age,Platform,ReleaseDate,Center.Name, sex,Tissue) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample_Name, Breed = BREED, Organization = Center.Name, Sex = sex,) %>%
  mutate(Tumor_type = "Non-Tumor (Inflammatory Colorectal Polyp ICRP)", Status = "Normal", Symbol = NA, Sample_status = NA, Disease_stage = NA, Survival_status = NA, Tumor_grade = NA)

# # PRJEB57227 sanger
# Bioproject <- "PRJEB57227"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
# df_PRJEB57227 <- read.csv(srarunlist, stringsAsFactors = F) %>%
#   filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
#   select(Run,BioProject,BioSample, Subject_ID, Platform,ReleaseDate, gender,PHENOTYPE) %>%
#   rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Subject_ID, Sex = gender,Tissue = PHENOTYPE) %>%
#   mutate(Tumor_type = "Urothelial carcinoma (UC)", Organization = "SC", Status = NA, Breed = NA, Symbol = NA, Sample_status = NA, Age = NA, Disease_stage = NA, Survival_status = NA, Tumor_grade = NA) %>%
#   mutate(Age = as.character(Age)) %>%
#   group_by(BioProject,Sample_ID,Case_ID,Status,Breed,Platform,ReleaseDate,Age,Sex,Tissue,Tumor_type,Organization,Sample_status,Symbol,Survival_status,Tumor_grade,Disease_stage) %>% 
#   arrange(Run_ID) %>% summarise(Run_ID = paste(Run_ID, collapse="-")) %>%
#   mutate(Case_ID = gsub("(a|b|c)$","_\\1", Case_ID, perl = T)) %>%
#   separate(Case_ID, c("Case_ID","Status"), sep = "_")

# PRJNA891496 Ghent non-tumor
Bioproject <- "PRJNA891496"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJNA891496 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Sample.Name, BREED, Platform,ReleaseDate, sex,Tissue) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex) %>%
  mutate(Tumor_type = "Non-Tumor (Unclassified)",Organization = "Ghent University", Status = "Normal", Symbol = NA, Sample_status = NA, Age = NA, Disease_stage = NA, Survival_status = NA, Tumor_grade = NA)

# PRJEB53653
Bioproject <- "PRJEB53653"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJEB53653 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "SINGLE" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Sample_Name, Platform, ReleaseDate,tissue_type) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample_Name, Tissue = tissue_type) %>%
  mutate(Tumor_type = "Mammary neoplasia (MT)", Organization = "Uppsala University", Status = NA, Sex = NA, Breed = NA, Symbol = NA, Sample_status = NA, Age = NA, Disease_stage = NA, Survival_status = NA, Tumor_grade = NA) %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Case_ID = gsub("CMT_","",Case_ID)) %>%  mutate(Case_ID = gsub("_DUP_BQSR.bam","",Case_ID)) %>% mutate(Case_ID = gsub("_MERGED_Dup_BQSR.bam","",Case_ID)) %>%
  separate(Case_ID, c("Status","Case_ID"), sep = "_") %>%
  mutate(Status = case_when(Status == "TUMOR" ~ "Tumor", Status == "NORMAL" ~ "Normal")) %>%
  mutate(Tissue = gsub("DNA from ","",Tissue))
  

# PRJDB10211
Bioproject <- "PRJDB10211"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJDB10211 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Library.Name, Platform, ReleaseDate) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Library.Name) %>%
  mutate(Tumor_type = "histiocytic sarcoma", Organization = "The University of Tokyo (UT-VMS)", Status = NA, Sex = NA, Tissue = NA, Breed = NA, Symbol = NA, Sample_status = NA, Age = NA, Disease_stage = NA, Survival_status = NA, Tumor_grade = NA) %>%
  mutate(Age = as.character(Age)) %>%
  separate(Case_ID, c("Case_ID", "Status")) %>%
  mutate(Status = case_when(Status == "T" ~ "Tumor", Status == "N" ~ "Normal"))
  
# PRJNA786469 oral melanoma
Bioproject <- "PRJNA786469"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJNA786469 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Sample.Name, BREED, Platform, ReleaseDate, Age, sex, Tissue, disease_stage, disease) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex, Disease_stage = disease_stage, Tumor_type = disease) %>%
  mutate(Organization = "CNRS - University Rennes", Status = NA, Symbol = NA, Sample_status = NA, Survival_status = NA, Tumor_grade = NA) %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Status = case_when(Tumor_type == "" ~ "Normal", TRUE ~ "Tumor")) %>%
  mutate(Case_ID = gsub("(A|B)$","A-B",Case_ID, perl = T), Case_ID = gsub("(C|D)$","C-D",Case_ID, perl = T), Case_ID = gsub("(E|F)$","E-F",Case_ID, perl = T), Case_ID = gsub("(G|H)$","G-H",Case_ID, perl = T),
         Case_ID = gsub("(I|J)$","I-J",Case_ID, perl = T), Case_ID = gsub("(K|L)$","K-L",Case_ID, perl = T), Case_ID = gsub("(M|N)$","M-N",Case_ID, perl = T), Case_ID = gsub("(O|P)$","O-P",Case_ID, perl = T),
         Case_ID = gsub("(Q|R)$","Q-R",Case_ID, perl = T), Case_ID = gsub("(S|T)$","S-T",Case_ID, perl = T), Case_ID = gsub("(U|V)$","U-V",Case_ID, perl = T), Case_ID = gsub("(W|X)$","W-X",Case_ID, perl = T),
         Case_ID = gsub("(Y|Z)$","Y-Z",Case_ID, perl = T), Case_ID = gsub("(0|1)$","0-1",Case_ID, perl = T), Case_ID = gsub("(2|3)$","2-3",Case_ID, perl = T), Case_ID = gsub("(4|5)$","4-5",Case_ID, perl = T),
         Case_ID = gsub("(6|7)$","6-7",Case_ID, perl = T), Case_ID = gsub("(8|9)$","8-9",Case_ID, perl = T)) %>%
  mutate(Case_ID = gsub("_.*$","",Case_ID, perl = T)) %>%
  group_by(BioProject,Sample_ID,Case_ID,Status,Breed,Platform,ReleaseDate,Age,Sex,Tissue,Tumor_type,Organization,Sample_status,Symbol,Survival_status,Tumor_grade,Disease_stage) %>% 
  arrange(Run_ID) %>% summarise(Run_ID = paste(Run_ID, collapse="-")) %>%
  mutate(Tumor_type = "Oral Melanoma")

  
# PRJNA752630 Padua
Bioproject <- "PRJNA752630"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJNA752630 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Sample.Name, BREED, Platform, ReleaseDate, Age, sex, Tissue, Library.Name) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex, Status = Library.Name) %>%
  mutate(Tumor_type = "Diffuse large B-cell Lymphoma", Organization = "University of Padua", Symbol = NA, Sample_status = NA, Survival_status = NA, Tumor_grade = NA,Disease_stage = NA) %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Status = case_when(grepl("DLBCL", Status) ~ "Tumor", grepl("Punch", Status) ~ "Normal"))

# # PRJEB36323 sanger (target exome sequencing)
# Bioproject <- "PRJEB36323"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
# df_PRJEB36323 <- read.csv(srarunlist, stringsAsFactors = F) %>%
#   filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
#   select(Run,BioProject,BioSample, Subject_ID, Platform,ReleaseDate, gender,PHENOTYPE) %>%
#   rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Subject_ID, Sex = gender,Tissue = PHENOTYPE) %>%
  # mutate(Tumor_type = "Urothelial carcinoma", Organization = "SI", Status = NA, Breed = NA, Symbol = NA, Sample_status = "FFPE", Age = NA, Disease_stage = NA, Survival_status = NA, Tumor_grade = NA)

# PRJNA701141 Padua cell line
Bioproject <- "PRJNA701141"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJNA701141 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Sample.Name, BREED, Platform, ReleaseDate, Age, sex, Tissue, Library.Name) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex, Status = Library.Name) %>%
  mutate(Tumor_type = "osteosarcoma (cOSA)", Organization = "University of Padua", Symbol = NA, Sample_status = "FFPE", Survival_status = NA, Tumor_grade = NA,Disease_stage = NA) %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Status = case_when(grepl("CL", Status) ~ "Cell Line", grepl("Normal", Status) ~ "Normal", TRUE ~ "Tumor")) %>%
  mutate(Case_ID = case_when(grepl("FFPE_2", Case_ID) ~ "FFPE_2", TRUE ~ Case_ID))

# PRJNA695534 Cornell 
Bioproject <- "PRJNA695534"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJNA695534 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Sample.Name, BREED, Platform, ReleaseDate, Age, sex, Tissue, store_cond) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex, Sample_status = store_cond) %>%
  mutate(Tumor_type = "B-Cell Lymphoma (cBCL)", Organization = "Cornell University", Status = NA, Symbol = NA, Survival_status = NA, Tumor_grade = NA,Disease_stage = NA) %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Status = case_when(grepl("normal", Tissue) ~ "Normal", grepl("tumor", Tissue) ~ "Tumor")) %>%
  mutate(Tissue = gsub("matched normal DNA from EDTA ", "", Tissue), Tissue = gsub("tumor DNA from excisional biopsy of ","",Tissue)) %>%
  mutate(Case_ID = gsub("\\..*$","", Case_ID))

# PRJNA680382 Missouri 2 sample normal-primary-Metastatic
Bioproject <- "PRJNA680382"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJNA680382 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Sample.Name, BREED, Platform, ReleaseDate, Age, sex, Tissue, Sample_type) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex, Status = Tissue, Tissue = Sample_type) %>%
  mutate(Tumor_type = "Appendicular Osteosarcoma (OSA)", Organization = "University of Missouri", Sample_status = NA, Symbol = NA, Survival_status = NA, Tumor_grade = NA,Disease_stage = NA) %>%
  mutate(Age = as.character(Age)) %>%
  mutate(Case_ID = gsub("(M|T).*$","",Case_ID,perl = T)) %>%
  mutate(Status = case_when(grepl("tumor",Status) ~ "Tumor", TRUE ~ Status)) %>%
  mutate(Disease_stage = case_when(Tissue == "Lung" ~ "Metastatic", Tissue == "Bone" ~ "Primary") )
add_normal_row <- slice(df_PRJNA680382, rep(c(3,5), each = 1)) %>%
  mutate(Case_ID = case_when(Case_ID == "680" ~ "994", Case_ID == "1032" ~ "1166"))
df_PRJNA680382 <- rbind(df_PRJNA680382, add_normal_row)

# PRJNA677995 Translational Genomics Research Institute
Bioproject <- "PRJNA677995"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJNA677995 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Sample.Name, BREED, Platform, ReleaseDate, Age, sex, Tissue) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex) %>%
  mutate(Tumor_type = "splenic angiosarcoma (hemangiosarcoma)", Organization = "Translational Genomics Research Institute", Status = NA, Sample_status = NA, Symbol = NA, Survival_status = NA, Tumor_grade = NA,Disease_stage = NA) %>%
  mutate(Age = as.character(Age)) %>%
  separate(Case_ID, c("Case_ID", "Status"), sep = "_") %>%
  mutate(Status = case_when(Status == "normal" ~ "Normal", Status == "tumor" ~ "Tumor")) %>%
  group_by(BioProject,Sample_ID,Case_ID,Status,Breed,Platform,ReleaseDate,Age,Sex,Tissue,Tumor_type,Organization,Sample_status,Symbol,Survival_status,Tumor_grade,Disease_stage) %>% 
  arrange(Run_ID) %>% summarise(Run_ID = paste(Run_ID, collapse="-"))

# PRJNA630029 Helsinki not tumor study
Bioproject <- "PRJNA630029"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
df_PRJNA630029 <- read.csv(srarunlist, stringsAsFactors = F) %>%
  filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
  select(Run,BioProject,BioSample, Sample.Name, BREED, Platform, ReleaseDate, Age, sex, Tissue) %>%
  rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex) %>%
  mutate(Tumor_type = "Non-Tumor", Organization = "University of Helsinki", Status = "Normal", Sample_status = NA, Symbol = NA, Survival_status = NA, Tumor_grade = NA,Disease_stage = NA)

# # PRJNA616374 Colorado State, only tumor samples
# Bioproject <- "PRJNA616374"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
# df_PRJNA616374 <- read.csv(srarunlist, stringsAsFactors = F) %>%
#   filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
#   select(Run,BioProject,BioSample, Sample.Name, BREED, Platform, ReleaseDate, Age, sex, Tissue) %>%
#   rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex) %>%
#   mutate(Age = as.character(Age)) %>%
#   mutate(Tumor_type = "Primary bladder tumor", Organization = "Colorado State University", Status = "Tumor", Sample_status = NA, Symbol = NA, Survival_status = NA, Tumor_grade = NA,Disease_stage = NA)
# 
# # PRJNA613479 Colorado State, only tumor samples
# Bioproject <- "PRJNA613479"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
# df_PRJNA613479 <- read.csv(srarunlist, stringsAsFactors = F) %>%
#   filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
#   select(Run,BioProject,BioSample, Sample.Name, BREED, Platform, ReleaseDate, Age, sex, Tissue, disease_stage) %>%
#   rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Sample.Name, Breed = BREED, Sex = sex, Disease_stage = disease_stage) %>%
#   mutate(Age = as.character(Age)) %>%
#   mutate(Tumor_type = "Osteosarcoma (OSA)", Organization = "Colorado State University", Status = "Tumor", Sample_status = NA, Symbol = NA, Survival_status = NA, Tumor_grade = NA)

# # PRJEB24200 sanger
# Bioproject <- "PRJEB24200"; srarunlist <- paste0(meta_dir, "/", Bioproject, "_SraRunTable.txt")
# df_PRJEB24200 <- read.csv(srarunlist, stringsAsFactors = F) %>%
#   filter(Assay.Type == "WXS" & LibraryLayout == "PAIRED" & Organism == "Canis lupus familiaris") %>%
#   select(Run,BioProject,BioSample, Subject_ID, Platform,ReleaseDate,PHENOTYPE) %>%
#   rename(Run_ID = Run, Sample_ID = BioSample, Case_ID = Subject_ID, Tissue = PHENOTYPE) %>%
#   mutate(Tumor_type = "Mast Cell Tumours", Organization = "SC", Status = NA, Breed = NA, Symbol = NA, Sample_status = NA, Age = NA, Sex = NA , Disease_stage = NA, Survival_status = NA, Tumor_grade = NA) %>%
#   group_by(BioProject,Sample_ID,Case_ID,Status,Breed,Platform,ReleaseDate,Age,Sex,Tissue,Tumor_type,Organization,Sample_status,Symbol,Survival_status,Tumor_grade,Disease_stage) %>% 
#   arrange(Run_ID) %>% summarise(Run_ID = paste(Run_ID, collapse="-")) %>%
#   mutate(Case_ID = gsub("(a|b|c|d)$","_\\1", Case_ID, perl = T)) %>%
#   separate(Case_ID, c("Case_ID","Status"), sep = "_")

all_new_datasets <- rbind(df_PRJDB16014,df_PRJNA891496,df_PRJEB53653,df_PRJDB10211,df_PRJNA786469,df_PRJNA752630,df_PRJNA701141,df_PRJNA695534,df_PRJNA680382,df_PRJNA677995,df_PRJNA630029)
all_new_datasets <- all_new_datasets %>% mutate(Breed = str_to_title(Breed)) %>%
  mutate(Age = gsub(" years","",Age)) %>% mutate(Age = round(as.numeric(Age), digits = 1))

write.csv(all_new_datasets, file = paste0(meta_dir,"/data_collection_new_pre.csv"), quote = TRUE, na = "NA", row.names = F)

# meta data for breed prediction
df_old_meta <- read.csv(file = paste0(meta_dir,"/data_collection_old.csv"), stringsAsFactors = F)

df_new_meta_normal <- all_new_datasets %>% 
  filter(Status != "Cell Line") %>% 
  filter(Status == "Normal")
df_new_meta_pairedtumor <- all_new_datasets %>% 
  filter(Status != "Cell Line") %>% 
  filter(Status == "Tumor") %>% filter(Case_ID %in% df_new_meta_normal$Case_ID)
df_new_meta <- rbind(df_new_meta_normal, df_new_meta_pairedtumor) %>%
  select(c("Case_ID","Run_ID","Tumor_type","Status","Breed","Sample_status","BioProject","Organization")) %>%
  mutate(DiseaseAcronym2 = case_when(grepl("Non-Tumor", Tumor_type, ignore.case = T) ~ "NonT",
                                    grepl("Mammary", Tumor_type, ignore.case = T) ~ "MT",
                                    grepl("histiocytic sarcoma", Tumor_type, ignore.case = T) ~ "HS",
                                    grepl("Oral Melanoma", Tumor_type, ignore.case = T) ~ "OM",
                                    grepl("B-cell Lymphoma", Tumor_type, ignore.case = T) ~ "BCL",
                                    grepl("Osteosarcoma", Tumor_type, ignore.case = T) ~ "OSA",
                                    grepl("hemangiosarcoma", Tumor_type, ignore.case = T) ~ "HSA"
                                    )) %>%
  mutate(Tumor_Type = case_when(grepl("Non-Tumor", Tumor_type, ignore.case = T) ~ "NonT",
                                grepl("Mammary", Tumor_type, ignore.case = T) ~ "MT",
                                grepl("histiocytic sarcoma", Tumor_type, ignore.case = T) ~ "HS",
                                grepl("Oral Melanoma", Tumor_type, ignore.case = T) ~ "OM",
                                grepl("B-cell Lymphoma", Tumor_type, ignore.case = T) ~ "LYM",
                                grepl("Osteosarcoma", Tumor_type, ignore.case = T) ~ "OSA",
                                grepl("hemangiosarcoma", Tumor_type, ignore.case = T) ~ "HSA"
                                )) %>%
  mutate(Breed = gsub("\\\\,.*","",Breed,perl = T)) %>%
  mutate(Breed = gsub(" Dog","",Breed,perl = T)) %>%
  mutate(Breed_info = case_when(grepl("mix",Breed,ignore.case = T) ~ "Mixed", 
                                grepl("Labrador",Breed) ~ "Labrador Retriever",
                                is.na(Breed) ~ "No breed provided",
                                grepl("American Cocker Spaniel",Breed) ~ "Cocker Spaniel",
                                TRUE ~ Breed)) %>%
  mutate(Code = case_when(Organization == "Ghent University" ~ "GU",
                          Organization == "The University of Tokyo (UT-VMS)" ~ "UTok",
                          Organization == "Uppsala University" ~ "UU",
                          Organization == "CNRS - University Rennes" ~ "URen",
                          Organization == "University of Padua" ~ "UPad",
                          Organization == "Cornell University" ~ "CU",
                          Organization == "University of Missouri" ~ "UMis",
                          Organization == "Translational Genomics Research Institute" ~ "TGen",
                          Organization == "University of Helsinki" ~ "UHel",
                          )) %>%
  mutate(Symbol = paste(DiseaseAcronym2,Code,sep = " ")) %>%
  rename(Bioproject = BioProject, Sample_ID = Run_ID) %>%
  select(-c("Tumor_type","Breed","Organization","Code"))

# check breed frequency
# table(df_new_meta$Breed_info)[order(table(df_new_meta$Breed_info), decreasing = T)]

write.csv(df_new_meta, file = paste0(meta_dir,"/data_collection_new.csv"), quote = TRUE, na = "NA", row.names = F)

## merge new datasets with old datasets

all_datasets <- rbind(df_new_meta,df_old_meta)
table(all_datasets$Breed_info)[order(table(all_datasets$Breed_info), decreasing = T)]

# 26 Labrador in old datasets but not included in the paper??

