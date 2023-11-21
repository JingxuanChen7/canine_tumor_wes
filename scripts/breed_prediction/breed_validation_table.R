library(dplyr)

df_assignment <- read.table("/scratch/jc33471/canine_tumor/breed_prediction/output_exclude_WGS/assignment_clusters.txt", header = T, sep = "\t",stringsAsFactors = F)

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
    as.numeric(rownames(.)) >= 651 & as.numeric(rownames(.)) <= 680 ~ "Boxer",
    TRUE ~ NA
  ))
