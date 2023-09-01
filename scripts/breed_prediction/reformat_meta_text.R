
df <- read.csv("/home/jc33471/canine_tumor_wes/metadata/data_collection_old.csv", stringsAsFactors = T)
df_out <- df[c(2,1,4,7,5)]
colnames(df_out) <- c("Sample_id","SampleName","DiseaseAcronym","Breed","Status")
table_SampleName <- table(df_out$SampleName)
df_out <- df_out[df_out$SampleName %in% names(table_SampleName)[table_SampleName>1],]
write.table(df_out, file = "/scratch/jc33471/canine_tumor_test/breed_prediction/breed_prediction_metadata.txt",row.names = F,sep = "\t",quote = F)
