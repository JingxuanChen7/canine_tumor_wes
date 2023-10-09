#!/usr/bin/env Rscript
# parse arguments from command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript --vanilla reformat_meta_text.R <in> <out>", call.=FALSE)
} else if (length(args) == 2) {
    in_meta <- args[1]
    out_meta <- args[2]

}

# df <- read.csv("/home/jc33471/canine_tumor_wes/metadata/data_collection_old.csv", stringsAsFactors = T)
df <- read.csv(in_meta, stringsAsFactors = T)
df_out <- df[c(2,1,4,7,5,10)]
colnames(df_out) <- c("Sample_id","SampleName","DiseaseAcronym","Breed","Status","The_reason_to_exclude")
table_SampleName <- table(df_out$SampleName)
df_out <- df_out[df_out$SampleName %in% names(table_SampleName)[table_SampleName>1],]
# write.table(df_out, file = "/scratch/jc33471/canine_tumor_test/breed_prediction/breed_prediction_metadata.txt",row.names = F,sep = "\t",quote = F)
write.table(df_out, file = out_meta,row.names = F,sep = "\t",quote = F)
