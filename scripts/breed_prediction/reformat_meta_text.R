#!/usr/bin/env Rscript
# parse arguments from command line
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript --vanilla reformat_meta_text.R <in> <out> <out_cases>", call.=FALSE)
} else if (length(args) == 3) {
    in_meta <- args[1]
    out_meta <- args[2]
    out_cases <- args[3]

}
## metatable to identify breed specific variants
# in_meta <- "/scratch/jc33471/canine_tumor/breed_prediction/this_wunpaired_meta.csv"
df <- read.csv(in_meta, stringsAsFactors = T)
df_out <- df[c(2,1,4,7,5,10)]
colnames(df_out) <- c("Sample_id","SampleName","DiseaseAcronym","Breed","Status","The_reason_to_exclude")
df_out$Breed <- gsub("Boston Terrier (JAXH)","Boston Terrier",df_out$Breed,fixed = T)
df_out$Breed <- gsub("Labrador$","Labrador Retriever",df_out$Breed,perl = T)

# write.table(df_out, file = "/scratch/jc33471/canine_tumor_test/breed_prediction/breed_prediction_metadata.txt",row.names = F,sep = "\t",quote = F)
write.table(df_out, file = out_meta,row.names = F,sep = "\t",quote = F)

## metatable for cases to identify breeds with n > 10
table_SampleName <- table(df_out$SampleName)
paired <- df_out[df_out$SampleName %in% names(table_SampleName)[table_SampleName>1],]
paired <- paired[paired$Status=="Normal",]
unpaired <- df_out[df_out$SampleName %in% names(table_SampleName)[table_SampleName==1],]
all_cases <- rbind(paired, unpaired)
write.table(all_cases, file = out_cases,row.names = F,sep = "\t",quote = F)
