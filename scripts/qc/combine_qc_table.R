library(dplyr)
library(ggplot2)

project_dir <- "/home/jc33471/canine_tumor_wes"
meta <- read.csv("/home/jc33471/canine_tumor_wes/metadata/data_collection_new.csv", stringsAsFactors = F)
qc_dir <- "/home/jc33471/canine_tumor_wes"

bwa_out <- read.table(paste0(qc_dir,"/results/QC/Total_WES_BWA_CDS.combined.txt"), header = T, stringsAsFactors = F, sep = "\t") %>%
  distinct()
mapq_out <- read.table(paste0(qc_dir,"/results/QC/WES_Total_Mapping_quality.combined.txt"), header = T, stringsAsFactors = F, sep = "\t")%>%
  distinct()
randomness_out <- read.table(paste0(qc_dir,"/results/QC/Total_WES_Randomness_Summary.combined.txt"), header = T, stringsAsFactors = F, sep = "\t")%>%
  distinct()

qc_combine <- left_join(bwa_out, mapq_out, by = c("File_name", "CancerType"="Cancer_Type","Status")) %>%
  left_join(randomness_out, by = c("File_name", "CancerType"="Cancer_type","Status")) %>%
  select(c("File_name","Total_Pairs","Uniq_mapped_rate","fra30","Uniq_Exonic_region_mapped_rate","average","rmse"))

# long branch samples from phylogenetic analaysis
sample_to_exclude <- c("SRR17139289","SRR17139287","SRR17139267","SRR17139221","SRR17139298",
  "SRR17139236","SRR17139256","SRR17139258","SRR17139253","SRR17139249",
  "SRR17139240","SRR17139247","SRR17139342","SRR17139320","SRR17139242")

wes_qc_meta <- left_join(meta, qc_combine, by = c("Sample_ID" = "File_name")) %>%
  mutate(Reason_to_exclude = case_when( Total_Pairs < 5000000 ~ "Total_Pair_reads < 5M",
                                        Uniq_mapped_rate < 0.5 ~ "Uniq mapping rate < 0.5",
                                        Uniq_Exonic_region_mapped_rate < 0.3 ~ "Unique CDS mapping rate < 0.3",
                                        average < 30 ~ "coverage < 30",
                                        #Sample_ID %in% sample_to_exclude ~ "phylogeny outlier",
                                        is.na(Total_Pairs) ~ "fail pipeline",
                                        TRUE ~ "Pass QC"
                                      )) 


out_table <- wes_qc_meta %>% select(colnames(meta)) %>%
  mutate(Case_ID = case_when(
    (Bioproject == "PRJNA786469" & grepl("^CC",Case_ID)) ~ paste0(Case_ID,"_URen"),
    TRUE ~ Case_ID
  ))


write.csv(out_table, "/home/jc33471/canine_tumor_wes/metadata/data_collection_new_renamed.csv", row.names = F)


######## plot dot plot distributions ######## 

# set color palette, text size etc.
file_color <- c("Normal"="darkblue","Tumor"="red3")
regular.text <- element_text(colour="black",size=20);

# Sequence_read_pairs
wes_qc_meta %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Symbol),
             y=as.numeric(Total_Pairs)/1000000,fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Sequence read pairs\nin millions")+
  scale_y_continuous(breaks = c(0,100,200))+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.position="none", 
        legend.title=regular.text, 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=5, linetype="longdash", color = "yellow4", size = 0.7)
ggsave(file = paste0(qc_dir,"/results/QC/Sequence_read_pairs.pdf"), device = "pdf", width = 8, height = 8)

# Uniq mapping rate
wes_qc_meta %>% 
  ggplot(aes(x=factor(Symbol),
             y=as.numeric(Uniq_mapped_rate),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Uniquely & concordant\nmapped rate")+
  scale_y_continuous(breaks = c(0.5,0.75,1.0))+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+  
  coord_cartesian(ylim=c(0,1))+  
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=0.6, linetype="longdash", color = "yellow4", size = 0.7)
ggsave(file = paste0(qc_dir,"/results/QC/Uniq_mapping_rate.pdf"), device = "pdf", width = 8, height = 8)

# fraction of mapping quality > 30
wes_qc_meta %>% 
  ggplot(aes(x=factor(Symbol),
             y=as.numeric(fra30),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Fraction of\nmapping quality > 30")+
  scale_y_continuous(breaks = c(0.25,0.5,0.75,1.0))+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+  
  coord_cartesian(ylim=c(0,1))+  
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))
ggsave(file = paste0(qc_dir,"/results/QC/Fraction_mapQ30.pdf"), device = "pdf", width = 8, height = 8)

# Uniq CDS mapping rate
wes_qc_meta %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  #filter(as.numeric(Uniquely_mapped_rate) >=0.6)%>% 
  ggplot(aes(x=factor(Symbol),
             y=as.numeric(Uniq_Exonic_region_mapped_rate),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("CDS targeting rate")+
  scale_y_continuous(breaks = c(0,0.4,0.8))+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+  coord_cartesian(ylim=c(0,0.8))+  
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=0.3, linetype="longdash", color = "yellow4", size = 0.7)#+
  #geom_hline(yintercept=0.2, linetype="longdash", color = "yellow3", size = 0.7)
ggsave(file = paste0(qc_dir,"/results/QC/CDS_mapping_rate.pdf"), device = "pdf", width = 8, height = 8)


# mean coverage
wes_qc_meta %>% 
  ggplot(aes(x=factor(Symbol),
             y=as.numeric(average),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Mean coverage")+
  scale_y_continuous(breaks = c(0,100,200))+
  coord_cartesian(ylim=c(0,200))+ 
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+   
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=30, linetype="longdash", color = "yellow4", size = 0.7)
ggsave(file = paste0(qc_dir,"/results/QC/mean_coverage.pdf"), device = "pdf", width = 8, height = 8)


# RMSE
wes_qc_meta %>% 
  ggplot(aes(x=factor(Symbol),
             y=as.numeric(rmse),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("RMSE of\nsequence coverage")+
  scale_y_continuous(breaks = c(0,0.005,0.01))+
  coord_cartesian(ylim=c(0,0.01))+ 
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+   
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=0.01, linetype="longdash", color = "yellow4", size = 0.7)
ggsave(file = paste0(qc_dir,"/results/QC/rmse.pdf"), device = "pdf", width = 8, height = 8)


