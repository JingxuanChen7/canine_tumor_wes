library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)
library(data.table)
#colors <- brewer.pal(3, "Set1"); # red, blue, green

## whole criteria 
#Sequence read pairs < 5M	
#Unique CDS mapping rate < 0.3	
#Mean_coverage < 30	
#Callable bases < 10Mb (samples/cases)

########### Photoshop use res = 450 #################

# PAIR <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Figure1\\Original_Data_summary.xlsx",
#                    sheet ="PAIRS")
# 
# callable <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Mutation_rate\\Mutation_rate.xlsx",
#                     sheet ='All_samples_Retro_non_retro')

new_uniq <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan-Cancer-Manuscript/Methods_legends_tables/V6Combined_Summary_Public_Data.xlsx",
                       sheet ='NGS_Data_summary', skip = 3)




exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                      sheet ="Before_Matching_excluded")

total_file <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan-Cancer-Manuscript/Methods_legends_tables/V6Combined_Summary_Public_Data.xlsx",
                         sheet ='NGS_Data_summary', skip = 3)


MC <- total_file %>% 
  filter(Tumor_Type == 'MT') %>% 
  filter(RMSE > 0.002)

  MC %>% 
ggplot(aes(x=0,
           y=as.numeric(RMSE), +
  
  geom_point(size=1.6,shape=20)))


# 
# Combined_total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Pan-Cancer-Manuscript\\Combined_Summary_Public_Data.xlsx",
#                                   sheet ='NGS_Data_summary', skip =0)

Mut_rate <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\final_mut_rate.xlsx",
                       sheet ='Sheet1')



pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\V13F1_and_supplementaryF1.pdf"
    , height=4.98, width=4.84);

regular.text <- element_text(colour="black",size=20);
file_color <- c("darkblue","red3")

# # 
# tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\square-sequence_pairs.tiff",
#      width = 3500, height =3000, units = "px", res = 400)
total_file <- data.table(total_file)
PairMedian <- total_file[,.(Median=median(Total_pairs)),by = .(Tumor_Type)]
UniqMapMedian <- total_file[,.(Median=median(Uniquely_mapped_rate)),by = .(Tumor_Type)]
CDSMedian <- total_file[,.(Median=median(uniq_CDS_region_paris_rates)),by = .(Tumor_Type)]
meanMedian <- total_file[,.(Median=median(mean)),by = .(Tumor_Type)]

total_file %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","LYM","OM", "OSA","HSA","UCL")),
             y=as.numeric(Total_pairs)/1000000,fill=Status,color=Status)) +
  
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Sequence read pairs\nin millions")+
  scale_y_continuous(breaks = c(0,100,200))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.position="top", 
        legend.title=regular.text, 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 20)+
  coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=5, linetype="longdash", color = "yellow4", size = 0.7)

# dev.off()

### Unique Mapping Rates ###
#plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\unique-mapping_pairs.png",
#                   width=2800,height=1800,res=350)
# 

# tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\unique-mapping_pairs.tiff",
#     width = 3500, height =3000, units = "px", res = 400)

#### Unique CDS ####
# tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\unique-CDS-mapping.tiff",
#      width = 3500, height =1900, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(as.numeric(Uniquely_mapped_rate) >=0.6)%>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","LYM","OM", "OSA","HSA","UCL")),
             y=as.numeric(uniq_CDS_region_paris_rates),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  
  ylab("CDS targeting rate")+
  scale_y_continuous(breaks = c(0,0.4,0.8))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 20)+
  coord_cartesian(ylim=c(0,0.8))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=0.3, linetype="longdash", color = "yellow4", size = 0.7)

# dev.off()

### Mean Coverage ###

# tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\coverage-Mean.tiff"
#      ,width = 3500, height =1900, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(as.numeric(Uniquely_mapped_rate) >=0.6) %>% 
  filter(as.numeric(uniq_CDS_region_paris_rates) >=0.3)%>%
  #filter(mean >30) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","LYM","OM", "OSA","HSA","UCL")),
             y=as.numeric(mean),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  scale_y_continuous(breaks = c(0,100,200))+
  ylab("Mean coverage")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=30, linetype="longdash", color = "yellow4", size = 0.7)

# 
# dev.off()



#### RMSE ####

# tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Coverage-RMSE.tiff", 
#      width = 3500, height =1900, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(as.numeric(Uniquely_mapped_rate) >=0.6) %>% 
  filter(as.numeric(uniq_CDS_region_paris_rates) >=0.3)%>% 
  filter(mean>30) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","LYM","OM", "OSA","HSA","UCL")),
             y=as.numeric(RMSE),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  scale_y_continuous(breaks = c(0,0.005,0.01))+
  ylab("RMSE of \nsequence coverage")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.4, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = file_color)+
  scale_fill_manual(values=file_color)+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0,0.01))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))


dev.off()


### change the color of the point ###
# tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\callable-bases.tiff", 
#      width = 3500, height =1900, units = "px", res = 400)
# pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\callablebases.pdf"
#     , height=3.6, width=6.2);

pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\V11callable_bases.pdf"
    , height=4.0, width=4.84);

total_file %>% 
  filter(Status=='Tumor') %>% 
  filter(!Case_ID %in% exclude$Cases) %>% 
  filter(!Callable_bases =='No-Pair')%>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","LYM","OM", "OSA","HSA","UCL")),
             y=as.numeric(Callable_bases)/1000000,color='black')) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Callable bases in millions")+
  scale_y_continuous(breaks = c(0,10,20,30))+
  theme(axis.text=regular.text, 
        axis.title.y=element_text(colour="black",size=20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.4, width = .7,colour = "black")+
  #scale_shape_manual(values = 19)+
  scale_color_manual(values = 'black')+
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))+
  coord_cartesian(ylim=c(5,30))+
  geom_hline(yintercept=10, linetype="longdash", color = "yellow4", size = 0.7)
dev.off()
#scale_x_discrete(labels=panel2label)+
#scale_color_manual(values=wes_palette(n=5, name="Set3")) ## for scatter plot Wes Anderson color palettes
#scale_fill_manual(values=wes_palette(n=5, name="Royal2")) ## for box plot, bar plot Wes Anderson color palettes
#scale_shape_manual(values = 20)





