# ARG_subtype_venn_diagram.R

library(tidyverse)
library(openxlsx)
library(RColorBrewer)
library(ggvenn)

# Read ARG_subtype file (Something wrong with read.table, so read by )
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
# Merge sample to sample type
arg_sampleType <- arg_subtype %>% mutate(AT = AT1 + AT2 + AT3 + AT4 + AT5)
arg_sampleType <- arg_sampleType %>% mutate(ARP = ARP1 + ARP2 + ARP3 + ARP4 + ARP5)
arg_sampleType <- arg_sampleType %>% mutate(ODP = ODP1 + ODP2 + ODP3 + ODP4 + ODP5)
arg_sampleType <- arg_sampleType %>% select(subtype, AT, ARP, ODP)
# Statistic
ARG_NotIn_AT_ButIn_ARP <- arg_sampleType %>% # Pick ARG not in AT but in ARP
                          filter(AT == 0) %>% filter(!(ARP == 0))
subsample_ARP_not_AT <- arg_subtype %>%  # filter by ARG_NotIn_AT_ButIn_ARP
                        filter(subtype %in% ARG_NotIn_AT_ButIn_ARP[,1])
g_subsample_ARP_not_AT <- gather(subsample_ARP_not_AT, key = "sample",
                                 value = "copy_per_cell",ARP1:ODP5)
g_subsample_ARP_not_AT <- g_subsample_ARP_not_AT %>% group_by(sample) %>% 
                          mutate(sum = sum(copy_per_cell)) %>% 
                          select(sample,sum) %>% unique() # Calculate sum ARG
## Add sample type
g_subsample_ARP_not_AT$sample_type <- g_subsample_ARP_not_AT$sample
g_subsample_ARP_not_AT$sample_type <- gsub("1|2|3|4|5","",g_subsample_ARP_not_AT$sample_type)
## Calculate mean and sum
g_subsample_ARP_not_AT <- g_subsample_ARP_not_AT %>% group_by(sample_type) %>% 
                          mutate(mean = mean(sum)) %>%  mutate(sd = sd(sum)) 
sum_ARP_not_AT  <- g_subsample_ARP_not_AT %>% filter(sample_type == "ARP") %>% 
                   summarise(sum(sum))
sum_ARP_not_AT[1,2]/sum(arg_sampleType$ARP) # ratio of ARG in ARP but not in AT 
# First column to row name
arg_subtype <- arg_sampleType
row.names(arg_subtype) <- arg_subtype$subtype
arg_subtype <- arg_subtype[,-1]
# Replace value to arg
arg_subtype <- arg_subtype %>% arrange(desc(AT))
arg_subtype$AT[!(arg_subtype$AT == 0)] <- row.names(arg_subtype)
arg_subtype <- arg_subtype %>% arrange(desc(ARP))
arg_subtype$ARP[!(arg_subtype$ARP == 0)] <- row.names(arg_subtype)
arg_subtype <- arg_subtype %>% arrange(desc(ODP))
arg_subtype$ODP[!(arg_subtype$ODP == 0)] <- row.names(arg_subtype)
arg_subtype[arg_subtype == 0] <- NA # Replace 0 to NA
arg_subtype <- as.list(arg_subtype) # Transform to list

# Select fill color
display.brewer.all()
brewer.pal(6, "Set3")

# Plot
p <- ggvenn(arg_subtype,
       fill_color = c("#FB8072", "#80B1D3", "#B3DE69"), 
       fill_alpha = 0.45,
       stroke_size = 0.3, stroke_alpha = 1,
       set_name_color = "white", set_name_size = 0.00000001, text_size = 4)
p <- p + annotate("text",x=-0.7,y=1.8,label="Aeration~tank",parse=T)
p <- p + annotate("text",x=0.75,y=1.785,label="Aeration~tank~PM[2.5]",parse=T)
p <- p + annotate("text",x=0,y=-1.8,label="Outdoor~PM[2.5]",parse=T)
print(p)

# ggsave("ARG_venn.png", p, path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 3.5, height = 3.2, units = "in") # save to png format
