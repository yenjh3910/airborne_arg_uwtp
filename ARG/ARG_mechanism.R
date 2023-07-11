# ARG_mechanism.R

# Import library and structure file
library(tidyverse)
library(openxlsx)
library(stringr)

######### Structure do not contain mechanism column, so fail #################
# # Bind structure file
# single_component <- read.xlsx("D:/ARG_project/airborne_arg_uwtp/ARG/single-component_structure.xlsx", 
#                                sheet = 1)
# two_component <- read.xlsx("D:/ARG_project/airborne_arg_uwtp/ARG/two-component_structure.xlsx", 
#                             sheet = 1)
# multi_component <- read.xlsx("D:/ARG_project/airborne_arg_uwtp/ARG/multi-component_structure.xlsx", 
#                            sheet = 1)
# all_component <- rbind(single_component, two_component, multi_component)
# # Read ARG_subtype file (Something wrong with read.table, so read by )
# arg_subtype <- read.table("D:/ARG_project/Shell/args_oap/ARG/stage_two_output/normalized_cell.gene.txt",
#                           sep = "\t", header = TRUE, quote = '')
# gather_arg_subtype <- gather(arg_subtype, key = "sample", value = "copy_per_cell", 
#                              ARP1:ODP5) # Transform to gather format
# #Join structure file and ARG_subtype file
# gather_arg_subtype <- left_join(gather_arg_subtype, all_component, by = "gene")




################################ Merage by gene id ####################################
# Read ARG_subtype file (Something wrong with read.table, so read by )
# arg_subtype <- read.table("D:/ARG_project/Shell/args_oap/ARG/stage_two_output/normalized_cell.gene.txt",
#                           sep = "\t", header = TRUE, quote = '')
# gather_arg_subtype <- gather(arg_subtype, key = "sample", value = "copy_per_cell", 
#                              ARP1:ODP5) # Transform to gather format
# # Read structure file
# single <- read.xlsx("C:/Users/Yen/Desktop/single-component_structure.xlsx", 
#                               sheet = 1)
# two <- read.xlsx("C:/Users/Yen/Desktop/two-component_structure.xlsx", 
#                            sheet = 1)
# multi <- read.xlsx("C:/Users/Yen/Desktop/multi-component_structure.xlsx", 
#                              sheet = 1)
# all <- rbind(single, two, multi)
# all <- all[,1:5]
# all <- all[,c(1,5)]
# colnames(all)[1] <- "gene"
# gather_arg_subtype <- left_join(gather_arg_subtype, all, by = "gene")
# na_gene <- gather_arg_subtype[!complete.cases(gather_arg_subtype), ]
# na_gene <- spread(na_gene, key = "sample", value = "copy_per_cell")





# Read ARG_subtype file (Something wrong with read.table, so read by )
arg_subtype <- read.table("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.txt",
                          sep = "\t", header = TRUE, quote = '')
gather_arg_subtype <- gather(arg_subtype, key = "sample", value = "copy_per_cell", 
                             ARP1:ODP5) # Transform to gather format
single <- read.xlsx("./ARGs_OAP_beta_structure/single-component_structure.xlsx", 
                    sheet = 1)
two <- read.xlsx("./ARGs_OAP_beta_structure/two-component_structure.xlsx", 
                 sheet = 1)
multi <- read.xlsx("./ARGs_OAP_beta_structure/multi-component_structure.xlsx", 
                   sheet = 1)
# Bind structure file
all <- rbind(single, two, multi)
all <- all[,c(3,5)]
colnames(all)[1] <- "subtype"
all <- unique(all)
# Join structure and subtype file
gather_arg_subtype <- left_join(gather_arg_subtype, all, by = "subtype")
spread_arg_subtype <- spread(gather_arg_subtype, key = "sample", value = "copy_per_cell")
colnames(spread_arg_subtype)[2] <- "mechasinm"
unique(spread_arg_subtype$mechasinm)
# See NA mechanism
na_gene <- gather_arg_subtype[!complete.cases(gather_arg_subtype), ]
na_gene <- spread(na_gene, key = "sample", value = "copy_per_cell")

#####################    Warning     ###########################
# ##############################################################
# # Output arg_mechanism file and edit manually!!!!
# write.csv(spread_arg_subtype, file = "D:/ARG_project/Shell/args_oap/ARG/stage_two_output/arg_mechasism.csv", 
#           row.names = FALSE, quote = FALSE,)
# ##############################################################


# Input back curated arg_mechanism file and adjust mechanism name
arg_mechanism <- read.csv("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/arg_mechanism.csv")
unique(arg_mechanism$mechanism)
arg_mechanism$mechanism <- str_to_sentence(arg_mechanism$mechanism) 
arg_mechanism$mechanism[arg_mechanism$mechanism == 'Enzymatic inactivation'] <- 'Antibiotic inactivation'
arg_mechanism$mechanism[arg_mechanism$mechanism == 'Efflux pump'] <- 'Antibiotic efflux'
arg_mechanism$mechanism[arg_mechanism$mechanism == 'Efflux pump rnd family'] <- 'Antibiotic efflux'
arg_mechanism$mechanism[arg_mechanism$mechanism == 'Reduced permeability'] <- 'Reduced permeability to antibiotic'
# Split subtype by "_"
arg_mechanism <- arg_mechanism %>% 
  separate(subtype, c("type","subtype"), sep = "__")
# Transform to gather format
gather_arg_mechanism <- gather(arg_mechanism, key = "sample", value = "copy_per_cell", 
                             ARP1:ODP5)
# Add sample type column
gather_arg_mechanism$sample_type <- gather_arg_mechanism$sample
gather_arg_mechanism$sample_type <- gsub("1|2|3|4|5","",gather_arg_mechanism$sample_type)

## Order
gather_arg_mechanism$sample_type <- factor(gather_arg_mechanism$sample_type, 
                               levels = c("AT","ARP","ODP"))
gather_arg_mechanism$mechanism <- factor(gather_arg_mechanism$mechanism, 
                                           levels = c("Antibiotic efflux",
                                                      "Antibiotic inactivation",
                                                      "Antibiotic target alteration",
                                                      "Antibiotic target replacement",
                                                      "Antibiotic target protection",
                                                      "Reduced permeability to antibiotic",
                                                      "Others"))
## Merge minor mechanism to other mechanism
gather_arg_mechanism$mechanism <- str_replace(gather_arg_mechanism$mechanism, 
                                       "Reduced permeability to antibiotic", 
                                       "Others")
# Calculate sum abundance of each sample
tmp <- gather_arg_mechanism %>% group_by(mechanism, sample) %>% 
  summarise(abundance = sum(copy_per_cell))
# Add sample type column
tmp$sample_type <- tmp$sample
tmp$sample_type <- gsub("1|2|3|4|5","",tmp$sample)
# Calculate mean and sd
tmp <- tmp %>% group_by(mechanism, sample_type) %>% mutate(mean = mean(abundance))
tmp <- tmp %>% group_by(mechanism, sample_type) %>% mutate(sd = sd(abundance))
# Filter necesssary column for figure
arg_mechanism <- tmp %>% select(mechanism,sample_type,mean,sd) %>% unique()
# Add st_mean for figure
arg_mechanism <- arg_mechanism %>% group_by(sample_type) %>% 
  mutate(sd_mean = mean)
arg_mechanism <- arg_mechanism %>% arrange(sample_type)
## sum ARP st_mean
for (i in c(1:5)) {
  arg_mechanism[{i},5] <-  arg_mechanism[{i},5] + 
    sum(arg_mechanism[{i+1}:6,5])}
## sum AT st_mean
for (i in c(7:11)) {
  arg_mechanism[{i},5] <-  arg_mechanism[{i},5] + 
    sum(arg_mechanism[{i+1}:12,5])}
## sum ODP st_mean
for (i in c(13:17)) {
  arg_mechanism[{i},5] <-  arg_mechanism[{i},5] + 
    sum(arg_mechanism[{i+1}:18,5])}
# Calculate the increase ratio of mechanism between AT & ARP
ARP_mech <- arg_mechanism %>% filter(sample_type == "ARP")
AT_mech <- arg_mechanism %>% filter(sample_type == "AT")
(ARP_mech$mean-AT_mech$mean)/AT_mech$mean
# Select color
library(RColorBrewer)
RColorBrewer::display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")
# Order
arg_mechanism$sample_type <- factor(arg_mechanism$sample_type, 
                                    levels = c("ODP",
                                               "ARP",
                                               "AT"))
# Plot
## Absoulte abundance
p <- ggplot(arg_mechanism, aes(x = sample_type, y = mean, fill = mechanism))+
  geom_bar(stat="identity")+
  geom_errorbar( aes(x=sample_type, ymin=sd_mean-sd, ymax=sd_mean+sd), 
                 width=0.4, colour="black", alpha=0.9, size=0.5) + 
  coord_flip()+
  theme_bw()+
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) +
  xlab("")+ ylab("Relative abundance (ARGs/cell)") +
  guides(fill=guide_legend(title="Mechanism",nrow=3,byrow=TRUE)) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.margin=unit(c(5,1,1,1),"cm"),
        legend.position = c(0.5,1.45),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.8, 'cm'),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) + #transparent legend bg)
  scale_fill_manual(values=c("#FB8072", "#80B1D3", "#FDB462",
                             "#B3DE69","#BEBADA","#FCCDE5"))

print(p)
# ggsave("ARG_mechanism_absoulte_abudnace.png", p, path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 11, height = 5, units = "in") # save to png format


## Relative abundance
p <- ggplot(arg_mechanism, aes(x = sample_type, y = mean, fill = mechanism))+
  geom_bar(stat="identity", position = "fill") +
  coord_flip()+
  theme_bw()+
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) +
  xlab("")+ ylab("Relative abundance (ARGs/cell)") +
  guides(fill=guide_legend(title="Mechanism",nrow=3,byrow=TRUE)) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.margin=unit(c(5,1,1,1),"cm"),
        legend.position = c(0.5,1.45),
        legend.spacing.x = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.8, 'cm'),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) + #transparent legend bg)
  scale_fill_manual(values=c("#FB8072", "#80B1D3", "#FDB462",
                             "#B3DE69","#BEBADA","#FCCDE5","#D9D9D9"))

print(p)
# ggsave("ARG_mechanism_relative_abundance.png", p, path = "../../airborne_arg_uwtp_result/Figure/ARG",
#         width = 11, height = 5, units = "in") # save to png format


# mechanism percentage
arg_proportion <- arg_mechanism %>% group_by(sample_type) %>% mutate(percentage = mean/sum(mean)*100) %>% 
                  select(mechanism,sample_type,percentage)
arg_proportion$percentage[1]-arg_proportion$percentage[7]



# Output mechanism reference
mechanism_reference <- gather_arg_mechanism %>% select(subtype,mechanism) %>% unique()
# write.csv(mechanism_reference, file = "./arg_mechanism_reference.csv",
#           row.names = FALSE, quote = FALSE,)
