# ARG_subtype_compare.R

# Import library
library(tidyverse)
library(openxlsx)
library(RColorBrewer)

# Read ARG_subtype file (Something wrong with read.table, so read by )
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
# Remove ourdoor sample
arg_subtype <- arg_subtype[,-(12:16)]
# Gather dataframe
gather_arg_subtype <- gather(arg_subtype, key = "sample", value = "copy_per_cell", 
                             ARP1:AT5) # Transform to gather format
# Add sample type column
gather_arg_subtype$sample_type <- gather_arg_subtype$sample
gather_arg_subtype$sample_type <- gsub("1|2|3|4|5","",gather_arg_subtype$sample_type)
# Calculate mean value of ARG subtype in sample type
average_arg <- gather_arg_subtype %>%
  group_by(sample_type, subtype) %>%
  summarise_at(vars(copy_per_cell), funs(mean))
# Spread dataframe
spread_average_arg <- spread(average_arg, key = "sample_type", value = "copy_per_cell")


# Filter ARG found more abundance in ARP than AT
higher_arg_ARP <- spread_average_arg %>% filter(ARP > AT)
# Remove 0 ARG in AT
higher_arg_ARP <- higher_arg_ARP %>% filter(!(AT == 0))
# Add ratio column
higher_arg_ARP$ratio <- higher_arg_ARP$ARP/higher_arg_ARP$AT
# Add difference coulmn
higher_arg_ARP$difference <- higher_arg_ARP$ARP - higher_arg_ARP$AT
# Descending order
higher_arg_ARP <- higher_arg_ARP %>% arrange(desc(ratio))
# Seperate type and subtype
higher_arg_ARP <- higher_arg_ARP %>% separate(subtype, c("type","subtype"), 
           sep = "__")
# Covert first letter to uppercase
higher_arg_ARP$type <- str_to_title(higher_arg_ARP$type)
## Change specific ARG type
higher_arg_ARP$type[higher_arg_ARP$type == "Macrolide-Lincosamide-Streptogramin"] <- "MLS"
higher_arg_ARP$type[higher_arg_ARP$type == "Beta_lactam"] <- "Beta-lactam"

#####################################################
### Fiter by difference (Select by yourself) ###
higher_arg_ARP <- higher_arg_ARP %>% filter(difference > 1e-3)
### Filter ARG on the figure (Select by yourself) ###
higher_arg_ARP <- higher_arg_ARP %>% filter(!(ratio > 300))
#higher_arg_ARP <- higher_arg_ARP %>% filter(ratio >= 10)
#####################################################

# Change the ARG subtype order on the figure
higher_arg_ARP$subtype <- factor(higher_arg_ARP$subtype,
                  levels = higher_arg_ARP$subtype[order(higher_arg_ARP$ratio, 
                                                        decreasing = TRUE)])
# Set color
library(RColorBrewer)
display.brewer.all()
nb.cols <- 13
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
# Visualization
p <- ggplot(higher_arg_ARP, aes(x=subtype, y=ratio, fill=type)) +
  geom_bar(stat="identity") + 
  xlab("") + ylab(expression(ARG~ratio~(Aeration~tank~PM[2.5]/Aeration~tank))) + 
  theme_bw() +
  scale_fill_manual(values = mycolors) +
  guides(fill=guide_legend(title="ARG type",ncol=3,byrow=TRUE)) + 
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = c(0.7, 0.65),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        # panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent')) #transparent legend bg

print(p)

# ggsave("ARG_compare.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 14, height = 9,
#        units = "in", bg='transparent') # save to png format




###################################################
####### Combine higher_arg_ARP to mechanism #######
###################################################
# Get "gather_arg_mechanism"
source("ARG_mechanism.R", 
       encoding = 'utf-8', # If windows CP950
       echo = T # 回傳在該程式執行的程式碼與結果
      )
# Select necessary column from "gather_arg_mechanism"
gather_arg_mechanism <- gather_arg_mechanism %>% select(subtype,mechanism)
gather_arg_mechanism <-unique(gather_arg_mechanism )
# Join "higher_arg_ARP" and "gather_arg_mechanism"
higher_arg_ARP <- left_join(higher_arg_ARP,gather_arg_mechanism, by = "subtype")
# Change the ARG subtype order on the figure
higher_arg_ARP$subtype <- factor(higher_arg_ARP$subtype,
                                 levels = higher_arg_ARP$subtype[order(higher_arg_ARP$ratio, 
                                                                       decreasing = TRUE)])
# Select color fill
RColorBrewer::display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")
# Visualization
p <- ggplot(higher_arg_ARP, aes(x=subtype, y=ratio, fill=mechanism)) +
  geom_bar(stat="identity") + 
  xlab("") + ylab(expression(ARG~ratio~(Aeration~tank~PM[2.5]/Aeration~tank))) + 
  theme_bw() +
  scale_fill_manual(values=c("#FB8072", "#80B1D3", "#FDB462",
                             "#B3DE69","#BEBADA","#FCCDE5","#D9D9D9")) +
  guides(fill=guide_legend(title="Mechanism",ncol=2,byrow=TRUE)) + 
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = c(0.7, 0.7),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        # panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent')) #transparent legend bg

print(p)

# ggsave("ARG_mechanism_compare.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 14, height = 9,
#        units = "in", bg='transparent') # save to png format
