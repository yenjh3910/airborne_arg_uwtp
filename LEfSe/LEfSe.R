# ARG_lefse.R
# Prepare ARG input file for LEfSe analysis on galaxy platform

# Import library
library(tidyverse)
library(openxlsx)
# Read ARG_subtype file (Something wrong with read.table, so read by )
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
# Split subtype by "_"
arg_subtype <- arg_subtype %>% 
  separate(subtype, c("type","sample_type"), sep = "__")
# Column to first row
new_row <- colnames(arg_subtype)
arg_subtype <- rbind(new_row,arg_subtype)
arg_subtype <- arg_subtype[,-1]
# Change first row
arg_subtype[1,1] <- "sample_type"
arg_subtype[1,2:6] <- "ARP"
arg_subtype[1,7:11] <- "AT"
arg_subtype[1,12:16] <- "ODP"
# Replace first row back to column name
colnames(arg_subtype) <- arg_subtype[1,]
arg_subtype <- arg_subtype[-1,]
# # Save LEfSE input file
# write.table(arg_subtype,"../../airborne_arg_uwtp_result/LEfSe/arg_lefse_3sampletype_input.csv",
#           quote= FALSE, sep = "\t",row.names=FALSE)

######## Upward is three sample type involing AT,ARP,ODP ########## 
######## Downward is two sample type involing aeration_area,outdoor_area ########## 

colnames(arg_subtype)[2:11] <- "aeration_area"
colnames(arg_subtype)[12:16] <- "outdoor_area"
arg_subtype_2sampletype <- arg_subtype

# # Save LEfSE input file
# write.table(arg_subtype_2sampletype,"../../airborne_arg_uwtp_result/LEfSe/arg_lefse_2sampletype_input.csv",
#             quote= FALSE, sep = "\t",row.names=FALSE)

###################################################################
###################################################################
# Preprocess of plotting
## 3 sample type
all_sampletype_lda <- read.table("../../airborne_arg_uwtp_result/LEfSe/Galaxy_LEfSe_3sampletype.txt",
                                 sep = "\t")
colnames(all_sampletype_lda) <- c("ARG","LogMaxMean","Class","LDA","pValue")

p <- all_sampletype_lda %>% drop_na(LDA) %>% 
  filter(LDA > 3.75) %>% 
  mutate(LDA = if_else(Class == "AT", -1*LDA, LDA),
         ARG = fct_reorder(ARG,LDA),
         ARG = fct_reorder(ARG,Class,.desc = TRUE)) %>% 
  ggplot(aes(x=LDA, y=ARG, fill = Class)) +
  geom_col(alpha = 0.6) +
  theme_bw() +
  labs(y = NULL, x = expression(LDA~Score~(log[10]))) +
  guides(fill=guide_legend(title="")) +
  scale_fill_discrete(name = 'Sample', 
                      labels = c(expression(Aeration~tank~PM[2.5]),
                                 expression(Aeration~tank),
                                 expression(Outdoor~PM[2.5]))) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position="top",
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) #transparent legend 
  
print(p)
# ggsave("ARG_LeFSe_all_sample_type.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/LEfSe",
#        width = 7, height = 5,
#        units = "in", bg='transparent') # save to png format


## 2 sample type
two_sampletype_lda <- read.table("../../airborne_arg_uwtp_result/LEfSe/Galaxy_LEfSe_2sampletype.txt",
                                 sep = "\t")
colnames(two_sampletype_lda) <- c("ARG","LogMaxMean","Class","LDA","pValue")

p <- two_sampletype_lda %>% 
  drop_na(LDA) %>%
  filter(LDA > 3.75) %>% 
  mutate(LDA = if_else(Class == "aeration_area", -1*LDA, LDA),
         ARG = fct_reorder(ARG,LDA)) %>% 
  ggplot(aes(x=LDA, y=ARG, fill = Class)) +
  geom_col(alpha = 0.6) +
  theme_bw() +
  labs(y = NULL, x = expression(LDA~Score~(log[10]))) +
  guides(fill=guide_legend(title="")) +
  scale_fill_discrete(name = 'Sample', 
                      labels = c(expression(Aeration~tank~+~Aeration~tank~PM[2.5]),
                                 expression(Outdoor~PM[2.5]))) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position="top",
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) #transparent legend bg


print(p)
  
# ggsave("ARG_LeFSe_two_sample_type.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/LEfSe",
#        width = 7, height = 5,
#        units = "in", bg='transparent') # save to png format