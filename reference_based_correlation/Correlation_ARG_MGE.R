# Correlation_ARG_MGE

# Import library
library(tidyverse)
library(ggpubr)

# Read type file
arg <- read.table("D:/ARG_project/Shell/args_oap/ARG/stage_two_output/normalized_cell.type.txt",
                       header = TRUE, sep = "")
mge <- read.table("D:/ARG_project/Shell/args_oap/MGE/AA_stage_two_output/normalized_cell.type.txt",
                       header = TRUE, sep = "")

# Calculate Sum of ARG and MGE
## ARG sum
gather_arg <- gather(arg, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5) # Transform to gather format
### Add sample type column
gather_arg$sample_type <- gather_arg$sample
gather_arg$sample_type <- gsub("1|2|3|4|5","",gather_arg$sample_type)
### Sum 
sum_arg <- gather_arg %>%
  group_by(sample) %>%
  summarise_at(vars(copy_per_cell), funs(sum)) # Calculate mean value of ARG type in sample
### Change column name
colnames(sum_arg)[2] <- "ARG"
## MGE sum
gather_mge <- gather(mge, key = "sample", value = "copy_per_cell", 
                     ARP1:ODP5) # Transform to gather format
### Add sample type column
gather_mge$sample_type <- gather_mge$sample
gather_mge$sample_type <- gsub("1|2|3|4|5","",gather_mge$sample_type)
### Sum
sum_mge <- gather_mge %>%
  group_by(sample) %>%
  summarise_at(vars(copy_per_cell), funs(sum)) # Calculate mean value of ARG type in sample
### Add MGE tag
colnames(sum_mge)[2] <- "MGE"

# Join ARG and MGE dataframe
arg_mge_type <- full_join(sum_arg,sum_mge,by="sample")
# Add sample type column
arg_mge_type$sample_type <-arg_mge_type$sample
arg_mge_type$sample_type <- gsub("1|2|3|4|5","",arg_mge_type$sample_type)
# Order sample type
arg_mge_type$sample_type <- factor(arg_mge_type$sample_type, levels = c("AT","ARP","ODP"))

# Plot
p <- ggscatter(arg_mge_type, x = "ARG", y = "MGE",
          color = 'sample_type', shape = 'sample_type', 
          add = "reg.line", conf.int = TRUE) + 
  stat_cor(aes(color= sample_type), show.legend = FALSE) + 
  labs(x = "ARGs / cells", y = "MGEs / cells") + 
  theme(aspect.ratio=1, 
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) +  #transparent legend bg)
  scale_color_manual(name = "", 
                     labels = c(expression(Aeration~tank),
                                expression(Aeration~tank~PM[2.5]),
                                expression(Outdoor~PM[2.5])),
                     values = c("AT" = "#FB8072",
                               "ARP" = "#80B1D3", 
                               "ODP" = "#FDB462")) + 
  scale_fill_manual(name = "", 
                    labels = c(expression(Aeration~tank),
                               expression(Aeration~tank~PM[2.5]),
                               expression(Outdoor~PM[2.5])),
                    values = c("AT" = "#FB8072",
                               "ARP" = "#80B1D3", 
                               "ODP" = "#FDB462"))
print(p)

### Can't hide the original legend, so edit it on PPT afterward !!! ###

# ggsave("ARG_MGE_correlation.png", p, path = "D:/ARG_project/Figure/correlation",
#        width = 7, height = 5, units = "in", bg='transparent') # save to png format