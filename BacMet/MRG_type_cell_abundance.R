# MRG_type_cell_abundance.R
# Calculate and visualize MRG type/cell in each sample

# Import library
library(tidyverse)
library(stringr)
library(tibble)

# Read MGE_type file
mrg_type <- read.table("D:/ARG_project/Shell/args_oap/BacMet/stage_two_output_evalue-5_id70/normalized_cell.type.txt",
                       header = TRUE, sep = "\t")

# Preview raw data
gather_mrg_type <- gather(mrg_type, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5) # Transform to gather format
## Add sample type column
gather_mrg_type$sample_type <- NA
for (i in 1:nrow(gather_mrg_type)) {
  if (grepl("^AT",gather_mrg_type$sample[i])) {
    gather_mrg_type$sample_type[i] <- "AT"
  } 
  if (grepl("^ARP",gather_mrg_type$sample[i])) {
    gather_mrg_type$sample_type[i] <- "ARP"
  } 
  if (grepl("^ODP",gather_mrg_type$sample[i])) {
    gather_mrg_type$sample_type[i] <- "ODP"
  }
}

# See MGE order
average_mrg <- gather_mrg_type %>%
  group_by(sample_type, type) %>%
  summarise_at(vars(copy_per_cell), funs(mean)) # Calculate mean value of MRG type in sample type

mrg_order <- unique(
  average_mrg[
    order(average_mrg$copy_per_cell,
          decreasing = T),
  ]$type) # unique MGE type order
## Order MGE type
mrg_order # View
gather_mrg_type$type <- factor(gather_mrg_type$type, 
                               levels = c("Copper (Cu)", "Multi-metal",
                                          "Iron (Fe)", "Arsenic (As)",
                                          "Silver (Ag)", "Mercury (Hg)",
                                          "Chromium (Cr)"))


## Visualization
gather_mrg_type$sample_type <- factor(gather_mrg_type$sample_type, levels = c("AT", "ARP", "ODP"), 
                                      labels = c(expression(Aeration~tank), 
                                                 expression(Aeration~tank~PM[2.5]), 
                                                 expression(Outdoor~PM[2.5])
                                      )) # Change facet title

p <- ggplot(gather_mrg_type, aes(x = sample, y = copy_per_cell, fill = type)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~sample_type, scales = "free_x", labeller = label_parsed) + 
  theme_bw() + 
  xlab("") + ylab("Relative abundance (MRGs/cell)") + 
  scale_fill_brewer(palette="Set3") + 
  guides(fill=guide_legend(title="MRG type")) + 
  # scale_y_continuous(expand = c(0, 0)) + # y start at 0
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size=10.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        # panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent')) #transparent legend bg


print(p)

# ggsave("MRG_type_cell.png", p, path = "D:/ARG_project/Figure/MRG",
#        width = 7, height = 5, units = "in", bg='transparent') # save to png format