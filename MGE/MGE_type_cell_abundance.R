# MGE_type_cell_abundance.R
# Calculate and visualize MGE type/cell in each sample

# Import library
library(tidyverse)
library(stringr)
library(tibble)

# Read ARG_type file
mge_type <- read.table("D:/ARG_project/Shell/args_oap/MGE/AA_blast_stage_two_output/normalized_cell.type.txt",
                       header = TRUE, sep = "")

# Preview raw data
gather_mge_type <- gather(mge_type, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5) # Transform to gather format
## Add sample type column
gather_mge_type$sample_type <- NA
for (i in 1:nrow(gather_mge_type)) {
  if (grepl("^AT",gather_mge_type$sample[i])) {
    gather_mge_type$sample_type[i] <- "AT"
  } 
  if (grepl("^ARP",gather_mge_type$sample[i])) {
    gather_mge_type$sample_type[i] <- "ARP"
  } 
  if (grepl("^ODP",gather_mge_type$sample[i])) {
    gather_mge_type$sample_type[i] <- "ODP"
  }
}

## Visualization
gather_mge_type$sample_type <- factor(gather_mge_type$sample_type, levels = c("AT", "ARP", "ODP"), 
                                      labels = c(expression(Aeration~tank), 
                                                 expression(Aeration~tank~PM[2.5]), 
                                                 expression(Outdoor~PM[2.5])
                                      )) # Change facet title
ggplot(gather_mge_type, aes(x = sample, y = copy_per_cell, fill = type)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~sample_type, scales = "free_x", labeller = label_parsed) + 
  theme_bw() + 
  xlab("") + ylab("Relative abundance (MGEs/cell)")
