# MGE_type_cell_abundance.R
# Calculate and visualize MGE type/cell in each sample

# Import library
library(tidyverse)
library(stringr)
library(tibble)

# Read MGE_type file
mge_type <- read.table("../../airborne_arg_uwtp_result/args_oap/MGE/AA_stage_two_output/normalized_cell.type.txt",
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


# Classify minimum MGE type to others
average_mge <- gather_mge_type %>%
  group_by(sample_type, type) %>%
  summarise_at(vars(copy_per_cell), funs(mean)) # Calculate mean value of MGE type in sample type

mge_order <- unique(
  average_mge[
    order(average_mge$copy_per_cell,
          decreasing = T),
  ]$type) # unique MGE type order

mge_type <- mge_type %>% 
  arrange(factor(type, levels = mge_order)) # Order MGE type in spread format
## Convert the first column into row name
rownames(mge_type) <- mge_type[,1]
mge_type[,1] <- NULL
## Calculate others by summing  minimum mge
other_mge <- colSums(mge_type[11:40, ])
mge_type <- rbind(mge_type, other_mge)
rownames(mge_type)[rownames(mge_type) == "41"] <- "Others" 
mge_type <- mge_type[-(11:40),] # delete minimum mge type
mge_type <- tibble::rownames_to_column(mge_type, "type") # Convert row name back to first column



# Plot
## Transform to gather format
gather_mge_type <- gather(mge_type, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5)
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
## Change facet title
gather_mge_type$sample_type <- factor(gather_mge_type$sample_type, levels = c("AT", "ARP", "ODP"), 
                                      labels = c(expression(Aeration~tank), 
                                                 expression(Aeration~tank~PM[2.5]), 
                                                 expression(Outdoor~PM[2.5])
                                      ))
## Change specific MGE type
gather_mge_type$type[gather_mge_type$type == "insertion_element_IS91"] <- "IS91"
gather_mge_type$type[gather_mge_type$type == "transposase"] <- "Transposase"
gather_mge_type$type[gather_mge_type$type == "plasmid"] <- "Plasmid"
gather_mge_type$type[gather_mge_type$type == "integrase"] <- "Integrase"
## Order MGE type
mge_order
gather_mge_type$type <- factor(gather_mge_type$type, 
                               levels = c("Transposase", "IS91",
                                          "istB", "Integrase",
                                          "istA", "qacEdelta",
                                          "tniB", "tniA",
                                          "istB1", "Plasmid",
                                          "Others"))


p <- ggplot(gather_mge_type, aes(x = sample, y = copy_per_cell, fill = type)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~sample_type, scales = "free_x", labeller = label_parsed) + 
  theme_bw() + 
  xlab("") + ylab("Relative abundance (MGEs/cell)") + 
  guides(fill=guide_legend(title="MGE type")) + 
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
        legend.background = element_rect(fill='transparent')) + #transparent legend bg
        scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                                   "#FDB462","#B3DE69","#FCCDE5", "#FFED6F","#BC80BD",
                                   "#D9D9D9"))

print(p)

# ggsave("MGE_type_cell.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/MGE",
#        width = 7, height = 5,
#        units = "in", bg='transparent') # save to png format