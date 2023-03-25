# ARG_type_cell_abundance.R
# Calculate and visualize ARG type/cell in each sample

# Import library
library(tidyverse)
library(stringr)
library(tibble)

# Read ARG_type file
arg_type <- read.table("D:/ARG_project/Shell/args_oap/ARG/stage_two_output/normalized_cell.type.txt",
                   header = TRUE, sep = "")


# Preview raw data
gather_arg_type <- gather(arg_type, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5) # Transform to gather format
## Add sample type column
gather_arg_type$sample_type <- NA
for (i in 1:nrow(gather_arg_type)) {
  if (grepl("^AT",gather_arg_type$sample[i])) {
    gather_arg_type$sample_type[i] <- "AT"
  } 
  if (grepl("^ARP",gather_arg_type$sample[i])) {
      gather_arg_type$sample_type[i] <- "ARP"
  } 
  if (grepl("^ODP",gather_arg_type$sample[i])) {
      gather_arg_type$sample_type[i] <- "ODP"
  }
}
## Visualization
gather_arg_type$sample_type <- factor(gather_arg_type$sample_type, levels = c("AT", "ARP", "ODP"), 
                  labels = c(expression(Aeration~tank), 
                             expression(Aeration~tank~PM[2.5]), 
                             expression(Outdoor~PM[2.5])
                             )) # Change facet title
ggplot(gather_arg_type, aes(x = sample, y = copy_per_cell, fill = type)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~sample_type, scales = "free_x", labeller = label_parsed) + 
  theme_bw() + 
  xlab("") + ylab("Relative abundance (ARGs/cell)")


# Classify minimum ARG type to others
average_arg <- gather_arg_type %>%
  group_by(sample_type, type) %>%
  summarise_at(vars(copy_per_cell), funs(mean)) # Calculate mean value of ARG type in sample type

arg_order <- unique(
  average_arg[
    order(average_arg$copy_per_cell,
                    decreasing = T),
    ]$type) # unique ARG type order

arg_type <- arg_type %>% 
  arrange(factor(type, levels = arg_order)) # Order ARG type in spread format

## Convert the first column into row name
rownames(arg_type) <- arg_type[,1]
arg_type[,1] <- NULL

## Calculate others by summing  minimum arg
other_arg <- colSums(arg_type[11:26, ])
arg_type <- rbind(arg_type, other_arg)
rownames(arg_type)[rownames(arg_type) == "27"] <- "others" 
arg_type <- arg_type[-(11:26),] # delete minimum arg type
rownames(arg_type) <- str_to_title(rownames(arg_type)) # Covert first letter to uppercase
arg_type <- tibble::rownames_to_column(arg_type, "type") # Convert row name back to first column


# Plot
## Transform to gather format
gather_arg_type <- gather(arg_type, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5)
## Add sample type column
gather_arg_type$sample_type <- NA
for (i in 1:nrow(gather_arg_type)) {
  if (grepl("^AT",gather_arg_type$sample[i])) {
    gather_arg_type$sample_type[i] <- "AT"
  } 
  if (grepl("^ARP",gather_arg_type$sample[i])) {
    gather_arg_type$sample_type[i] <- "ARP"
  } 
  if (grepl("^ODP",gather_arg_type$sample[i])) {
    gather_arg_type$sample_type[i] <- "ODP"
  }
}
## Change facet title
gather_arg_type$sample_type <- factor(gather_arg_type$sample_type, levels = c("AT", "ARP", "ODP"), 
                                      labels = c(expression(Aeration~tank), 
                                                 expression(Aeration~tank~PM[2.5]), 
                                                 expression(Outdoor~PM[2.5])
                                      ))
## Change specific ARG type
gather_arg_type$type[gather_arg_type$type == "Macrolide-Lincosamide-Streptogramin"] <- "MLS"
gather_arg_type$type[gather_arg_type$type == "Beta_lactam"] <- "Beta-lactam"
## Order ARG type
gather_arg_type$type <- factor(gather_arg_type$type, 
                               levels = c("Multidrug", "Sulfonamide",
                               "MLS","Aminoglycoside", 
                               "Tetracycline", "Rifamycin", 
                               "Beta-lactam", "Bacitracin",
                               "Polymyxin", "Chloramphenicol", 
                               "Others"))
## ggplot
p <- ggplot(gather_arg_type, aes(x = sample, y = copy_per_cell, fill = type)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~sample_type, scales = "free_x", labeller = label_parsed) + 
  theme_bw() + 
  xlab("") + ylab("Relative abundance (ARGs/cell)") + 
  scale_fill_brewer(palette="Set3") + 
  guides(fill=guide_legend(title="ARG type")) + 
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

# ggsave("ARG_type_cell.png", p, path = "D:/ARG_project/Figure/ARG",
#        width = 7, height = 5, units = "in", bg='transparent') # save to png format