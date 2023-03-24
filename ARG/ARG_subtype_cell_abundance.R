# ARG_subtype_cell_abundance.R
# Calculate and visualize ARG subtype/cell in each sample

# Import library
library(tidyverse)
library(stringr)
library(tibble)
library(openxlsx)
library(pheatmap)

# Read ARG_type file (Something wrong with read.table, so read by )
arg_subtype <- read.xlsx("D:/ARG_project/Shell/args_oap/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)

# Preview raw data
gather_arg_subtype <- gather(arg_subtype, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5) # Transform to gather format
gather_arg_subtype_order <- gather_arg_subtype[order(gather_arg_subtype$copy_per_cell, 
                                                     decreasing = T),] # Order by abundance
top_subtype_list<- unique(gather_arg_subtype_order[,1])[1:100] #### Select top100 subtype (Choose)
## View top ARG percentage
top_subtype <- filter(arg_subtype, subtype %in% top_subtype_list) # Assign top subtype dataframe
for (i in c(2:16)) {
  print(paste(colnames(top_subtype)[i],': ',
              sum(top_subtype[,i])/sum(arg_subtype[,i]), 
              sep = ""))
  } # Print top ARG percentage in each sample


top_subtype <- top_subtype %>% 
  separate(subtype, c("type","subtype"), 
           sep = "__") # Split subtype by "_"

top_subtype <- top_subtype[,-1] # Remove ARG type

# Change the first column to row name
row.names(top_subtype) <- top_subtype[,1]
top_subtype <- top_subtype[,-1]


#Plot
top_subtype <- as.matrix(top_subtype) # Convert table to matrix 
log_top_subtype <- log10(top_subtype) # log10 conversion 
log_top_subtype[log_top_subtype == -Inf] <- NA

## Create annotation label
annotation_row = data.frame(sample_class = factor(rep(c("ARP", "AT", "ODP"), 
                                                      c(5, 5, 5))))
rownames(annotation_row) = colnames(log_top_subtype)


pheatmap(t(log_top_subtype), annotation_row = annotation_row, border=TRUE,  border_color = "black")
