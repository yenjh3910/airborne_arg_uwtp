# ARG_subtype_cell_abundance.R
# Calculate and visualize ARG subtype/cell in each sample

# Import library
library(tidyverse)
library(stringr)
library(tibble)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(scales)

# Read ARG_subtype file (Something wrong with read.table, so read by )
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)

# Preview raw data
gather_arg_subtype <- gather(arg_subtype, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5) # Transform to gather format
gather_arg_subtype_order <- gather_arg_subtype[order(gather_arg_subtype$copy_per_cell, 
                                                     decreasing = T),] # Order by abundance
top_subtype_list<- unique(gather_arg_subtype_order[,1])[1:75] #### Select top subtype (Choose)
## View top ARG percentage
top_subtype <- filter(arg_subtype, subtype %in% top_subtype_list) # Assign top subtype dataframe
top_arg_sum <- 0 # Initialization
total_arg_sum <- 0 # Initialization
arg_select<- paste("Top ",length(top_subtype_list),sep = "")
### Print statistic result 
for (i in c(2:16)) {
  print(paste(arg_select,' ARGs proportion in ',colnames(top_subtype)[i],': ',
              label_percent(accuracy = 0.01)(sum(top_subtype[,i])/sum(arg_subtype[,i])), 
              sep = ""))
  top_arg_sum <- sum(top_subtype[,i]) + top_arg_sum
  total_arg_sum <- sum(arg_subtype[,i]) + total_arg_sum
  } # Print top ARG percentage in each sample
total_top_proportion <- label_percent(accuracy = 0.01)(top_arg_sum/total_arg_sum)
print(paste(arg_select,'ARGs proportion in all samples:', 
            total_top_proportion), sep="") # Print top ARG percentage in all sample

## Adjust "_" in dataframe
top_subtype <- top_subtype %>% 
  separate(subtype, c("type","subtype"), 
           sep = "__") # Split subtype by "_"
top_subtype$type <- gsub("_", " ", top_subtype$type) # Change _ to " " (e.g. beta_lactam to beta lactam)

## Create subtype-type table
type_subtype <- top_subtype %>% select(type, subtype)
row.names(type_subtype) <- type_subtype$subtype
type_subtype <- type_subtype %>% select(!subtype)
colnames(type_subtype) <- "ARG" # Change column name

top_subtype <- top_subtype[,-1] # Remove ARG type

# Change the first column to row name
row.names(top_subtype) <- top_subtype[,1]
top_subtype <- top_subtype[,-1]
row.names(top_subtype) <- gsub("Bifidobacteria intrinsic ileS conferring resistance to mupirocin", 
                                  "Bifidobacteria intrinsic ileS", 
                               row.names(top_subtype)) # Adjust subtype name

# Usage pheatmap
top_subtype <- as.matrix(top_subtype) # Convert table to matrix 
log_top_subtype <- log10(top_subtype) # log10 conversion 
log_top_subtype[log_top_subtype == -Inf] <- NA

## Create annotation label (Note annotation_row was already created as type_subtype)
annotation_row = data.frame(Sample = factor(rep(c("ARP", "AT", "ODP"), 
                                                      c(5, 5, 5))))
rownames(annotation_row) = colnames(log_top_subtype)
annotation_col <- type_subtype
annotation_col$ARG <- str_to_sentence(annotation_col$ARG) # Covert first letter to uppercase
annotation_col$ARG <- gsub("Macrolide-lincosamide-streptogramin", "MLS", annotation_col$ARG)
row.names(annotation_col) <- gsub("Bifidobacteria intrinsic ileS conferring resistance to mupirocin", 
                           "Bifidobacteria intrinsic ileS", row.names(annotation_col))

# Select color
display.brewer.all()
brewer.pal(3, "Set1") # However, I did not select these color finally
library("ggsci")
ann_colors = list(Sample = c(AT = "#FC4E07", ARP = "#00AFBB", ODP = "#E7B800"))

## Plot
phtmap <- pheatmap(t(log_top_subtype), 
         annotation_row = annotation_row, annotation_col = annotation_col,
         cluster_cols = FALSE, clustering_distance_rows = "euclidean")

### Rotate the cluster
library(seriation)
library(dendextend)
row_dend <- phtmap[[1]]
row_dend <- dendextend::rotate(row_dend, order = c("ODP5","ODP1","ODP3","ODP2","ODP4",
                                                   "ARP4","ARP2","ARP3","ARP1","ARP5",
                                                   "AT2","AT5","AT1","AT3","AT4") )
p <- pheatmap(t(log_top_subtype), 
         annotation_row = annotation_row, annotation_col = annotation_col,
         cluster_cols = FALSE, annotation_colors = ann_colors,
         clustering_distance_rows = "euclidean", 
         fontsize = 12, fontsize_row = 12, fontsize_col = 10,
         cellwidth = 12, cellheight = 17, bg = "transparent", 
         cluster_rows=as.hclust(row_dend))

# ggsave("ARG_subtype_cell.png", p, 
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 17, height = 8, 
#        units = "in", bg='transparent') # save to png format