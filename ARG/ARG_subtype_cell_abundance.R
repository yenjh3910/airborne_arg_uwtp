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
## Change the first column to row name
#row.names(arg_subtype) <- arg_subtype[,1]
#arg_subtype <- arg_subtype[,-1]




#Plot
## Convert table to matrix 
arg_subtype <- as.matrix(arg_subtype)
pheatmap(arg_subtype)
