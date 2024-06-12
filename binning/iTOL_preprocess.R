# iTOL_preprocess.R

# Import
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
bin_abundance <- read.table("../../airborne_arg_uwtp_result/metawrap_bin/bin_quant/bin_abundance_table.tab",
                            header = TRUE, sep="\t")
colnames(bin_abundance)[1] <- "MAG"
source("./bin_gtdbtk.R")
colnames(sep_bin_tax)[1] <- "MAG"
source("./bins_ARG_diamond.R")
source("./bins_MGE_diamond.R")
source("./bins_VF_diamond.R")