# bin_gtdbtk.R

library(tidyverse)
bin_tax <- read.table("../../airborne_arg_uwtp_result/bin_gtdbtk/gtdbtk.bac120.summary.tsv",
           sep = "\t", header = TRUE)
sep_bin_tax <- bin_tax %>% select(user_genome,classification) %>% 
                separate(classification, sep=";",
                into=c("Domain","Phylum","Class","Order","Family","Genus","Species"))