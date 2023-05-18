# bins_ARG_diamond

library(tidyverse)
library(openxlsx)
# SARG
## Import file & bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/bins_diamond/SARG/", pattern = "*.SARG.dmnd")
z<-read.table(paste("../../airborne_arg_uwtp_result/bins_diamond/SARG/", temp[1], sep = ''))
z$BinID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/bins_diamond/SARG/", temp[i], sep = ''))
  z2$BinID <- temp[i]
  z <- rbind(z, z2)
}
colnames(z) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                 "qstart","qend","sstart","send","evalue","bitscore","BinID") 
z$BinID <- gsub(".SARG.dmnd", "", z$BinID)
## Join with SARG database
multi_component_structure <- read.xlsx("../ARG/multi-component_structure.xlsx",sheet = 1)
single_component_structure <- read.xlsx("../ARG/single-component_structure.xlsx",sheet = 1)
two_component_structure <- read.xlsx("../ARG/two-component_structure.xlsx",sheet = 1)
component_structure <- rbind(multi_component_structure,single_component_structure, two_component_structure)
colnames(component_structure)[1] <- "sseqid"
z <- left_join(z, component_structure, by="sseqid")
## Add contigs column
z$ORF <- z$qseqid
z <- z %>% separate(qseqid, c("tmp1","tmp2","tmp3"), sep = "_")
z$contigs <- paste(z$tmp1,z$tmp2,sep = "_")
z <- z %>% select(!(tmp1)) %>% select(!(tmp2)) %>% select(!(tmp3))
## Filter by identity & evalue
SARG <- z %>% filter(pident >= 70) %>% filter(evalue <= 1e-10) %>% filter(evalue <= 1e-10)





# ## Import contigs kraken file
# source("./contigs_kraken2.R")
# contigs_kraken <- contigs_kraken %>% select(contigs, `taxonomy ID`,SampleID)
# contigs_kraken$contigs_SampleID <- paste(contigs_kraken$contigs,
#                                          contigs_kraken$SampleID,
#                                          sep = "_")
# contigs_kraken <- contigs_kraken %>% select(contigs_SampleID, `taxonomy ID`)
# # Join diamond & kraken
# SARG$contigs_SampleID <- paste(SARG$contigs,
#                                SARG$SampleID,
#                                sep = "_")
# SARG <- left_join(SARG, contigs_kraken, by = "contigs_SampleID")
# SARG <- SARG %>% select(SampleID,subtype,type,ORF,contigs,`taxonomy ID`)
