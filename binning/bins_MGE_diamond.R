# bins_MGE_diamond

library(tidyverse)
library(openxlsx)
# MGE
## Import file & bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/bins_diamond/MGE/", pattern = "*.MGE.dmnd")
z<-read.table(paste("../../airborne_arg_uwtp_result/bins_diamond/MGE/", temp[1], sep = ''))
z$BinID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/bins_diamond/MGE/", temp[i], sep = ''))
  z2$BinID <- temp[i]
  z <- rbind(z, z2)
}
colnames(z) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                 "qstart","qend","sstart","send","evalue","bitscore","BinID") 
z$BinID <- gsub(".MGE.dmnd", "", z$BinID)
## Join with MGE database
MGE_db <- read.table("../MGE/MGE_structure/MGE_curated_structure.txt", header = TRUE)
colnames(MGE_db)[1] <- "sseqid"
z <- left_join(z, MGE_db, by="sseqid")
## Add contigs column
z$ORF <- z$qseqid
z <- z %>% separate(qseqid, c("tmp1","tmp2","tmp3"), sep = "_")
z$contigs <- paste(z$tmp1,z$tmp2,sep = "_")
z <- z %>% select(!(tmp1)) %>% select(!(tmp2)) %>% select(!(tmp3))
## Filter by identity & evalue
MGE <- z %>% filter(pident >= 70) %>% filter(evalue <= 1e-10) %>% filter(evalue <= 1e-10)




# ## Import contigs kraken file
# source("./contigs_kraken2.R")
# contigs_kraken <- contigs_kraken %>% select(contigs, `taxonomy ID`,SampleID)
# contigs_kraken$contigs_SampleID <- paste(contigs_kraken$contigs,
#                                          contigs_kraken$SampleID,
#                                          sep = "_")
# contigs_kraken <- contigs_kraken %>% select(contigs_SampleID, `taxonomy ID`)
# # Join diamond & kraken
# MGE$contigs_SampleID <- paste(MGE$contigs,
#                               MGE$SampleID,
#                               sep = "_")
# MGE <- left_join(MGE, contigs_kraken, by = "contigs_SampleID")
# MGE <- MGE %>% select(SampleID,subtype,type,ORF,contigs,`taxonomy ID`)