# contigs_VF_diamond.R

library(tidyverse)
library(seqinr)
# VF
## Import file & bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_diamond/VF/", pattern = "*_contigs.VF.dmnd")
z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_diamond/VF/", temp[1], sep = ''))
z$SampleID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_diamond/VF/", temp[i], sep = ''))
  z2$SampleID <- temp[i]
  z <- rbind(z, z2)
}
colnames(z) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                 "qstart","qend","sstart","send","evalue","bitscore","SampleID") 
z$SampleID <- gsub("_contigs.VF.dmnd", "", z$SampleID)
## Join with VF database
fastafile <- read.fasta("./VFDB_setB_pro.fa",seqtype = "AA",
                   as.string = TRUE, set.attributes = FALSE,
                   whole.header = TRUE)
fastafile <- names(fastafile)
fastafile <- as.data.frame(fastafile)
VF_db <- fastafile %>% separate(fastafile, c("sseqid","gene"), sep = "^\\S*\\K\\s+")
z <- left_join(z, VF_db, by="sseqid")
## Add contigs column
z$ORF <- z$qseqid
z <- z %>% separate(qseqid, c("tmp1","tmp2","tmp3"), sep = "_")
z$contigs <- paste(z$tmp1,z$tmp2,sep = "_")
z <- z %>% select(!(tmp1)) %>% select(!(tmp2)) %>% select(!(tmp3))
## Filter by identity & evalue
VF <- z %>% filter(pident >= 70) %>% filter(evalue <= 1e-10) %>% filter(evalue <= 1e-10)
## Import contigs kraken file
source("./contigs_kraken2.R")
contigs_kraken <- contigs_kraken %>% select(contigs, `taxonomy ID`,SampleID)
contigs_kraken$contigs_SampleID <- paste(contigs_kraken$contigs,
                                         contigs_kraken$SampleID,
                                         sep = "_")
contigs_kraken <- contigs_kraken %>% select(contigs_SampleID, `taxonomy ID`)
# Join diamond & kraken
VF$contigs_SampleID <- paste(VF$contigs,
                             VF$SampleID,
                             sep = "_")
VF <- left_join(VF, contigs_kraken, by = "contigs_SampleID")
VF <- VF %>% select(SampleID,gene,ORF,contigs,`taxonomy ID`)
