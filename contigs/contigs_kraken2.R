# contigs_kraken2.R

library(tidyverse)
# Import file and bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_kraken2/", pattern = "*_contigs_kraken2.output")
z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_kraken2/", temp[1], sep = ''),
              sep = "\t",quote = "")
z$SampleID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_kraken2/", temp[i], sep = ''),
                  sep = "\t",quote = "")
  z2$SampleID <- temp[i]
  z <- rbind(z, z2)
}
colnames(z) <- c("classified","contigs","taxonomy ID","length","LCA mapping","SampleID")
z$SampleID <- gsub("_contigs_kraken2.output", "", z$SampleID)
contigs_kraken <- z