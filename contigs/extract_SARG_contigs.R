# extract_SARG_contigs.R

########## contigs.fa is too big to save in airborne_arg_uwtp_result #########
########## so import from another place ###########
library(seqinr)
source("./contigs_ARG_diamond.R")
temp <- list.files(path = "D:/Mirror/ARG_project/Shell/contigs_prodigal/", pattern = "*_contigs.nucl")
sample_list <- c("ARP1","ARP2","ARP3","ARP4","ARP5",
                 "AT1","AT2","AT3","AT4","AT5",
                 "ODP1","ODP2","ODP3","ODP4","ODP5")
for (i in 1:length(temp)){
  fastafile <- read.fasta(paste("D:/Mirror/ARG_project/Shell/contigs_prodigal/", temp[i], sep = ''), 
                seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
  SARG_subset <- SARG %>% filter(SampleID == sample_list[i])
  fastafile <- fastafile[names(fastafile) %in% SARG_subset$ORF]
  print(length(fastafile[names(fastafile)]))
  # write.fasta(sequences = fastafile, names = names(fastafile),
  #             file.out = paste("D:/Mirror/ARG_project/Shell/extract_contigs/SARG/",
  #                              gsub("contigs.nucl", "SARG.fa", temp[i]),
  #                              sep = ""),
  #             open = "w")
}
