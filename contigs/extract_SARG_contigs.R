# extract_SARG_contigs.R

########## contigs.fa is too big to save in airborne_arg_uwtp_result #########
########## so import from another place ###########
library(seqinr)
source("./contigs_ARG_diamond.R")
temp <- list.files(path = "D:/Mirror/ARG_project/Shell/contigs_prodigal/", pattern = "*_contigs.nucl")
for (i in 1:length(temp)){
  fastafile <- read.fasta(paste("D:/Mirror/ARG_project/Shell/contigs_prodigal/", temp[i], sep = ''), 
                seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
  fastafile <- fastafile[names(fastafile) %in% SARG$ORF]
  write.fasta(sequences = fastafile, names = names(fastafile), 
              file.out = paste("D:/Mirror/ARG_project/Shell/extract_contigs/",
                               gsub("contigs.nucl", "SARG.fa", temp[i]),
                               sep = ""),
              open = "w")
}
