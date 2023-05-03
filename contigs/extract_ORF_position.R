# extract_ORF_position.R

########## contigs.nucl is too big to save in airborne_arg_uwtp_result #########
########## so import from another place ###########
library(seqinr)
temp <- list.files(path = "D:/Mirror/ARG_project/Shell/contigs_prodigal/", pattern = "*_contigs.nucl")
sample_list <- c("ARP1","ARP2","ARP3","ARP4","ARP5",
                 "AT1","AT2","AT3","AT4","AT5",
                 "ODP1","ODP2","ODP3","ODP4","ODP5")
# Create first dataframe to bind
fastafile <- read.fasta(paste("D:/Mirror/ARG_project/Shell/contigs_prodigal/", temp[1], sep = ''), 
                        seqtype = "DNA",as.string = TRUE, 
                        set.attributes = FALSE, whole.header =  TRUE)
FastaHeader <- as.data.frame(names(fastafile))
colnames(FastaHeader) <- "header"
FastaHeader <- FastaHeader %>%
  separate(header, sep="\\#",into=c("ORF","start","end","direction","detail"))
FastaHeader$ORF_SampleID <- paste(FastaHeader$ORF,
                                  sample_list[1],
                                  sep = "_")
FastaHeader <- FastaHeader %>% select(start,end,direction,ORF_SampleID)
z <- FastaHeader
# Bind each sample dataframe
for (i in 2:length(temp)){
  fastafile <- read.fasta(paste("D:/Mirror/ARG_project/Shell/contigs_prodigal/", temp[i], sep = ''), 
                          seqtype = "DNA",as.string = TRUE, 
                          set.attributes = FALSE, whole.header =  TRUE)
  FastaHeader <- as.data.frame(names(fastafile))
  colnames(FastaHeader) <- "header"
  FastaHeader <- FastaHeader %>%
    separate(header, sep="\\#",into=c("ORF","start","end","direction","detail"))
  FastaHeader$ORF_SampleID <- paste(FastaHeader$ORF,
                                    sample_list[i],
                                    sep = "_")
  FastaHeader <- FastaHeader %>% select(start,end,direction,ORF_SampleID)
  z <- rbind(z, FastaHeader)
}
ORF_position <- z

# write_csv(ORF_position, 
#           "../../airborne_arg_uwtp_result/contigs_ORF_position/contigs_ORF_position.csv")