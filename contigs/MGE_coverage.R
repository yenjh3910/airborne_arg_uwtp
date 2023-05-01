# MGE_coverage.R

library(tidyverse)
library(openxlsx)
# MGE
## Import coverage file & bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_bowtie2/MGE/coverage", pattern = "*_MGE.sam.map.txt")
z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_bowtie2/MGE/coverage/", temp[1], sep = ''))
z$SampleID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_bowtie2/MGE/coverage/", temp[i], sep = ''))
  z2$SampleID <- temp[i]
  z <- rbind(z, z2)
}
colnames(z) <- c("ID","Avg_fold","Length","Ref_GC","Covered_percent","Covered_based",
                 "Plus_reads","Minus_reads","Read_GC","Median_fold","Std_Dev","SampleID") 
z$SampleID <- gsub("_MGE.sam.map.txt", "", z$SampleID)
## Import base number file
nbase <- read.xlsx("../../airborne_arg_uwtp_result/read_base_count.xlsx", sheet = 2)
## Join coverage file & base file
MGE_coverage <- full_join(z, nbase, by = "SampleID")
## Calculate coverage
MGE_coverage$coverage <- MGE_coverage$Avg_fold/(MGE_coverage$base/1000000000)
# Import Diamond file
source("./contigs_MGE_diamond.R")
# Join with MGE and taxonomy
MGE_coverage$ORF_SampleID <- paste(MGE_coverage$ID,
                                   MGE_coverage$SampleID,
                                   sep = "_")
MGE$ORF_SampleID <- paste(MGE$ORF,
                          MGE$SampleID,
                          sep = "_")
MGE_coverage <- left_join(MGE_coverage, MGE, by = "ORF_SampleID")
MGE_coverage <- MGE_coverage %>% select(ORF,Avg_fold,Length,Std_Dev,
                                        SampleID.x,coverage,subtype,type,
                                        contigs,`taxonomy ID`)