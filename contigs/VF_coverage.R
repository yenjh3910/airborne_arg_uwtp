# VF_coverage.R

library(tidyverse)
library(openxlsx)
# VF
## Import coverage file & bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_bowtie2/VF/coverage", pattern = "*_VF.sam.map.txt")
z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_bowtie2/VF/coverage/", temp[1], sep = ''))
z$SampleID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_bowtie2/VF/coverage/", temp[i], sep = ''))
  z2$SampleID <- temp[i]
  z <- rbind(z, z2)
}
colnames(z) <- c("ID","Avg_fold","Length","Ref_GC","Covered_percent","Covered_based",
                 "Plus_reads","Minus_reads","Read_GC","Median_fold","Std_Dev","SampleID") 
z$SampleID <- gsub("_VF.sam.map.txt", "", z$SampleID)
## Import base number file
nbase <- read.xlsx("../../airborne_arg_uwtp_result/read_base_count.xlsx", sheet = 2)
## Join coverage file & base file
VF_coverage <- full_join(z, nbase, by = "SampleID")
## Calculate coverage
VF_coverage$coverage <- VF_coverage$Avg_fold/(VF_coverage$base/1000000000)
# Import Diamond file
source("./contigs_VF_diamond.R")
# Join with VF and taxonomy
VF_coverage$ORF_SampleID <- paste(VF_coverage$ID,
                                  VF_coverage$SampleID,
                                  sep = "_")
VF$ORF_SampleID <- paste(VF$ORF,
                         VF$SampleID,
                         sep = "_")
VF_coverage <- left_join(VF_coverage, VF, by = "ORF_SampleID")
VF_coverage <- VF_coverage %>% select(ORF,Avg_fold,Length,Std_Dev,
                                      SampleID.x,coverage,gene,
                                      contigs,`taxonomy ID`)