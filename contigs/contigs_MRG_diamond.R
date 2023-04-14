# contigs_MRG_diamond

library(tidyverse)
# MRG
## Import file & bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_diamond/MRG/", pattern = "*_contigs.MRG.dmnd")
z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_diamond/MRG/", temp[1], sep = ''))
z$SampleID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_diamond/MRG/", temp[i], sep = ''))
  z2$SampleID <- temp[i]
  z <- rbind(z, z2)
}
######### Some sample do not align any MRG ###########