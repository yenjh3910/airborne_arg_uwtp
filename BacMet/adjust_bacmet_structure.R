library(tidyverse)

# Modify structure file
bacmet_official_structure <- read.table("D:/ARG_project/airborne_arg_uwtp/BacMet/BacMet_structure/BacMet_official_structure.txt", 
                              sep = "\t", header = TRUE)
bacMet_structure <- read.table("D:/ARG_project/airborne_arg_uwtp/BacMet/BacMet_structure/BacMet_structure.txt", 
                                       sep = "\t", header = TRUE)
bacMet_structure <- bacMet_structure %>% 
  separate(level1, c('BacMet_ID', 'Gene_name','Unknown','Accession','Complete_accession'), sep = "|")





# Modify fasta
library(Biostrings)
fasta <- readDNAStringSet("D:/ARG_project/airborne_arg_uwtp/BacMet/BacMet_structure/BacMet_experimentally_confirmed.fasta")

# Modify the header using regular expressions
newHeaders <- gsub("(>.*?\\|abeM).*", "", names(fasta))

# assign the new headers to the object
names(fasta) <- newHeaders

# write out the modified fasta file
writeXStringSet(fasta, "modified.fasta", format="fasta")