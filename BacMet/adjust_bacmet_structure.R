library(tidyverse)

# Modify structure file
bacmet_official_structure <- read.table("./BacMet_structure/BacMet_official_structure.txt", 
                              sep = "\t", header = TRUE)
bacmet_curated_structure <- bacmet_official_structure %>% select(BacMet_ID, Gene_name, Compound)
colnames(bacmet_curated_structure) <- c("gene","subtype","type")
# Filter metal only
metal_list <- c("Aluminium (Al)","Antimony (Sb)","Arsenic (As)","Bismuth (Bi)",
                "Cadmium (Cd)","Chromium (Cr)","Cobalt (Co)","Copper (Cu)",
                "Gallium (Ga)","Gold (Au)","Iron (Fe)","Lead (Pb)",
                "Magnesium (Mg)","Manganese (Mn)","Mercury (Hg)","Molybdenum (Mo)",
                "Nickel (Ni)","Selenium (Se)","Silver (Ag)","Tellurium (Te)",
                "Vanadium (V)","Zinc (Zn)")

metal_only_structure <- bacmet_curated_structure %>% 
  filter(grepl("Aluminium|Antimony|Arsenic|Bismuth|
               Cadmium|Chromium|Cobalt|Copper|
               Gallium|Gold|Iron|Lead|
               Magnesium|Manganese|Mercury|Molybdenum|
               Nickel|Selenium|Silver|Tellurium|
               Vanadium|Zinc",type)) # Filter to metal only dataframe

metal_only_structure$type <- gsub(".*,.*" , "Multi-metal" , 
                                  metal_only_structure$type)  # Add multi-metal tag

# write.table(metal_only_structure, "./BacMet_structure/metal_only_structure.txt",
#             sep = "\t", quote = FALSE, append = FALSE, row.names=  FALSE)

######################################################################################

# Modify fasta
library(Biostrings)
library(stringr)
fasta <- readDNAStringSet("./BacMet_structure/BacMet_experimentally_confirmed.fasta")
names(fasta)
# Modify the header using regular expressions
names(fasta) <- str_extract_all(names(fasta),"^BAC....")
# Filter metal resistance gene
fasta<- fasta[names(fasta) %in% metal_only_structure$gene]
# write out the modified fasta file
#writeXStringSet(fasta, "./BacMet_structure/BacMet_exp_metal.fasta", format="fasta")
