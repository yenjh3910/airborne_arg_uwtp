# adjust_mge_database.R
# Structure file of MGE dastbase created by ourselves is different from the curated MGE DNA structure file 
# Therefore, we want to adjust our MGE database structure file

library(tidyverse)
library(compare)
# Import
MGE_AA_structure <- read.table("D:/ARG_project/airborne_arg_uwtp/MGE/MGE_structure/MGE_AA_structure.txt", 
                               sep ="\t", header = TRUE)
MGE_DNA_structure <- read.table("D:/ARG_project/airborne_arg_uwtp/MGE/MGE_structure/MGE_DNA_structure.txt", 
                               sep ="\t", header = TRUE)
# Arrange dataframe
colnames(MGE_AA_structure) <- "gene" # Change column name

# Print in MGE_AA_structure, but not in MGE_DNA_structure
for (i in MGE_AA_structure[,1]) {
  if (!(i %in% MGE_DNA_structure$gene)) {
    print(i)}
  } ### 8 genes were found here!!

# Print in MGE_DNA_structure, but not in MGE_AA_structure
for (i in MGE_DNA_structure[,1]) {
  if (!(i %in% MGE_AA_structure$gene)) {
    print(i)}
}  ### All genes in MGE_DNA_structure exist in MGE_AA_structure

# Add addition gene back to MGE_DNA_structure
gene <- c("398_tnpA*_EU287476.1", "345_[tnpA[_Tn3_like]KR822247.1",
          "2218_[tnpA[]KR822246.1", "1240_tnpA*_JX194160.1",
          "44_rep8_1_EP002*(pAM373)_NC002630", "755_IS26_KT334335.1",
          "2375_tnpA*_JX194160.1", "1534_tnpA*_JX194160.1")
subtype <- c("tnpA", "tnpA",
             "tnpA", "tnpA",
             "rep8", "IS26",
             "tnpA", "tnpA")
type <- c("transposase", "transposase",
          "transposase", "transposase",
          "plasmid", "IS26",
          "transposase", "transposase")
add_df <- data.frame(gene, subtype, type) 
curated_MGE_structure <- rbind(MGE_DNA_structure, add_df) # Created curated MGE database

# Save file (curated MGE database = MGE_AA_structure + MGE_DNA_structure)
# write.table(curated_MGE_structure, 
#             file = "D:/ARG_project/airborne_arg_uwtp/MGE/MGE_structure/MGE_curated_structure.txt",
#             sep = "\t", quote = FALSE, append = FALSE, row.names=  FALSE)