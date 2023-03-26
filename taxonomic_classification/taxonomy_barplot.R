# taxonomy_barplot.R

library(tidyverse)
library(vegan)

# Merge bracken file
## Import ARP
tmp <- read.table("D:/ARG_project/Shell/kraken2/ARP/ARP1.S.bracken", 
                                 header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ARP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("D:/ARG_project/Shell/kraken2/ARP/ARP",i,".S.bracken", sep = ""), 
  header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ARP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "Species"
ARP_species<- tmp
## Import AT
tmp <- read.table("D:/ARG_project/Shell/kraken2/AT/AT1.S.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "AT1"
for (i in 2:5) {
  tmp2 <- read.table(paste("D:/ARG_project/Shell/kraken2/AT/AT",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("AT", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "Species"
AT_species<- tmp
## Import ODP
tmp <- read.table("D:/ARG_project/Shell/kraken2/ODP/ODP1.S.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ODP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("D:/ARG_project/Shell/kraken2/ODP/ODP",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ODP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "Species"
ODP_species<- tmp
## Full join all sample type file
ARP_AT_species <- full_join(ARP_species, AT_species)
merge_species <- full_join(ARP_AT_species, ODP_species)
merge_species <- merge_species %>% replace(is.na(.), 0) # Replace NA with 0

# Filtering low abundance read
THRESHOLD <- 0.00005 # 0.005 % as filtering threshold
for (i in 2:16) {
  individual_threshold <- sum(merge_species[,i])*THRESHOLD
  merge_species[,i] <- replace(merge_species[,i], merge_species[,i] < individual_threshold, 0)
}
## Remove species if 0 in all sample
### First column to row name
row.names(merge_species) <- merge_species[,1]
merge_species <- merge_species[,-1]
### Romove all 0 row
merge_species <-  merge_species[rowSums(merge_species[])>0,]


# # Plot Rarefaction_curve
# merge_species_t=as.data.frame(t(merge_species))
# # Count the number of species
# S <- specnumber(merge_species_t)
# raremax <-min(rowSums(merge_species_t))
# # Rarefaction of the samples
# Srare <- rarefy(merge_species_t, raremax)
# Srare
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# abline(0, 1)
# pdf("Rarefaction_curve.pdf")
# rarecurve(merge_species_t, step =10000, sample = raremax, col = "blue", cex = 0.4)
# dev.off()


# Procrustes Analysis
## Import ARG subtype dataset
arg_subtype <- read.xlsx("D:/ARG_project/Shell/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
row.names(arg_subtype) <- arg_subtype[,1]
arg_subtype <- arg_subtype[,-1]


merge_species.dist<-vegdist(merge_species,method = "bray")
arg_subtype.dist<-vegdist(arg_subtype,method = "bray")

mds.species<-monoMDS(merge_species.dist)
mds.subtype<-monoMDS(arg_subtype.dist)

pro.species.subtype<-procrustes(mds.species,mds.subtype)
pro.species.subtype
protest(mds.species,mds.subtype)