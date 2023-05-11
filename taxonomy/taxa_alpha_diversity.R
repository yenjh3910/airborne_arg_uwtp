# taxa_alpha_dicersity.R

library(vegan)
library(tidyverse)
library(openxlsx)

# Merge bracken file
## Import ARP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ARP/ARP1.S.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ARP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ARP/ARP",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ARP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
ARP_taxa<- tmp
## Import AT
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/AT/AT1.S.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "AT1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/AT/AT",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("AT", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
AT_taxa<- tmp
## Import ODP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ODP/ODP1.S.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ODP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ODP/ODP",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ODP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
ODP_taxa<- tmp
## Full join all sample type file
ARP_AT_taxa <- full_join(ARP_taxa, AT_taxa)
merge_taxa <- full_join(ARP_AT_taxa, ODP_taxa)
merge_taxa <- merge_taxa %>% replace(is.na(.), 0) # Replace NA with 0
# Filtering low abundance read
THRESHOLD <- 0.00005 # 0.005 % as filtering threshold
for (i in 2:16) {
  individual_threshold <- sum(merge_taxa[,i])*THRESHOLD
  merge_taxa[,i] <- replace(merge_taxa[,i], merge_taxa[,i] < individual_threshold, 0)
}
# Fitst column as row name
rownames(merge_taxa) <- merge_taxa$taxa
merge_taxa <-  merge_taxa[,-1]
# Remove all 0 row
merge_taxa <- merge_taxa[rowSums(merge_taxa[])>0,]
# Transform row name back to column one
merge_taxa$taxa <- rownames(merge_taxa)
# Gather transformation
gat_merge_taxa <- merge_taxa %>% gather(key = "sample", value = "abundance",ARP1:ODP5)
# Add sample type column
gat_merge_taxa$sample_type <- gat_merge_taxa$sample
gat_merge_taxa$sample_type <- gsub("1|2|3|4|5","",gat_merge_taxa$sample)


richness <- function(x){
  
  # r <- sum(x > 0)
  # return(r)
  
  sum(x>0)
}

shannon <- function(x){
  
  rabund <- x[x>0]/sum(x)
  -sum(rabund * log(rabund))
  
}

simpson <- function(x){
  
  n <- sum(x)
  
  # sum(x * (x-1) / (n * (n-1)))
  1 - sum((x/n)^2)
}

gat_merge_taxa %>%
  group_by(sample_type) %>%
  summarize(sobs = richness(value),
            shannon = shannon(value),
            simpson = simpson(value),
            invsimpson = 1/simpson,
            n = sum(value)) %>%
  pivot_longer(cols=c(sobs, shannon, invsimpson, simpson),
               names_to="metric") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")


gat_merge_taxa %>%
  group_by(sample) %>%
  summarize(sobs = specnumber(abundance),
            shannon = diversity(abundance, index="shannon"),
            simpson = diversity(abundance, index="simpson"),
            invsimpson = 1/simpson,
            n = sum(abundance))

  pivot_longer(cols=c(sobs, shannon, invsimpson, simpson),
               names_to="metric") %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~metric, nrow=4, scales="free_y")