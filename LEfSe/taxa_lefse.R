# taxa_lefse.R

# Import library
library(tidyverse)
# Merge bracken file
## Import ARP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ARP/ARP1.S.bracken", 
                  header = TRUE ,sep = "\t", quote = "")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ARP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ARP/ARP",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t", quote = "") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ARP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
ARP_taxa<- tmp
## Import AT
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/AT/AT1.S.bracken", 
                  header = TRUE ,sep = "\t", quote = "")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "AT1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/AT/AT",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t", quote = "") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("AT", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
AT_taxa<- tmp
## Import ODP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ODP/ODP1.S.bracken", 
                  header = TRUE ,sep = "\t", quote = "")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ODP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ODP/ODP",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t", quote = "") # Import
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

# Change column name
colnames(merge_taxa) <- c("sample_type",rep("aeration_aera",10),rep("outdoor_aera",5))

# Adjust taxa name to debug in lefse
merge_taxa$sample_type <- gsub('\\[','',merge_taxa$sample_type)
merge_taxa$sample_type <- gsub('\\]','',merge_taxa$sample_type)
merge_taxa$sample_type <- gsub(' ','__',merge_taxa$sample_type) # Transform back later
merge_taxa$sample_type <- gsub('\\.','___',merge_taxa$sample_type) # Transform back later

# # Save LEfSE input file
# write.table(merge_taxa,"../../airborne_arg_uwtp_result/LEfSe/species_lefse_2sampletype_input.csv",
#             quote= FALSE, sep = "\t",row.names=FALSE)


###################################################################
###################################################################
# Preprocess of plotting
## 2 sample type
two_sampletype_lda <- read.table("../../airborne_arg_uwtp_result/LEfSe/Galaxy_LEfSe_species_2sampletype.txt",
                                 sep = "\t")
colnames(two_sampletype_lda) <- c("Taxa","LogMaxMean","Class","LDA","pValue")
# Transform the text back
two_sampletype_lda$Taxa <- gsub('___','.',two_sampletype_lda$Taxa) # Transform back later
two_sampletype_lda$Taxa <- gsub('__',' ',two_sampletype_lda$Taxa) # Transform back later

### Revise NTM five genera name ###
two_sampletype_lda$Taxa <- gsub("Mycobacteroides", "Mycobacterium", two_sampletype_lda$Taxa)
two_sampletype_lda$Taxa <- gsub("Mycolicibacillus", "Mycobacterium", two_sampletype_lda$Taxa)
two_sampletype_lda$Taxa <- gsub("Mycolicibacterium", "Mycobacterium", two_sampletype_lda$Taxa)
two_sampletype_lda$Taxa <- gsub("Mycolicibacter", "Mycobacterium", two_sampletype_lda$Taxa)

p <- two_sampletype_lda %>% 
  drop_na(LDA) %>%
  filter(LDA > 3.4) %>% 
  mutate(LDA = if_else(Class == "aeration_aera", -1*LDA, LDA),
         Taxa = fct_reorder(Taxa,LDA)) %>% 
  ggplot(aes(x=LDA, y=Taxa, fill = Class)) +
  geom_col(alpha = 0.6) +
  theme_bw() +
  labs(y = NULL, x = expression(LDA~Score~(log[10]))) +
  guides(fill=guide_legend(title="")) +
  scale_fill_discrete(name = 'Sample', 
                      labels = c(expression(Aeration~tank~+~Aeration~tank~PM[2.5]),
                                 expression(Outdoor~PM[2.5]))) +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position="top",
        legend.text = element_text(size = 14),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) #transparent legend bg


print(p)

# ggsave("species_LeFSe_two_sample_type.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/LEfSe",
#        width = 8, height = 5.5,
#        units = "in", bg='transparent') # save to png format
