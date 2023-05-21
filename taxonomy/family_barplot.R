# family_barplot.R

library(tidyverse)
library(vegan)
library(openxlsx)

# Merge bracken file
## Import ARP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ARP/ARP1.F.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ARP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ARP/ARP",i,".F.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ARP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
ARP_taxa<- tmp
## Import AT
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/AT/AT1.F.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "AT1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/AT/AT",i,".F.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("AT", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
AT_taxa<- tmp
## Import ODP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ODP/ODP1.F.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ODP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ODP/ODP",i,".F.bracken", sep = ""), 
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
## Transform to percentage
for (i in 1:ncol(merge_taxa)) {
  sum_read<- sum(merge_taxa[,i])
  for (j in 1:nrow(merge_taxa))
    merge_taxa[j,i]<- merge_taxa[j,i]/sum_read
}



# Remove all 0 row
merge_taxa <- merge_taxa[rowSums(merge_taxa[])>0,]
# Transform row name back to column one
merge_taxa$taxa <- rownames(merge_taxa)
# Gather transformation
gat_merge_taxa <- merge_taxa %>% gather(key = "sample", value = "abundance",ARP1:ODP5)
# Unique taxa family
taxa_family <- unique(
  gat_merge_taxa[
    order(gat_merge_taxa$abundance,
          decreasing = T),
  ]$taxa)

# Statistic for all taxa
all_taxa <- gat_merge_taxa
all_taxa$sample_type <- gsub("1|2|3|4|5","",all_taxa$sample)
all_taxa <- all_taxa %>% group_by(taxa,sample_type) %>% 
                         mutate(mean_percent = mean(abundance*100)) %>%
                         mutate(sd_percent = sd(abundance*100)) %>%
                         select(taxa,sample_type,mean_percent,sd_percent) %>% unique()
spread_all_taxa <- all_taxa %>% select(!(sd_percent)) %>% spread(key = "sample_type", value = "mean_percent")

# family taxa as their highest abundance
merge_taxa <- merge_taxa[,-16]
merge_taxa <- merge_taxa %>% 
  arrange(factor(rownames(merge_taxa), levels = taxa_family))
# Calculate others by summing  minimum taxa
## Select other taxa
other_taxa <- colSums(
  merge_taxa[
    12:length(merge_taxa[,1]),])
merge_taxa  <- merge_taxa [-(12:length(merge_taxa[,1])),] # delete minimum taxa
merge_taxa <- rbind(merge_taxa, other_taxa)
rownames(merge_taxa)[length(rownames(merge_taxa))] <- "Others" 
merge_taxa <- tibble::rownames_to_column(merge_taxa, "taxa") # Convert row name back to first column
# Gather transformation
gat_merge_taxa <- merge_taxa %>% gather(key = "sample", value = "abundance",ARP1:ODP5)
# Add sample type column
gat_merge_taxa$sample_type <- gat_merge_taxa$sample
gat_merge_taxa$sample_type <- gsub("1|2|3|4|5","",gat_merge_taxa$sample)
# Change facet title
gat_merge_taxa$sample_type <- factor(gat_merge_taxa$sample_type, levels = c("AT", "ARP", "ODP"), 
                                     labels = c(expression(Aeration~tank), 
                                                expression(Aeration~tank~PM[2.5]), 
                                                expression(Outdoor~PM[2.5])
                                     ))
## family taxa
gat_merge_taxa$taxa <- factor(gat_merge_taxa$taxa , 
                              levels = c(taxa_family[1:11], 
                                         "Others"))


# Plot
p<- ggplot(gat_merge_taxa, aes(x = sample, y = abundance, fill = taxa)) + 
  geom_bar(stat="identity", position = "fill") + 
  facet_wrap(~sample_type, scales = "free_x", labeller = label_parsed) + 
  theme_bw() + 
  xlab("") + ylab("Relative abundance (%)") + 
  guides(fill=guide_legend(title="Family")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size=10.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) + #transparent legend bg
        scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                                   "#FDB462","#B3DE69","#FCCDE5", "#FFED6F","#BC80BD",
                                   "#CCEBC5","#D9D9D9"))
  
print(p)

# ggsave("family.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/taxonomy",
#        width = 7.75, height = 5,
#        units = "in", bg='transparent') # save to png format
