# Rarefaction.R

library(tidyverse)
library(iNEXT)
library(vegan)
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
# First column as row name
rownames(merge_taxa) <- merge_taxa$taxa
merge_taxa <-  merge_taxa[,-1]
# Remove all 0 row
merge_taxa <- merge_taxa[rowSums(merge_taxa[])>0,]

# 
out <- iNEXT(merge_taxa, q=0, datatype="abundance")
p<-ggiNEXT(out, type=1)
# Choose color
library(RColorBrewer)
RColorBrewer::display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")

final <- p+scale_color_manual(values=c(rep("#80B1D3",5),rep("#FB8072",5),rep("#FDB462",5)))+
  scale_shape_manual(values=c(rep(17,5),rep(15,5),rep(18,5))) +
  scale_fill_manual(values=c(rep("#D9D9D9",15)))+
  xlab("Number of sequence") + ylab("Observed species richness") +
  theme_bw()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'))+
  geom_line(size=0.00001)
print(final)
# ggsave("rarefaction_curve.png", final,
#        path = "../../airborne_arg_uwtp_result/Figure/taxonomy",
#        width = 9, height = 4,
#        units = "in", bg='transparent') # save to png format

legend <- p+scale_color_manual(values=c(rep("#80B1D3",5),rep("#FB8072",5),rep("#FDB462",5)))+
  scale_shape_manual(values=c(rep(17,5),rep(15,5),rep(18,5))) +
  scale_fill_manual(values=c(rep("#D9D9D9",15)))+
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'))

# ggsave("rarefaction_legend.png", legend,
#        path = "../../airborne_arg_uwtp_result/Figure/taxonomy",
#        width = 9, height = 4,
#        units = "in", bg='transparent') # save to png format

# Speed up
library(parallel)
library(pbmcapply)
Max_CPU_Cores = detectCores()
Upper_Limit_CPU_Cores = 2*round((Max_CPU_Cores*0.8)/2)
# Parallel Rareification Function
# This is a working parallelized function of iNEXT. 5x faster than previously
parallel_rarefaction <- function(shuffled_data){
  out_df <- iNEXT(as.vector(shuffled_data), q=0, datatype="abundance")
  df <- fortify(out_df, type=1)
  return(df)
}
#iNEXT_output <- pbmclapply(merge_taxa, parallel_rarefaction, mc.cores = Upper_Limit_CPU_Cores)
