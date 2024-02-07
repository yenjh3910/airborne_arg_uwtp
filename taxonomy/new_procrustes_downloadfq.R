# ARG_species_procrustes_downloadfq.R

library(tidyverse)
library(vegan)
library(openxlsx)

# Merge bracken file
## Import ARP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ARP/ARP1.S.bracken", 
                                 header = TRUE ,sep = "\t",quote = "")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ARP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ARP/ARP",i,".S.bracken", sep = ""), 
  header = TRUE ,sep = "\t",quote = "") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ARP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "Species"
ARP_species<- tmp
## Import AT
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/AT/AT1.S.bracken", 
                  header = TRUE,sep = "\t",quote = "")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "AT1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/AT/AT",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t",quote = "") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("AT", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "Species"
AT_species<- tmp
## Import ODP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ODP/ODP1.S.bracken", 
                  header = TRUE ,sep = "\t",quote = "")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ODP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ODP/ODP",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t",quote = "") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ODP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "Species"
ODP_species<- tmp

## Import nineyear
tmp <- read.table("../../airborne_arg_uwtp_result/download_fq/nineyear_args_oap/SRR6747711.S.bracken",
                  header = TRUE ,sep = "\t",quote = "")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "SRR6747711"
for (i in 12:19) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/download_fq/nineyear_args_oap/SRR67477",i,".S.bracken", sep = ""),
                     header = TRUE ,sep = "\t",quote = "",fill = TRUE) # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("SRR67477", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
for (i in 21:25) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/download_fq/nineyear_args_oap/SRR67477",i,".S.bracken", sep = ""),
                     header = TRUE ,sep = "\t",quote = "",fill = TRUE) # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("SRR67477", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "Species"
nineyear_species<- tmp

# ## Import SRP107015
# tmp <- read.table("../../airborne_arg_uwtp_result/download_fq/SRP107015/SRR5571011.S.bracken",
#                   header = TRUE ,sep = "\t",quote = "")
# tmp <- tmp %>% select(name, new_est_reads)
# colnames(tmp)[2] <- "SRR5571011"
# for (i in 8:9) {
#   tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/download_fq/SRP107015/SRR557100",i,".S.bracken", sep = ""),
#                      header = TRUE ,sep = "\t",quote = "",fill = TRUE) # Import
#   tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
#   colnames(tmp2)[2] <- paste("SRR557100", i, sep = "") # Change column headers to sample name
#   tmp <- full_join(tmp ,tmp2) # Merge each sample file
#   tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
# }
# colnames(tmp)[1] <- "Species"
# SRP107015_species<- tmp
# 
# ## Import PRJEB26809
# tmp <- read.table("../../airborne_arg_uwtp_result/download_fq/PRJEB26809/ERR2586218.S.bracken",
#                   header = TRUE ,sep = "\t",quote = "")
# tmp <- tmp %>% select(name, new_est_reads)
# colnames(tmp)[2] <- "ERR2586218"
# for (i in 19:20) {
#   tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/download_fq/PRJEB26809/ERR25862",i,".S.bracken", sep = ""),
#                      header = TRUE ,sep = "\t",quote = "",fill = TRUE) # Import
#   tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
#   colnames(tmp2)[2] <- paste("ERR25862", i, sep = "") # Change column headers to sample name
#   tmp <- full_join(tmp ,tmp2) # Merge each sample file
#   tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
# }
# colnames(tmp)[1] <- "Species"
# PRJEB26809_species<- tmp

## Import PRJNA178295
tmp <- read.table("../../airborne_arg_uwtp_result/download_fq/PRJNA178295/SRR609293.S.bracken",
                  header = TRUE ,sep = "\t",quote = "")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "SRR609293"
tmp2 <- read.table("../../airborne_arg_uwtp_result/download_fq/PRJNA178295/SRR609435.S.bracken",
                  header = TRUE ,sep = "\t",quote = "")
tmp2 <- tmp2 %>% select(name, new_est_reads)
colnames(tmp2)[2] <- "SRR609435"
tmp <- full_join(tmp ,tmp2)
tmp2 <- read.table("../../airborne_arg_uwtp_result/download_fq/PRJNA178295/SRR609449.S.bracken",
                   header = TRUE ,sep = "\t",quote = "")
tmp2 <- tmp2 %>% select(name, new_est_reads)
colnames(tmp2)[2] <- "SRR609449"
tmp <- full_join(tmp ,tmp2)
colnames(tmp)[1] <- "Species"
PRJNA178295_species<- tmp

## Full join all sample type file
ARP_AT_species <- full_join(ARP_species, AT_species)
merge_species <- full_join(ARP_AT_species, ODP_species)
merge_species <- full_join(merge_species, nineyear_species)
# merge_species <- full_join(merge_species, SRP107015_species)
# merge_species <- full_join(merge_species, PRJEB26809_species)
merge_species <- full_join(merge_species, PRJNA178295_species)
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

# ## Transform to percentage
# for (i in 1:ncol(merge_species)) {
#   sum_read<- sum(merge_species[,i])
#   for (j in 1:nrow(merge_species))
#   merge_species[j,i]<- merge_species[j,i]/sum_read
#   }



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
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
### nineyear
download_fq <- read.table("../../airborne_arg_uwtp_result/download_fq/nineyear_args_oap/normalized_cell.subtype.txt",
                          header = TRUE, sep = "\t", quote = "")
download_fq <- download_fq %>% select(!(SRR6747720)) ### Remove because fastq has error
arg_subtype <- full_join(arg_subtype, download_fq)

# ### SRP107015
# download_fq <- read.table("../../airborne_arg_uwtp_result/download_fq/SRP107015/normalized_cell.subtype.txt",
#                           header = TRUE, sep = "\t", quote = "")
# arg_subtype <- full_join(arg_subtype, download_fq)
# 
# ### PRJEB26809
# download_fq <- read.table("../../airborne_arg_uwtp_result/download_fq/PRJEB26809/normalized_cell.subtype.txt",
#                           header = TRUE, sep = "\t", quote = "")
# arg_subtype <- full_join(arg_subtype, download_fq)

### PRJEB26809
download_fq <- read.table("../../airborne_arg_uwtp_result/download_fq/PRJNA178295/normalized_cell.subtype.txt",
                          header = TRUE, sep = "\t", quote = "")
arg_subtype <- full_join(arg_subtype, download_fq)

row.names(arg_subtype) <- arg_subtype[,1]
arg_subtype <- arg_subtype[,-1]
arg_subtype[is.na(arg_subtype)] <- 0
## Transform dataframe
merge_species <-as.data.frame(t(merge_species))
arg_subtype <-as.data.frame(t(arg_subtype))
## Normalization
arg_subtype <- decostand(arg_subtype, method = 'hellinger')
merge_species <- decostand(merge_species, method = 'hellinger')
taxa_bray<-vegdist(merge_species, method="bray")
arg_bray<-vegdist(arg_subtype, method="bray")
## Choose pcoa,pca,and nmds....
pcoa1 = cmdscale(taxa_bray, eig=TRUE)
pcoa2 = cmdscale(arg_bray, eig=TRUE)
pro.g.s<-procrustes(pcoa1,pcoa2,symmetric = F)
## Check statistic result
summary(pro.g.s)
plot(pro.g.s, kind = 1,type="text")
plot(pro.g.s, kind = 2)
## Check statistic result by protest
prot <- protest(X = pcoa1, Y = pcoa2, permutations = 999)
prot
names(prot)
prot$signif  # p value
prot$ss  # M2

Y<-cbind(data.frame(pro.g.s$Yrot),data.frame(pro.g.s$X))
X<-data.frame(pro.g.s$rotation)
Y$sample_type<-rownames(Y)
Y$sample_type <- gsub("SRR67477","nineyear",Y$sample_type)
# Y$sample_type <- gsub("SRR55710","SRR",Y$sample_type)
# Y$sample_type <- gsub("ERR25862","ERR_A",Y$sample_type)
Y$sample_type <- gsub("SRR609","SRR_B",Y$sample_type)
Y$sample_type <- gsub("[0-9]","",Y$sample_type)

#color
library(RColorBrewer)
RColorBrewer::display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")
color<-c("#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5")
# Legend Order
Y$sample_type <- factor(Y$sample_type, levels = c("AT","ARP","ODP",
                                                  "nineyear","SRR_B"))

p <- ggplot(Y) +
  geom_segment(aes(x = X1, y = X2,
                   xend = Dim1, yend = Dim2, color=sample_type),
               # geom_segment 绘制两点间的直线
               size = 0.75,linetype="dashed",alpha=0.7) +
  geom_point(aes(X1, X2, color =sample_type),shape=16,size = 3,alpha=0.5) +
  geom_point(aes(Dim1,Dim2,color = sample_type),shape=17,size = 3,alpha=0.5)+
  scale_color_manual("Sample type", values = color, 
                     labels = c(c(expression(Aeration~tank), 
                                  expression(Aeration~tank~PM[2.5]), 
                                  expression(Outdoor~PM[2.5]),
                                  "Nineyear",
                                  "SRR_B")
                                )) +
  guides(color = guide_legend(label.hjust = 0)) + 
  theme_bw() + labs(title="Procrustes analysis") +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  # annotate('text', label = sprintf('M^2 == 0.0348 '),
  #          x = 0.15, y = 0.29, size = 5, parse = TRUE) +
  # annotate('text', label = 'P<0.01',
  #          x = 0.15, y = 0.26, size = 5, parse = TRUE) +
  annotate('text', label = sprintf(''),
           x = 0.15, y = 0.29, size = 5, parse = TRUE) +
  annotate('text', label = '',
           x = 0.15, y = 0.26, size = 5, parse = TRUE) +
  theme(title = element_text(size=14),
        axis.title = element_text(size=15),
        legend.title= element_text(size=13),
        legend.text = element_text(size=13),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) #transparent legend bg)
print(p)

# # Save
# ggsave("ARG_species_procrustes.png", p, path = "../../airborne_arg_uwtp_result/Figure/procrustes",
#        width = 7, height = 5, units = "in", bg='transparent') # save to png format

# Create legend manually
Sample <- c("AT", "ARP", "ODP")
legend_two <- c("Antibiotic resistome","Kraken2 species","Kraken2 species")
value <- c(1,2,3)
df <- data.frame(Sample = Sample, legend_two = legend_two, value = value)
p <- ggplot(df, aes(x=Sample, y=value, fill=Sample)) + 
  geom_bar(stat="identity",alpha=0.5)+
  scale_fill_manual("Sample", values = color)+
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))
print(p)
p<-ggplot(df) + 
  geom_point(aes(value, value, shape=legend_two,size = 3,alpha=0.5))+
  scale_fill_manual("Sample", values = color)+
  theme_bw() +
  theme(legend.text = element_text(size=8),
        #legend.key.size = unit(0.5, 'cm'),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'))
print(p)
# # Save
# ggsave("procrustes_legend.png", p, path = "../../airborne_arg_uwtp_result/Figure/procrustes",
#        width = 7, height = 5, units = "in", bg='transparent') # save to png format
