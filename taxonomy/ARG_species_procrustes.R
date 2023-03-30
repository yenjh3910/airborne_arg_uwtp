# taxonomy_barplot.R

library(tidyverse)
library(vegan)
library(openxlsx)

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
arg_subtype <- read.xlsx("D:/ARG_project/Shell/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
row.names(arg_subtype) <- arg_subtype[,1]
arg_subtype <- arg_subtype[,-1]
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
Y$sample_type <- gsub("1|2|3|4|5","",Y$sample_type)
#color
library(RColorBrewer)
RColorBrewer::display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")
color<-c("#FB8072","#80B1D3","#FDB462")


p <- ggplot(Y) +
  geom_segment(aes(x = X1, y = X2,
                   xend = Dim1, yend = Dim2, color=sample_type),
               # geom_segment 绘制两点间的直线
               size = 0.75,linetype="dashed",alpha=0.7) +
  geom_point(aes(X1, X2, color =sample_type),shape=16,size = 3,alpha=0.5) +
  geom_point(aes(Dim1,Dim2,color = sample_type),shape=17,size = 3,alpha=0.5) +
  scale_color_manual("Sample type", values = color, 
                     labels = c(c(expression(Aeration~tank), 
                                  expression(Aeration~tank~PM[2.5]), 
                                  expression(Outdoor~PM[2.5])))) +
  guides(color = guide_legend(label.hjust = 0)) + 
  theme_bw() + labs(title="Procrustes analysis") +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  annotate('text', label = sprintf('M^2 == 0.2799 '),
           x = 0.20, y = 0.29, size = 5, parse = TRUE) +
  annotate('text', label = 'P==0.001',
           x = 0.20, y = 0.26, size = 5, parse = TRUE) +
  theme(title = element_text(size=14),
        axis.title = element_text(size=15),
        legend.title= element_text(size=13),
        legend.text = element_text(size=13),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) #transparent legend bg)
print(p)

# # Save
ggsave("ARG_species_procrustes.png", p, path = "D:/ARG_project/Figure/procrustes",
       width = 7, height = 5, units = "in", bg='transparent') # save to png format