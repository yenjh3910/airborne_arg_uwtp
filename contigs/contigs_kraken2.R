# contigs_kraken2.R

library(tidyverse)
# Import file and bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_kraken2/", pattern = "*_contigs_kraken2.output")
z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_kraken2/", temp[1], sep = ''),
              sep = "\t",quote = "")
z$SampleID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_kraken2/", temp[i], sep = ''),
                  sep = "\t",quote = "")
  z2$SampleID <- temp[i]
  z <- rbind(z, z2)
}
colnames(z) <- c("classified","contigs","taxonomy ID","length","LCA mapping","SampleID")
z$SampleID <- gsub("_contigs_kraken2.output", "", z$SampleID)
contigs_kraken <- z

# # Import mpa format report
# temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_kraken2/", pattern = "*_contigs_kraken2.report")
# z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_kraken2/", temp[1], sep = ''),
#               sep = "\t",quote = "")
# z$SampleID <- temp[1]
# for (i in 2:length(temp)){
#   z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_kraken2/", temp[i], sep = ''),
#                   sep = "\t",quote = "")
#   z2$SampleID <- temp[i]
#   z <- rbind(z, z2)
# }
# z <- z %>% select(V1)
# z <- unique(z)
# z <- separate(z, col=V1, into=c('Domain','Phylum'), sep='\\|p__')
# z <- separate(z, col=Domain, into=c('Domain','Class'), sep='\\|c__')
# z <- separate(z, col=Domain, into=c('Domain','Family'), sep='\\|f__')
# z <- separate(z, col=Domain, into=c('Domain','Species'), sep='\\|s__')
# z <- separate(z, col=Phylum, into=c('Phylum','Class'), sep='\\|c__')
# z <- separate(z, col=Phylum, into=c('Phylum','Order'), sep='\\|o__')
# z <- separate(z, col=Phylum, into=c('Phylum','Family'), sep='\\|f__')
# z <- separate(z, col=Phylum, into=c('Phylum','Genus'), sep='\\|g__')
# z <- separate(z, col=Phylum, into=c('Phylum','Species'), sep='\\|s__')
# z <- separate(z, col=Class, into=c('Class','Order'), sep='\\|o__')
# z <- separate(z, col=Class, into=c('Class','Family'), sep='\\|f__')
# z <- separate(z, col=Class, into=c('Class','Genus'), sep='\\|g__')
# z <- separate(z, col=Class, into=c('Class','Species'), sep='\\|s__')
# z <- separate(z, col=Order, into=c('Order','Family'), sep='\\|f__')
# z <- separate(z, col=Order, into=c('Order','Genus'), sep='\\|g__')
# z <- separate(z, col=Order, into=c('Order','Species'), sep='\\|s__')
# z <- separate(z, col=Family, into=c('Family','Genus'), sep='\\|g__')
# z <- separate(z, col=Family, into=c('Family','Species'), sep='\\|s__')
# z <- separate(z, col=Genus, into=c('Genus','Species'), sep='\\|s__')
# z$Domain <- gsub("d__", "", z$Domain)
# 
# # z <- separate(z, col=V1, into=c('Domain','Phylum','Class','Order','Family','Genus','Species'), sep='\\|')
# # z$Domain <- gsub("d__", "", z$Domain)
# # z$Phylum <- gsub("p__", "", z$Phylum)
# # z$Class <- gsub("c__", "", z$Class)
# # z$Order <- gsub("o__", "", z$Order)
# # z$Family <- gsub("f__", "", z$Family)
# # z$Genus <- gsub("g__", "", z$Genus)
# # z$Species <- gsub("s__", "", z$Species)
# ## Split output format column to fit mpa report format
# contigs_kraken <- separate(contigs_kraken, col='taxonomy ID', into=c('taxonomy','ID'), sep=' \\(taxid')
# contigs_kraken$ID <- gsub(")", "", contigs_kraken$ID)
# contigs_kraken$ID <- gsub(" ", "", contigs_kraken$ID)
# ## Join output format and report format
# ### Species
# sp <- contigs_kraken %>% filter(taxonomy %in% z$Species)
# colnames(sp)[which(names(sp) == 'taxonomy')] <- 'Species' 
# sp <- left_join(sp, z, by='Species')
# ### Genus
# gen <- contigs_kraken %>% filter(taxonomy %in% z$Genus)
# colnames(gen)[which(names(gen) == 'taxonomy')] <- 'Genus'
# z <- z %>% select(!(Species))
# z <- unique(z)
# gen <- left_join(gen, z, by='Genus')
# ### Family
# fam <- contigs_kraken %>% filter(taxonomy %in% z$Family)
# colnames(fam)[which(names(fam) == 'taxonomy')] <- 'Family'
# z <- z %>% select(!(Genus))
# z <- unique(z)
# fam <- left_join(fam, z, by='Family')
# ### Order
# ord <- contigs_kraken %>% filter(taxonomy %in% z$Order)
# colnames(ord)[which(names(ord) == 'taxonomy')] <- 'Order'
# z <- z %>% select(!(Family))
# z <- unique(z)
# ord <- left_join(ord, z, by='Order')
# ### Class
# cl <- contigs_kraken %>% filter(taxonomy %in% z$Class)
# colnames(cl)[which(names(cl) == 'taxonomy')] <- 'Class'
# z <- z %>% select(!(Order))
# z <- unique(z)
# cl <- left_join(cl, z, by='Class')
# ### Phylum
# phy <- contigs_kraken %>% filter(taxonomy %in% z$Phylum)
# colnames(phy)[which(names(phy) == 'taxonomy')] <- 'Phylum'
# z <- z %>% select(!(Class))
# z <- unique(z)
# phy <- left_join(phy, z, by='Phylum')
# ### Domain
# dom <- contigs_kraken %>% filter(taxonomy %in% z$Domain)
# colnames(dom)[which(names(dom) == 'taxonomy')] <- 'Domain'
# z <- z %>% select(!(Phylum))
# z <- unique(z)
# dom <- left_join(dom, z, by='Domain')
# 
# tmp <- rbind.fill(dom,phy,cl,ord,fam,gen,sp)
# x <- contigs_kraken %>% filter(classified == "U")
# y <- contigs_kraken %>% filter(classified == "C")
# unique(contigs_kraken$classified)
# notin <- setdiff(y$ID,tmp$ID)
# test2 <- contigs_kraken %>% filter(ID %in% notin)
# unique(test2$ID)
# 
