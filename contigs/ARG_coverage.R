# ARG_coverage.R

library(tidyverse)
library(openxlsx)
# SARG
## Import coverage file & bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_bowtie2/SARG/coverage", pattern = "*_SARG.sam.map.txt")
z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_bowtie2/SARG/coverage/", temp[1], sep = ''))
z$SampleID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_bowtie2/SARG/coverage/", temp[i], sep = ''))
  z2$SampleID <- temp[i]
  z <- rbind(z, z2)
}
colnames(z) <- c("ID","Avg_fold","Length","Ref_GC","Covered_percent","Covered_based",
                 "Plus_reads","Minus_reads","Read_GC","Median_fold","Std_Dev","SampleID") 
z$SampleID <- gsub("_SARG.sam.map.txt", "", z$SampleID)
## Import base number file
nbase <- read.xlsx("../../airborne_arg_uwtp_result/read_base_count.xlsx", sheet = 2)
## Join coverage file & base file
SARG_coverage <- full_join(z, nbase, by = "SampleID")
## Calculate coverage
SARG_coverage$coverage <- SARG_coverage$Avg_fold/(SARG_coverage$base/1000000000)
# Import Diamond file
source("./contigs_ARG_diamond.R")
# Join with ARG and taxonomy
SARG_coverage$ORF_SampleID <- paste(SARG_coverage$ID,
                                    SARG_coverage$SampleID,
                                    sep = "_")
SARG$ORF_SampleID <- paste(SARG$ORF,
                           SARG$SampleID,
                           sep = "_")
SARG_coverage <- left_join(SARG_coverage, SARG, by = "ORF_SampleID")
SARG_coverage <- SARG_coverage %>% select(ORF,Avg_fold,Length,Std_Dev,
                                          SampleID.x,coverage,subtype,type,
                                          contigs,`taxonomy ID`)


# Import mpa format report
temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_kraken2/", pattern = "*_contigs_kraken2.report")
z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_kraken2/", temp[1], sep = ''),
              sep = "\t",quote = "")
z$SampleID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_kraken2/", temp[i], sep = ''),
                  sep = "\t",quote = "")
  z2$SampleID <- temp[i]
  z <- rbind(z, z2)
}
z <- z %>% select(V1)
z <- unique(z)

z <- z %>%
  separate(V1, sep="\\|",into=c("Domain","Phylum","Class","Order","Family","Genus","Species"))
z_na <- z[!complete.cases(z), ]
z_na <- z_na %>% filter(!grepl("g__", Genus))
z_na <- z_na %>% filter(!(grepl("f__", Family)&
                           (grepl(NA, Genus))))
z_na <- z_na %>% filter(!(grepl("o__", Order)&
                           (grepl(NA, Family))&
                           (grepl(NA, Genus))
                           ))
z_na <- z_na %>% filter(!(grepl("c__", Class)&
                           (grepl(NA, Order))&
                           (grepl(NA, Family))&
                           (grepl(NA, Genus))
                           ))
z_na <- z_na %>% filter(!(grepl("p__", Phylum)&
                           (grepl(NA, Class))&
                           (grepl(NA, Order))&
                           (grepl(NA, Family))&
                           (grepl(NA, Genus))
                           ))
z_na <- z_na[-1,]

z_na2 <- z_na[grep("c__",z_na$Phylum),]
z_na2 <- z_na2[,-7]
colnames(z_na2) <- c("Domain","Class","Order","Family","Genus","Species")


z_na3 <- z_na[grep("f__",z_na$Phylum),]
z_na3 <- z_na3[,-(5:7)]
colnames(z_na3) <- c("Domain","Family","Genus","Species")

z_na4 <- z_na[grep("s__",z_na$Phylum),]
z_na4 <- z_na4[,-(3:7)]
colnames(z_na4) <- c("Domain","Species")
## bind dataframe
z <- z %>% anti_join(z_na)
z <- rbind.fill(z,z_na2,z_na3,z_na4)
## Remove taxonomy header
z$Domain <- gsub("d__", "", z$Domain)
z$Phylum <- gsub("p__", "", z$Phylum)
z$Class <- gsub("c__", "", z$Class)
z$Order <- gsub("o__", "", z$Order)
z$Family <- gsub("f__", "", z$Family)
z$Genus <- gsub("g__", "", z$Genus)
z$Species <- gsub("s__", "", z$Species)
## Split output format column to fit mpa report format
SARG_coverage <- separate(SARG_coverage, col='taxonomy ID', into=c('taxonomy','taxID'), sep=' \\(taxid')
SARG_coverage$taxID <- gsub(")", "", SARG_coverage$taxID)
SARG_coverage$taxID <- gsub(" ", "", SARG_coverage$taxID)
## Join mpa format with SARG_coverage
### Species
sp <- SARG_coverage %>% filter(taxonomy %in% z$Species)
colnames(sp)[which(names(sp) == 'taxonomy')] <- 'Species'
sp <- sp %>% left_join(z, by='Species')
### Genus
gen <- SARG_coverage %>% filter(taxonomy %in% z$Genus)
colnames(gen)[which(names(gen) == 'taxonomy')] <- 'Genus'
z <- z %>% select(!(Species))
z <- unique(z)
gen <- gen %>% left_join(z, by='Genus')
### Family
fam <- SARG_coverage %>% filter(taxonomy %in% z$Family)
colnames(fam)[which(names(fam) == 'taxonomy')] <- 'Family'
z <- z %>% select(!(Genus))
z <- unique(z)
fam <- left_join(fam, z, by='Family')
### Order
ord <- SARG_coverage %>% filter(taxonomy %in% z$Order)
colnames(ord)[which(names(ord) == 'taxonomy')] <- 'Order'
z <- z %>% select(!(Family))
z <- unique(z)
ord <- left_join(ord, z, by='Order')
### Class
cl <- SARG_coverage %>% filter(taxonomy %in% z$Class)
colnames(cl)[which(names(cl) == 'taxonomy')] <- 'Class'
z <- z %>% select(!(Order))
z <- unique(z)
cl <- left_join(cl, z, by='Class')
### Phylum
phy <- SARG_coverage %>% filter(taxonomy %in% z$Phylum)
colnames(phy)[which(names(phy) == 'taxonomy')] <- 'Phylum'
z <- z %>% select(!(Class))
z <- unique(z)
phy <- left_join(phy, z, by='Phylum')
### Domain
dom <- SARG_coverage %>% filter(taxonomy %in% z$Domain)
colnames(dom)[which(names(dom) == 'taxonomy')] <- 'Domain'
z <- z %>% select(!(Phylum))
z <- unique(z)
dom <- left_join(dom, z, by='Domain')
# Bind each taxonomy level
tmp <- rbind.fill(dom,phy,cl,ord,fam,gen,sp)
tmp2 <- SARG_coverage %>% anti_join(tmp)
unique(tmp2$taxonomy)
