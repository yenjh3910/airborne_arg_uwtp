# Assembly_ARG_Risk.R

library(tidyverse)
library(openxlsx)
# SARG
## Import file & bind
temp <- list.files(path = "../../airborne_arg_uwtp_result/contigs_diamond/SARG/", pattern = "*_contigs.SARG.dmnd")
z<-read.table(paste("../../airborne_arg_uwtp_result/contigs_diamond/SARG/", temp[1], sep = ''))
z$SampleID <- temp[1]
for (i in 2:length(temp)){
  z2 = read.table(paste("../../airborne_arg_uwtp_result/contigs_diamond/SARG/", temp[i], sep = ''))
  z2$SampleID <- temp[i]
  z <- rbind(z, z2)
}
colnames(z) <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                 "qstart","qend","sstart","send","evalue","bitscore","SampleID") 
z$SampleID <- gsub("_contigs.SARG.dmnd", "", z$SampleID)
## Join with SARG Risk database
rank_db <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/ARG_rank_db.xlsx",
                     sheet = 1)
rank_db <- rank_db %>% select(variant,rank,Type,Subtype)
colnames(rank_db)[1] <- "sseqid"
z <- left_join(z, rank_db, by="sseqid")
## Add contigs column
z$ORF <- z$qseqid
z <- z %>% separate(qseqid, c("tmp1","tmp2","tmp3"), sep = "_")
z$contigs <- paste(z$tmp1,z$tmp2,sep = "_")
z <- z %>% select(!(tmp1)) %>% select(!(tmp2)) %>% select(!(tmp3))
## Filter by identity & evalue
SARG <- z %>% filter(pident >= 70) %>% filter(evalue <= 1e-10) %>% filter(evalue <= 1e-10)
## Import contigs kraken file
source("./contigs_kraken2.R")
contigs_kraken <- contigs_kraken %>% select(contigs, `taxonomy ID`,SampleID)
contigs_kraken$contigs_SampleID <- paste(contigs_kraken$contigs,
                                         contigs_kraken$SampleID,
                                         sep = "_")
contigs_kraken <- contigs_kraken %>% select(contigs_SampleID, `taxonomy ID`)
# Join diamond & kraken
SARG$contigs_SampleID <- paste(SARG$contigs,
                               SARG$SampleID,
                               sep = "_")
SARG <- left_join(SARG, contigs_kraken, by = "contigs_SampleID")
SARG <- SARG %>% select(SampleID,rank,Subtype,Type,ORF,contigs,`taxonomy ID`)

##### Until here, we finish diamond, but coverage not yet #####


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
### OP1 sample do not match any ARG, so remove ODP1 row ###
SARG_coverage <- SARG_coverage %>% filter(!(SampleID == "ODP1"))
## Calculate coverage
SARG_coverage$coverage <- SARG_coverage$Avg_fold/(SARG_coverage$base/1000000000)

# Join with ARG and taxonomy
SARG_coverage$ORF_SampleID <- paste(SARG_coverage$ID,
                                    SARG_coverage$SampleID,
                                    sep = "_")
SARG$ORF_SampleID <- paste(SARG$ORF,
                           SARG$SampleID,
                           sep = "_")
SARG_coverage <- left_join(SARG_coverage, SARG, by = "ORF_SampleID")
SARG_coverage <- SARG_coverage %>% select(ORF,Avg_fold,Length,Std_Dev,
                                          SampleID.x,coverage,rank,Subtype,Type,
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
# Find out wrong order row
z_na <- z[!complete.cases(z), ]
### Filter by wrong order in phylum
phy_na <- z_na %>% filter(!grepl("p__", Phylum))
phy_na <- phy_na %>% filter(!(is.na(Phylum)))
### Filter by wrong order in class
cla_na <- z_na %>% filter(!grepl("c__", Class))
cla_na <- cla_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (is.na(Class))
))
cla_na <- cla_na %>% filter(!(is.na(Phylum)))
### Filter by wrong order in order
ord_na <- z_na %>% filter(!grepl("o__", Order))
ord_na <- ord_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (grepl("c__", Class))&
                                (is.na(Order))
))
ord_na <- ord_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (is.na(Class))
))
ord_na <- ord_na %>% filter(!(is.na(Phylum)))
### Filter by wrong order in family
fam_na <- z_na %>% filter(!grepl("f__", Family))
fam_na <- fam_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (grepl("c__", Class))&
                                (grepl("o__", Order))&
                                (is.na(Family))
))
fam_na <- fam_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (grepl("c__", Class))&
                                (is.na(Order))
))
fam_na <- fam_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (is.na(Class))
))
fam_na <- fam_na %>% filter(!(is.na(Phylum)))
### Filter by wrong order in genus
gen_na <- z_na %>% filter(!grepl("g__", Genus))
gen_na <- gen_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (grepl("c__", Class))&
                                (grepl("o__", Order))&
                                (grepl("f__", Family))&
                                (is.na(Genus))
))
gen_na <- gen_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (grepl("c__", Class))&
                                (grepl("o__", Order))&
                                (is.na(Family))
))
gen_na <- gen_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (grepl("c__", Class))&
                                (is.na(Order))
))
gen_na <- gen_na %>% filter(!(grepl("d__", Domain)&
                                (grepl("p__", Phylum))&
                                (is.na(Class))
))
gen_na <- gen_na %>% filter(!(grepl("d__", Domain)&
                                (is.na(Phylum))
))
# Bind wrong order df
z_na <- bind_rows(phy_na,cla_na,ord_na,fam_na,gen_na) %>% unique()
# Correct order df
z <- z %>% anti_join(z_na)
# Edit manually in excel
#### write.csv(z_na, "./contigs_taxa_WrongOrder.csv", quote = FALSE, row.names = FALSE) ####
correct_z_na <- read_csv("./contigs_taxa_WrongOrder.csv")
## bind dataframe
z <- rbind(z,correct_z_na) %>% unique()
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
ARG_bind_coverage <- bind_rows(dom,phy,cl,ord,fam,gen,sp)


manual_mpa <- read.xlsx("../../airborne_arg_uwtp_result/contigs_bowtie2/SARG/ARG_TaxonomyAnnotateManually.xlsx",
                        sheet =1)
manual_mpa <- manual_mpa %>%
  separate(mpa, sep=";",into=c("Domain","Phylum","Class","Order","Family","Genus","Species"))
manual_mpa[manual_mpa == ''] <- NA
manual_mpa <- manual_mpa %>% select(!(taxonomy))
# Finish binding manual mpa with origin df
ARG_bind_coverage$taxID <- as.numeric(ARG_bind_coverage$taxID)
final_ARG_coverage <- bind_rows(ARG_bind_coverage,manual_mpa)
final_ARG_coverage <- final_ARG_coverage %>% select(!(taxID))
## Finetune format
## Remove blank space
final_ARG_coverage$Domain <- trimws(final_ARG_coverage$Domain)
final_ARG_coverage$Phylum <- trimws(final_ARG_coverage$Phylum)
final_ARG_coverage$Class <- trimws(final_ARG_coverage$Class)
final_ARG_coverage$Order <- trimws(final_ARG_coverage$Order)
final_ARG_coverage$Family <- trimws(final_ARG_coverage$Family)
final_ARG_coverage$Genus <- trimws(final_ARG_coverage$Genus)
final_ARG_coverage$Species <- trimws(final_ARG_coverage$Species)
## Remove \u00a0 type space
utf8::utf8_print(unique(final_ARG_coverage$Phylum), utf8 = FALSE)
final_ARG_coverage$Phylum <- gsub("^\u00a0", "", final_ARG_coverage$Phylum)
final_ARG_coverage$Class <- gsub("^\u00a0", "", final_ARG_coverage$Class)
final_ARG_coverage$Order <- gsub("^\u00a0", "", final_ARG_coverage$Order)
final_ARG_coverage$Family <- gsub("^\u00a0", "", final_ARG_coverage$Family)
final_ARG_coverage$Genus <- gsub("^\u00a0", "", final_ARG_coverage$Genus)
final_ARG_coverage$Species <- gsub("^\u00a0", "", final_ARG_coverage$Species)
## Check
unique(final_ARG_coverage$Domain)
unique(final_ARG_coverage$Phylum)
unique(final_ARG_coverage$Class)
unique(final_ARG_coverage$Order)
unique(final_ARG_coverage$Family)
unique(final_ARG_coverage$Genus)
unique(final_ARG_coverage$Species)
# NA to Unclassfied
final_ARG_coverage$Domain <- replace(final_ARG_coverage$Domain, final_ARG_coverage$Domain == "unclassified", "Unclassified")
final_ARG_coverage$Phylum <- replace_na(final_ARG_coverage$Phylum, "Unclassified")
final_ARG_coverage$Class <- replace_na(final_ARG_coverage$Class, "Unclassified")
final_ARG_coverage$Order <- replace_na(final_ARG_coverage$Order, "Unclassified")
final_ARG_coverage$Family <- replace_na(final_ARG_coverage$Family, "Unclassified")
final_ARG_coverage$Genus <- replace_na(final_ARG_coverage$Genus, "Unclassified")
final_ARG_coverage$Species <- replace_na(final_ARG_coverage$Species, "Unclassified")
# Simplify table
final_ARG_coverage <- final_ARG_coverage %>% select(!(Avg_fold)) %>% select(!(Length)) %>% select(!(Std_Dev))
colnames(final_ARG_coverage)[2] <- "Sample" 
final_ARG_coverage$Sample_type <- gsub("1|2|3|4|5","",final_ARG_coverage$Sample)
final_ARG_coverage <- final_ARG_coverage %>% select(!(subtype)) %>% select(!(type))

## Fix bug
final_ARG_coverage <- final_ARG_coverage %>% select(!(rank)) %>% select(!(Subtype)) %>% select((!Type)) 
SARG_coverage <- SARG_coverage %>% select(ORF,rank,Subtype,Type)
final_ARG_coverage <- left_join(final_ARG_coverage,SARG_coverage,by='ORF')


###### Data analysis ########
####### Phylum #########
# Sum ARG coverage in same sample type
ARG_Phylum <- final_ARG_coverage  %>% ungroup() %>% 
  group_by(Sample_type, rank, Phylum) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
## Filter notassessed rank
ARG_Phylum <- ARG_Phylum %>% filter(!(rank=='notassessed'))
## Order phylum
ARG_Phylum$Phylum <- factor(ARG_Phylum$Phylum, 
                            levels = c("Pseudomonadota","Actinomycetota","Bacteroidota","Bacillota","Thermodesulfobacteriota",
                                       "Chlorobiota","Fusobacteriota","Campylobacterota","Unclassified"))
#### Order sample_type
ARG_Phylum$Sample_type <- factor(ARG_Phylum$Sample_type, 
                                 levels = c("AT","ARP","ODP"))
## Facet with sample type
ggplot(ARG_Phylum, aes(x = rank, y = coverage, fill = Phylum)) + 
  geom_bar(stat="identity",position = 'fill') + 
  facet_wrap(~Sample_type)+
  theme_bw()+ 
  xlab("") + ylab("Proportion of ARG Coverage") + 
  scale_fill_brewer(palette="Set3")
## Facet with rank
ggplot(ARG_Phylum, aes(x = Sample_type, y = coverage, fill = Phylum)) + 
  geom_bar(stat="identity",position = 'fill') + 
  facet_wrap(~rank)+
  theme_bw()+ 
  xlab("") + ylab("Proportion of ARG Coverage") + 
  scale_fill_brewer(palette="Set3")+
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5])))+
  theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent')) #transparent legend bg)

# Filter ODP
aer_ARG_Phylum <- ARG_Phylum %>% filter(!(Sample_type=='ODP'))
p<-ggplot(aer_ARG_Phylum, aes(x = Sample_type, y = coverage, fill = Phylum)) + 
  geom_bar(stat="identity",position = 'fill') + 
  facet_wrap(~rank,nrow = 1)+
  theme_bw()+ 
  xlab("") + ylab("Proportion of ARG Coverage") + 
  scale_fill_brewer(palette="Set3")+
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5])))+
  theme(axis.text.x = element_text(size = 8, angle = 70, hjust=1),
        #title.text.y = element_text(size = 10),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent')) #transparent legend bg)
print(p)
# ggsave("risk_rank_phylum.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG_coverage",
#        width = 6, height = 3.7,
#        units = "in", bg='transparent') # save to png format

aer_ARG_Phylum <- aer_ARG_Phylum %>% group_by(Sample_type,rank) %>% mutate(group_sum=sum(coverage))
aer_ARG_Phylum<- aer_ARG_Phylum %>% mutate(proportion=(coverage/group_sum))
####### Class #########
# Sum ARG coverage in same sample type
ARG_Class <- final_ARG_coverage  %>% ungroup() %>% 
  group_by(Sample_type, rank, Class) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
## Filter notassessed rank
ARG_Class <- ARG_Class %>% filter(!(rank=='notassessed'))
#### Order sample_type
ARG_Class$Sample_type <- factor(ARG_Class$Sample_type, 
                                 levels = c("AT","ARP","ODP"))
## Facet with sample type
ggplot(ARG_Class, aes(x = rank, y = coverage, fill = Class)) + 
  geom_bar(stat="identity",position = 'fill') + 
  facet_wrap(~Sample_type)+
  theme_bw()+ 
  xlab("") + ylab("Proportion of ARG Coverage")
## Facet with rank
ggplot(ARG_Class, aes(x = Sample_type, y = coverage, fill = Class)) + 
  geom_bar(stat="identity",position = 'fill') + 
  facet_wrap(~rank)+
  theme_bw()+ 
  xlab("") + ylab("Proportion of ARG Coverage")

####### Genus #########
# Sum ARG coverage in same sample type
ARG_Genus <- final_ARG_coverage  %>% ungroup() %>% 
  group_by(Sample_type, rank, Genus) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
## Filter notassessed rank
ARG_Genus <- ARG_Genus %>% filter(!(rank=='notassessed'))
#### Order sample_type
ARG_Genus$Sample_type <- factor(ARG_Genus$Sample_type, 
                                levels = c("AT","ARP","ODP"))
## Facet with sample type
ggplot(ARG_Genus, aes(x = rank, y = coverage, fill = Genus)) + 
  geom_bar(stat="identity",position = 'fill') + 
  facet_wrap(~Sample_type)+
  theme_bw()+ 
  xlab("") + ylab("Proportion of ARG Coverage")
## Facet with rank
ggplot(ARG_Genus, aes(x = Sample_type, y = coverage, fill = Genus)) + 
  geom_bar(stat="identity",position = 'fill') + 
  facet_wrap(~rank)+
  theme_bw()+ 
  xlab("") + ylab("Proportion of ARG Coverage")


# View
rank_I_AT_Pseudomonadota <- final_ARG_coverage %>% filter(rank=='I') %>% 
                                    filter(Sample_type == 'AT') %>% 
                                    filter(Phylum == 'Pseudomonadota')
rank_I_ARP_Pseudomonadota <- final_ARG_coverage %>% filter(rank=='I') %>% 
                                                   filter(Sample_type == 'ARP') %>% 
                                                   filter(Phylum == 'Pseudomonadota')
rank_I_AT <- final_ARG_coverage %>% filter(rank=='I') %>% 
                                    filter(Sample_type == 'AT') %>% 
                                    filter(Phylum == 'Bacillota')
rank_I_ARP <- final_ARG_coverage %>% filter(rank=='I') %>% 
                                     filter(Sample_type == 'ARP') %>% 
                                     filter(Phylum == 'Bacillota')
