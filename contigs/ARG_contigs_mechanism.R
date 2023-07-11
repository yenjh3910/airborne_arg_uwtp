# ARG_contigs_mechanism.R

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
### OP1 sample do not match any ARG, so remove ODP1 row ###
SARG_coverage <- SARG_coverage %>% filter(!(SampleID == "ODP1"))
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
coverage_without_annotation <- SARG_coverage %>% anti_join(ARG_bind_coverage)
length(unique(coverage_without_annotation$taxonomy)) #### Stop here so far ###
########## Then annotate other tax in excel ############
# write_csv(coverage_without_annotation,
# "../../airborne_arg_uwtp_result/contigs_bowtie2/SARG/ARG_TaxonomyAnnotateManually.csv")
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


# Seperate type-suntype column
final_ARG_coverage <- final_ARG_coverage %>%
  separate(subtype, sep="__",into=c('type2','subtype')) %>% 
  select(!(type2))
# Merge with mechanism data
arg_mechanism_reference <- read.table("../ARG/arg_mechanism_reference.csv", sep = ",",quote = "", header = TRUE)
###arg_mechanism_reference$mechanism[arg_mechanism_reference$subtype == 'Klebsiella pneumoniae OmpK37'] <- 'Antibiotic permeability reduction'
final_ARG_coverage <- left_join(final_ARG_coverage,arg_mechanism_reference)
# Save to result
#write.csv(final_ARG_coverage, "../../airborne_arg_uwtp_result/contigs_bowtie2/SARG/final_ARG_coverage.csv", quote = FALSE, row.names = FALSE)
