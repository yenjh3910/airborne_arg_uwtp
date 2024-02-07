# ARB_compare.R

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

############### Phylum ##############
# Extract phylum
ARG_Phylum <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Phylum) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Phylum$Sample_type <- gsub("1|2|3|4|5","",ARG_Phylum$Sample)
ARG_Phylum <- ARG_Phylum %>% ungroup() %>% group_by(Phylum,Sample_type)%>% mutate(total_coverage = sum(coverage)) %>% 
                    select(!(Sample)) %>% select(!(coverage)) %>% filter(!(Sample_type == 'ODP')) %>% unique()
# Spread dataframe
spread_ARG_Phylum <- spread(ARG_Phylum, key = "Sample_type", value = "total_coverage")
# Filter NA row
spread_ARG_Phylum <- spread_ARG_Phylum %>% filter(!(if_any(everything(), is.na)))
# Add ratio column
spread_ARG_Phylum$ratio <- spread_ARG_Phylum$ARP/(spread_ARG_Phylum$AT+spread_ARG_Phylum$ARP)
# Arrange order
spread_ARG_Phylum <- spread_ARG_Phylum %>% 
  arrange(desc(ratio))
spread_ARG_Phylum$Phylum <- factor(spread_ARG_Phylum$Phylum , levels=spread_ARG_Phylum$Phylum)
# Plot
ggplot(spread_ARG_Phylum, aes(x=Phylum, y=ratio, fill=Phylum)) +
  geom_bar(stat="identity") +
  xlab("Phylum") + ylab('Enriched taxa proportion (ARP/(ARP+AT))') + 
  theme_bw() +
  scale_fill_manual(values=c("#FB8072", "#80B1D3", "#FDB462",
                             "#B3DE69","#BEBADA"))+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1))+
  theme(axis.title = element_text(size = 18),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 12),
        #legend.key.size = unit(1.1, 'cm'),
        legend.title = element_text(size = 17),
        legend.text = element_text(size = 17),
        legend.position = 'none',
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'))



############### Class ##############
# Extract Class
ARG_Class <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Class) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Class$Sample_type <- gsub("1|2|3|4|5","",ARG_Class$Sample)
ARG_Class <- ARG_Class %>% ungroup() %>% group_by(Class,Sample_type)%>% mutate(total_coverage = sum(coverage)) %>% 
  select(!(Sample)) %>% select(!(coverage)) %>% filter(!(Sample_type == 'ODP')) %>% unique()
# Spread dataframe
spread_ARG_Class <- spread(ARG_Class, key = "Sample_type", value = "total_coverage")
# Filter NA row
spread_ARG_Class <- spread_ARG_Class %>% filter(!(if_any(everything(), is.na)))
# Add ratio column
spread_ARG_Class$ratio <- spread_ARG_Class$ARP/(spread_ARG_Class$AT+spread_ARG_Class$ARP)
# Arrange order
spread_ARG_Class <- spread_ARG_Class %>% 
  arrange(desc(ratio))
spread_ARG_Class$Class <- factor(spread_ARG_Class$Class , levels=spread_ARG_Class$Class)
# Plot
ggplot(spread_ARG_Class, aes(x=Class, y=ratio, fill=Class)) +
  geom_bar(stat="identity") +
  xlab("Class") + ylab('Enriched taxa proportion (ARP/(ARP+AT))') + 
  theme_bw()+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        legend.position = 'none')

############### Order ##############
# Extract Order
ARG_Order <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Order) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Order$Sample_type <- gsub("1|2|3|4|5","",ARG_Order$Sample)
ARG_Order <- ARG_Order %>% ungroup() %>% group_by(Order,Sample_type)%>% mutate(total_coverage = sum(coverage)) %>% 
  select(!(Sample)) %>% select(!(coverage)) %>% filter(!(Sample_type == 'ODP')) %>% unique()
# Spread dataframe
spread_ARG_Order <- spread(ARG_Order, key = "Sample_type", value = "total_coverage")
# Filter NA row
spread_ARG_Order <- spread_ARG_Order %>% filter(!(if_any(everything(), is.na)))
# Add ratio column
spread_ARG_Order$ratio <- spread_ARG_Order$ARP/(spread_ARG_Order$AT+spread_ARG_Order$ARP)
# Arrange order
spread_ARG_Order <- spread_ARG_Order %>% 
  arrange(desc(ratio))
spread_ARG_Order$Order <- factor(spread_ARG_Order$Order , levels=spread_ARG_Order$Order)
# Plot
ggplot(spread_ARG_Order, aes(x=Order, y=ratio, fill=Order)) +
  geom_bar(stat="identity") +
  xlab("Order") + ylab('Enriched taxa proportion (ARP/(ARP+AT))') + 
  theme_bw()+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        legend.position = 'none')


############### Family ##############
# Extract Family
ARG_Family <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Family) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Family$Sample_type <- gsub("1|2|3|4|5","",ARG_Family$Sample)
ARG_Family <- ARG_Family %>% ungroup() %>% group_by(Family,Sample_type)%>% mutate(total_coverage = sum(coverage)) %>% 
  select(!(Sample)) %>% select(!(coverage)) %>% filter(!(Sample_type == 'ODP')) %>% unique()
# Spread dataframe
spread_ARG_Family <- spread(ARG_Family, key = "Sample_type", value = "total_coverage")
# Filter NA row
spread_ARG_Family <- spread_ARG_Family %>% filter(!(if_any(everything(), is.na)))
# Add ratio column
spread_ARG_Family$ratio <- spread_ARG_Family$ARP/(spread_ARG_Family$AT+spread_ARG_Family$ARP)
# Arrange Family
spread_ARG_Family <- spread_ARG_Family %>% 
  arrange(desc(ratio))
spread_ARG_Family$Family <- factor(spread_ARG_Family$Family , levels=spread_ARG_Family$Family)
# Plot
ggplot(spread_ARG_Family, aes(x=Family, y=ratio, fill=Family)) +
  geom_bar(stat="identity") +
  xlab("Family") + ylab('Enriched taxa proportion (ARP/(ARP+AT))') + 
  theme_bw()+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        legend.position = 'none')

############### Genus ##############
# Extract Genus
ARG_Genus <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Genus) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Genus$Sample_type <- gsub("1|2|3|4|5","",ARG_Genus$Sample)
ARG_Genus <- ARG_Genus %>% ungroup() %>% group_by(Genus,Sample_type)%>% mutate(total_coverage = sum(coverage)) %>% 
  select(!(Sample)) %>% select(!(coverage)) %>% filter(!(Sample_type == 'ODP')) %>% unique()
# Spread dataframe
spread_ARG_Genus <- spread(ARG_Genus, key = "Sample_type", value = "total_coverage")
# Filter NA row
spread_ARG_Genus <- spread_ARG_Genus %>% filter(!(if_any(everything(), is.na)))
# Add ratio column
spread_ARG_Genus$ratio <- spread_ARG_Genus$ARP/(spread_ARG_Genus$AT+spread_ARG_Genus$ARP)
# Arrange Genus
spread_ARG_Genus <- spread_ARG_Genus %>% 
  arrange(desc(ratio))
spread_ARG_Genus$Genus <- factor(spread_ARG_Genus$Genus , levels=spread_ARG_Genus$Genus)
# Plot
ggplot(spread_ARG_Genus, aes(x=Genus, y=ratio, fill=Genus)) +
  geom_bar(stat="identity") +
  xlab("Genus") + ylab('Enriched taxa proportion (ARP/(ARP+AT))') + 
  theme_bw()+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        legend.position = 'none')


############### Species ##############
# Extract Species
ARG_Species <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Species) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Species$Sample_type <- gsub("1|2|3|4|5","",ARG_Species$Sample)
ARG_Species <- ARG_Species %>% ungroup() %>% group_by(Species,Sample_type)%>% mutate(total_coverage = sum(coverage)) %>% 
  select(!(Sample)) %>% select(!(coverage)) %>% filter(!(Sample_type == 'ODP')) %>% unique()
# Spread dataframe
spread_ARG_Species <- spread(ARG_Species, key = "Sample_type", value = "total_coverage")
# Filter NA row
spread_ARG_Species <- spread_ARG_Species %>% filter(!(if_any(everything(), is.na)))
# Add ratio column
spread_ARG_Species$ratio <- spread_ARG_Species$ARP/(spread_ARG_Species$AT+spread_ARG_Species$ARP)
# Arrange Species
spread_ARG_Species <- spread_ARG_Species %>% 
  arrange(desc(ratio))
spread_ARG_Species$Species <- factor(spread_ARG_Species$Species , levels=spread_ARG_Species$Species)
# Plot
ggplot(spread_ARG_Species, aes(x=Species, y=ratio, fill=Species)) +
  geom_bar(stat="identity") +
  xlab("Species") + ylab('Enriched taxa proportion (ARP/(ARP+AT))') + 
  theme_bw()+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        legend.position = 'none')




# Fill with order color
reference <- final_ARG_coverage %>% select(Order,Genus) %>% 
             unique() %>% filter(!(Genus=='Unclassified'))
spread_ARG_Genus_fill <- left_join(spread_ARG_Genus,reference,)
spread_ARG_Genus_fill <- spread_ARG_Genus_fill %>% filter(!(Genus=='Unclassified'))
# Arrange Genus
spread_ARG_Genus_fill <- spread_ARG_Genus_fill %>% 
  arrange(desc(ratio))
spread_ARG_Genus_fill$Genus <- factor(spread_ARG_Genus_fill$Genus , levels=spread_ARG_Genus_fill$Genus)
# Plot
ggplot(spread_ARG_Genus_fill, aes(x=Genus, y=ratio, fill=Order)) +
  geom_bar(stat="identity") +
  xlab("Genus") + ylab('Enriched taxa proportion (ARP/(ARP+AT))') + 
  theme_bw()+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        legend.position = 'right')
# Replace to other
spread_ARG_Genus_fill$Order[spread_ARG_Genus_fill$Order=='Aeromonadales']<-'Others'
spread_ARG_Genus_fill$Order[spread_ARG_Genus_fill$Order=='Desulfovibrionales']<-'Others'
spread_ARG_Genus_fill$Order[spread_ARG_Genus_fill$Order=='Hyphomicrobiales']<-'Others'
spread_ARG_Genus_fill$Order[spread_ARG_Genus_fill$Order=='Kitasatosporales']<-'Others'
spread_ARG_Genus_fill$Order[spread_ARG_Genus_fill$Order=='Propionibacteriales']<-'Others'
spread_ARG_Genus_fill$Order[spread_ARG_Genus_fill$Order=='Rhodobacterales']<-'Others'
spread_ARG_Genus_fill$Order[spread_ARG_Genus_fill$Order=='Thiotrichales']<-'Others'
# Rearrange
spread_ARG_Genus_fill$Order <- factor(spread_ARG_Genus_fill$Order, levels = c("Mycobacteriales" ,    "Burkholderiales"  ,   "Xanthomonadales"   , 
                                                      "Flavobacteriales",    "Pseudomonadales"  ,   "Moraxellales"      ,  "Bacillales"         ,
                                                      "Rhodocyclales"   ,    "Aeromonadales"    ,   "Enterobacterales"  ,  "Lactobacillales"    ,
                                                      "Unclassified","Others"))
# Plot again
ggplot(spread_ARG_Genus_fill, aes(x=Genus, y=ratio, fill=Order)) +
  geom_bar(stat="identity") +
  xlab("Genus") + ylab('Enriched taxa proportion (ARP/(ARP+AT))') + 
  theme_bw()+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        legend.position = 'right') +
  scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                             "#FDB462","#FCCDE5","#BC80BD","#CCEBC5","#D9D9D9","#D9D9D9"))
spread_ARG_Genus_fill <- spread_ARG_Genus_fill %>% filter(ratio>0.535)

# Rearrange again
spread_ARG_Genus_fill$Order <- factor(spread_ARG_Genus_fill$Order, levels = c("Mycobacteriales" ,    "Burkholderiales"  ,   "Xanthomonadales"   , 
                                                                          "Flavobacteriales",    "Pseudomonadales"  ,   "Moraxellales"      ,  "Bacillales"         ,
                                                                              "Rhodocyclales"   ,    "Aeromonadales"    ,   "Enterobacterales"  ,  "Lactobacillales"))
# Plot again
ggplot(spread_ARG_Genus_fill, aes(x=Genus, y=ratio, fill=Order)) +
  geom_bar(stat="identity") +
  xlab("Genus") + ylab('Enriched taxa proportion (ARP/(ARP+AT))') + 
  theme_bw()+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        legend.position = 'right') +
  scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072",
                             "#FDB462","#FCCDE5","#CCEBC5"))
# Plot without legend
p <- ggplot(spread_ARG_Genus_fill, aes(x=Genus, y=ratio, fill=Order)) +
  geom_bar(stat="identity") +
  xlab("Genus") + ylab('Enriched ARG-ORF proportion (ARP/(ARP+AT))') + 
  theme_bw()+
  scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072",
                             "#FDB462","#FCCDE5","#CCEBC5"))+
  coord_cartesian(ylim = c(0.5, 1)) +
  theme(axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 20),
        axis.text.x=element_text(size = 20, angle = 70, hjust=1),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 18),
        legend.position = 'none',
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'))
print(p)

# # Save
# ggsave("dominant_ARG-ORF_proportion.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG_coverage",
#        width = 16, height = 9,
#        units = "in", bg='transparent') # save to png format


# # Save as thinner version
# ggsave("dominant_ARG-ORF_proportion.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG_coverage",
#        width = 10, height = 9,
#        units = "in", bg='transparent') # save to png format
