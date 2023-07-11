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
# Statistic
ARG_Phylum <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Phylum) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Phylum$Sample_type <- gsub("1|2|3|4|5","",ARG_Phylum$Sample)
# Calculate mean and sd (Some microbes do not appear in all sample, so mean()&sd() must wrong)
ARG_Phylum <- ARG_Phylum %>% ungroup() %>% group_by(Phylum, Sample_type) %>% mutate(mean = sum(coverage)/5)
ARG_Phylum <- ARG_Phylum %>% ungroup() %>% group_by(Phylum, Sample_type) %>% 
  mutate(sd = sqrt(sum((coverage-mean)^2/(5-1))))
# Filter necessary column for figure
ARG_Phylum <- ARG_Phylum %>% select(Phylum,Sample_type,mean,sd) %>% unique()
# Add st_mean for figure
ARG_Phylum <- ARG_Phylum %>% group_by(Sample_type) %>% 
  mutate(sd_mean = mean)
ARG_Phylum <- ARG_Phylum %>% arrange(Sample_type)
## sum ARP st_mean
for (i in 1:(nrow(ARG_Phylum %>% filter(Sample_type == "ARP"))-1)) {
  ARG_Phylum[{i},5] <-  ARG_Phylum[{i},5] + 
    sum(ARG_Phylum[{i+1}:nrow(ARG_Phylum %>% filter(Sample_type == "ARP")),5])}
## sum AT st_mean
for (i in (nrow(ARG_Phylum %>% filter(Sample_type == "ARP"))+1):(nrow(ARG_Phylum %>% filter(!(Sample_type == "ODP")))-1)) {
  ARG_Phylum[{i},5] <-  ARG_Phylum[{i},5] + 
    sum(ARG_Phylum[{i+1}:nrow(ARG_Phylum %>% filter(!(Sample_type == "ODP"))),5])}
## sum ODP st_mean
for (i in (nrow(ARG_Phylum %>% filter(!(Sample_type == "ODP")))+1):(nrow(ARG_Phylum)-1)) {
  ARG_Phylum[{i},5] <-  ARG_Phylum[{i},5] + 
    sum(ARG_Phylum[{i+1}:nrow(ARG_Phylum) ,5])}
# Visualize with barplot ()
## Select color
library(RColorBrewer)
RColorBrewer::display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")
phylum_col <- c( "#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#FDB462","#B3DE69","#D9D9D9",
                 "#8DD3C7","#FFFFB3","#BEBADA","#80B1D3","#B3DE69","#FCCDE5","#D9D9D9",
                 "#8DD3C7","#FFFFB3","#B3DE69","#D9D9D9")
ggplot(ARG_Phylum, aes(x = Sample_type, y = mean, fill = Phylum))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(x=Sample_type, ymin=sd_mean, ymax=sd_mean+sd), 
                width=0.4, color=phylum_col, alpha=0.9, size=0.5) + 
  theme_bw()+ 
  scale_fill_brewer(palette="Set3")
### Not plot with errorbar
#### Order phylum
ARG_Phylum <- ARG_Phylum %>% arrange(desc(mean))
unique(ARG_Phylum$Phylum)
ARG_Phylum$Phylum <- factor(ARG_Phylum$Phylum, 
                            levels = c("Pseudomonadota","Actinomycetota","Bacteroidota","Bacillota","Thermodesulfobacteriota",
                                       "Chlorobiota","Fusobacteriota","Campylobacterota","Unclassified"))
#### Order sample_type
ARG_Phylum$Sample_type <- factor(ARG_Phylum$Sample_type, 
                            levels = c("AT","ARP","ODP"))
# Plot
p <- ggplot(ARG_Phylum, aes(x = Sample_type, y = mean, fill = Phylum)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)") + 
  scale_fill_brewer(palette="Set3") +
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
print(p)
# ggsave("ARG_coverage_phylum.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG_coverage",
#        width = 4, height = 5,
#        units = "in", bg='transparent') # save to png format



############### Class ################
ARG_Class <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Class) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Class$Sample_type <- gsub("1|2|3|4|5","",ARG_Class$Sample)
## Calculate mean
ARG_Class <- ARG_Class %>% ungroup() %>% group_by(Class, Sample_type) %>% mutate(mean = sum(coverage)/5)
## Filter necesssary column for figure
ARG_Class <- ARG_Class %>% select(Class,Sample_type,mean) %>% unique()
## Expand color
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
## Order Class
ARG_Class <- ARG_Class %>% arrange(desc(mean))
unique(ARG_Class$Class)
ARG_Class$Class <- factor(ARG_Class$Class, 
                          levels = c("Actinomycetes","Gammaproteobacteria","Betaproteobacteria","Flavobacteriia",  
                                     "Bacilli","Alphaproteobacteria","Deltaproteobacteria","Clostridia","Desulfuromonadia",
                                     "Chlorobiia","Epsilonproteobacteria","Bacteroidia","Fusobacteriia",
                                     "Cytophagia","Negativicutes","Tissierellia","Erysipelotrichia","Unclassified"))
## Plot
ggplot(ARG_Class, aes(x = Sample_type, y = mean, fill = Class)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)") +
  scale_fill_manual(values = mycolors)
## Calculate others by summing  minimum arg
ARG_Class <-ARG_Class %>% 
  arrange(factor(Class, levels = c("Actinomycetes","Gammaproteobacteria","Betaproteobacteria","Flavobacteriia",  
                                   "Bacilli","Alphaproteobacteria","Deltaproteobacteria","Clostridia","Desulfuromonadia",
                                   "Chlorobiia","Epsilonproteobacteria","Bacteroidia","Fusobacteriia",
                                   "Cytophagia","Negativicutes","Tissierellia","Erysipelotrichia","Unclassified")))
other_class <- ARG_Class[25:34,] %>% group_by(Sample_type) %>% summarise(mean = sum(mean))
other_class$Class <- "Unclassified/Others"
ARG_Class <- ARG_Class[1:24,]
ARG_Class <- rbind(ARG_Class, other_class)
ARG_Class$Class <- factor(ARG_Class$Class, 
                          levels = c("Actinomycetes","Gammaproteobacteria","Betaproteobacteria","Flavobacteriia",  
                                     "Bacilli","Alphaproteobacteria","Deltaproteobacteria","Clostridia","Desulfuromonadia",
                                     "Chlorobiia","Epsilonproteobacteria","Unclassified/Others"))
#### Order sample_type
ARG_Class$Sample_type <- factor(ARG_Class$Sample_type, 
                                 levels = c("AT","ARP","ODP"))
## Plot
p <- ggplot(ARG_Class, aes(x = Sample_type, y = mean, fill = Class)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)") +
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
        legend.background = element_rect(fill='transparent')) + #transparent legend bg)
        scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                             "#FDB462","#B3DE69","#FCCDE5", "#FFED6F","#BC80BD",
                             "#CCEBC5","#D9D9D9"))
print(p)
# ggsave("ARG_coverage_class.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG_coverage",
#        width = 3.835, height = 5,
#        units = "in", bg='transparent') # save to png format


############### Order ################
ARG_Order <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Order) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Order$Sample_type <- gsub("1|2|3|4|5","",ARG_Order$Sample)
## Calculate mean
ARG_Order <- ARG_Order %>% ungroup() %>% group_by(Order, Sample_type) %>% mutate(mean = sum(coverage)/5)
## Filter necessary column for figure
ARG_Order <- ARG_Order %>% select(Order,Sample_type,mean) %>% unique()
## Expand color
nb.cols <- 46
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
# Order order tax
ARG_Order <- ARG_Order %>% arrange(desc(mean))
unique(ARG_Order$Order)
ARG_Order$Order <- factor(ARG_Order$Order, levels = c(  "Mycobacteriales" ,    "Burkholderiales"  ,   "Xanthomonadales"   , 
                                                        "Flavobacteriales",    "Pseudomonadales"  ,   "Moraxellales"      ,  "Bacillales"         ,
                                                        "Rhodocyclales"   ,    "Aeromonadales"    ,   "Enterobacterales"  ,  "Lactobacillales"    ,
                                                        "Micrococcales"   ,    "Kitasatosporales" ,   "Hyphomicrobiales"  ,  "Chromatiales"       ,
                                                        "Rhodospirillales",    "Nitrosomonadales" ,   "Rhodobacterales"   ,  "Propionibacteriales",
                                                        "Neisseriales"    ,    "Thiotrichales"    ,   "Eubacteriales"     ,  "Desulfovibrionales" ,
                                                        "Geobacterales"   ,    "Chlorobiales"     ,   "Sphingomonadales"  ,  "Pseudonocardiales"  ,
                                                        "Actinomycetales" ,    "Campylobacterales",   "Bacteroidales"     ,  "Marinilabiliales"   ,
                                                        "Fusobacteriales" ,    "Cytophagales"     ,   "Alteromonadales"   ,  "Pasteurellales"     ,
                                                        "Myxococcales"    ,    "Tissierellales"   ,   "Methylococcales"   ,  "Selenomonadales"    ,
                                                        "Vibrionales"     ,    "Hyphomonadales"   ,   "Veillonellales"    ,  "Erysipelotrichales" ,
                                                        "Unclassified"))
## Plot
ggplot(ARG_Order, aes(x = Sample_type, y = mean, fill = Order)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)")+
  scale_fill_manual(values = mycolors)
# Calculate proportion
order_proportion <- ARG_Order %>% group_by(Sample_type) %>% 
                                  mutate(proportion = (mean/sum(mean))*100)
## Calculate others by summing  minimum arg
unique(ARG_Order$Order)
ARG_Order <-ARG_Order %>% 
  arrange(factor(Order, levels = c("Mycobacteriales" ,    "Burkholderiales"  ,   "Xanthomonadales"   , 
                                   "Flavobacteriales",    "Pseudomonadales"  ,   "Moraxellales"      ,  "Bacillales"         ,
                                   "Rhodocyclales"   ,    "Aeromonadales"    ,   "Enterobacterales"  ,  "Lactobacillales"    ,
                                   "Micrococcales"   ,    "Kitasatosporales" ,   "Hyphomicrobiales"  ,  "Chromatiales"       ,
                                   "Rhodospirillales",    "Nitrosomonadales" ,   "Rhodobacterales"   ,  "Propionibacteriales",
                                   "Neisseriales"    ,    "Thiotrichales"    ,   "Eubacteriales"     ,  "Desulfovibrionales" ,
                                   "Geobacterales"   ,    "Chlorobiales"     ,   "Sphingomonadales"  ,  "Pseudonocardiales"  ,
                                   "Actinomycetales" ,    "Campylobacterales",   "Bacteroidales"     ,  "Marinilabiliales"   ,
                                   "Fusobacteriales" ,    "Cytophagales"     ,   "Alteromonadales"   ,  "Pasteurellales"     ,
                                   "Myxococcales"    ,    "Tissierellales"   ,   "Methylococcales"   ,  "Selenomonadales"    ,
                                   "Vibrionales"     ,    "Hyphomonadales"   ,   "Veillonellales"    ,  "Erysipelotrichales" ,
                                   "Unclassified")))
other_order <- ARG_Order[26:77,] %>% group_by(Sample_type) %>% summarise(mean = sum(mean))
other_order$Order <- "Unclassified/Others"
ARG_Order <- ARG_Order[1:25,]
ARG_Order <- rbind(ARG_Order, other_order)
## Plot
ARG_Order$Order <- factor(ARG_Order$Order, levels = c("Mycobacteriales" ,    "Burkholderiales"  ,   "Xanthomonadales"   , 
                                                      "Flavobacteriales",    "Pseudomonadales"  ,   "Moraxellales"      ,  "Bacillales"         ,
                                                      "Rhodocyclales"   ,    "Aeromonadales"    ,   "Enterobacterales"  ,  "Lactobacillales"    ,
                                                      "Unclassified/Others"))
#### Order sample_type
ARG_Order$Sample_type <- factor(ARG_Order$Sample_type, 
                                 levels = c("AT","ARP","ODP"))
p <- ggplot(ARG_Order, aes(x = Sample_type, y = mean, fill = Order)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)")+
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
        legend.background = element_rect(fill='transparent')) + #transparent legend bg)
        scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                             "#FDB462","#B3DE69","#FCCDE5", "#FFED6F","#BC80BD",
                             "#CCEBC5","#D9D9D9"))
print(p)
# ggsave("ARG_coverage_order.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG_coverage",
#        width = 3.65, height = 5,
#        units = "in", bg='transparent') # save to png format



############### Family ################
ARG_Family <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Family) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Family$Sample_type <- gsub("1|2|3|4|5","",ARG_Family$Sample)
## Calculate mean
ARG_Family <- ARG_Family %>% ungroup() %>% group_by(Family, Sample_type) %>% mutate(mean = sum(coverage)/5)
## Filter necesssary column for figure
ARG_Family <- ARG_Family %>% select(Family,Sample_type,mean) %>% unique()
## Expand color
nb.cols <- 93
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
# Order family tax
ARG_Family <- ARG_Family %>% arrange(desc(mean))
unique(ARG_Family$Family)
ARG_Family$Family <- factor(ARG_Family$Family, levels = c( "Mycobacteriaceae"     ,    "Xanthomonadaceae" ,    "Weeksellaceae"      ,   
                                                           "Comamonadaceae"       ,  "Pseudomonadaceae"   ,    "Moraxellaceae"      ,    "Azonexaceae"        ,   
                                                           "Gordoniaceae"         ,  "Aeromonadaceae"     ,    "Bacillaceae"        ,    "Enterobacteriaceae" ,   
                                                           "Microbacteriaceae"    ,  "Alcaligenaceae"     ,    "Streptomycetaceae"  ,    "Chromatiaceae"      ,   
                                                           "Staphylococcaceae"    ,  "Nitrosomonadaceae"  ,    "Micrococcaceae"     ,    "Nocardioidaceae"    ,   
                                                           "Enterococcaceae"      ,  "Chromobacteriaceae" ,    "Azospirillaceae"    ,    "Streptococcaceae"   ,   
                                                           "Thiotrichaceae"       ,  "Burkholderiaceae"   ,    "Sphaerotilaceae"    ,    "Brucellaceae"       ,   
                                                           "Flavobacteriaceae"    ,  "Paracoccaceae"      ,    "Desulfovibrionaceae",    "Casimicrobiaceae"   ,   
                                                           "Corynebacteriaceae"   ,  "Rhodospirillaceae"  ,    "Geobacteraceae"     ,    "Roseobacteraceae"   ,   
                                                           "Oxalobacteraceae"     ,  "Rhizobiaceae"       ,    "Chlorobiaceae"      ,    "Aerococcaceae"      ,   
                                                           "Morganellaceae"       ,  "Sphingomonadaceae"  ,    "Nocardiaceae"       ,    "Fluviibacteraceae"  ,   
                                                           "Clostridiaceae"       ,  "Pseudonocardiaceae" ,    "Hyphomicrobiaceae"  ,    "Intrasporangiaceae" ,   
                                                           "Actinomycetaceae"     ,  "Rhodocyclaceae"     ,    "Arcobacteraceae"    ,    "Acetobacteraceae"   ,   
                                                           "Erwiniaceae"          ,  "Pectobacteriaceae"  ,    "Marinilabiliaceae"  ,    "Rhodanobacteraceae" ,   
                                                           "Nitrobacteraceae"     ,  "Paenibacillaceae"   ,    "Sulfuricellaceae"   ,    "Lactobacillaceae"   ,   
                                                           "Fusobacteriaceae"     ,  "Hymenobacteraceae"  ,    "Devosiaceae"        ,    "Lachnospiraceae"    ,   
                                                           "Peptostreptococcaceae",  "Methylocystaceae"   ,    "Idiomarinaceae"     ,    "Carnobacteriaceae"  ,   
                                                           "Pasteurellaceae"      ,  "Oscillospiraceae"   ,    "Acidilutibacteraceae",   "Planococcaceae"     ,   
                                                           "Methylococcaceae"     ,  "Phyllobacteriaceae" ,    "Selenomonadaceae"   ,    "Vibrionaceae"       ,   
                                                           "Anaeromyxobacteraceae",  "Peptococcaceae"     ,    "Hyphomonadaceae"    ,    "Veillonellaceae"    ,   
                                                           "Bacteroidaceae"       ,  "Erysipelotrichaceae",    "Labilitrichaceae"   ,    "Pseudoalteromonadaceae",
                                                           "Unclassified"))
## Plot
ggplot(ARG_Family, aes(x = Sample_type, y = mean, fill = Family)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)")+
  scale_fill_manual(values = mycolors)
## Calculate others by summing  minimum arg
unique(ARG_Family$Family)
ARG_Family <-ARG_Family %>% 
  arrange(factor(Family, levels = c("Mycobacteriaceae"     ,    "Xanthomonadaceae" ,    "Weeksellaceae"      ,   
                                    "Comamonadaceae"       ,  "Pseudomonadaceae"   ,    "Moraxellaceae"      ,    "Azonexaceae"        ,   
                                    "Gordoniaceae"         ,  "Aeromonadaceae"     ,    "Bacillaceae"        ,    "Enterobacteriaceae" ,   
                                    "Microbacteriaceae"    ,  "Alcaligenaceae"     ,    "Streptomycetaceae"  ,    "Chromatiaceae"      ,   
                                    "Staphylococcaceae"    ,  "Nitrosomonadaceae"  ,    "Micrococcaceae"     ,    "Nocardioidaceae"    ,   
                                    "Enterococcaceae"      ,  "Chromobacteriaceae" ,    "Azospirillaceae"    ,    "Streptococcaceae"   ,   
                                    "Thiotrichaceae"       ,  "Burkholderiaceae"   ,    "Sphaerotilaceae"    ,    "Brucellaceae"       ,   
                                    "Flavobacteriaceae"    ,  "Paracoccaceae"      ,    "Desulfovibrionaceae",    "Casimicrobiaceae"   ,   
                                    "Corynebacteriaceae"   ,  "Rhodospirillaceae"  ,    "Geobacteraceae"     ,    "Roseobacteraceae"   ,   
                                    "Oxalobacteraceae"     ,  "Rhizobiaceae"       ,    "Chlorobiaceae"      ,    "Aerococcaceae"      ,   
                                    "Morganellaceae"       ,  "Sphingomonadaceae"  ,    "Nocardiaceae"       ,    "Fluviibacteraceae"  ,   
                                    "Clostridiaceae"       ,  "Pseudonocardiaceae" ,    "Hyphomicrobiaceae"  ,    "Intrasporangiaceae" ,   
                                    "Actinomycetaceae"     ,  "Rhodocyclaceae"     ,    "Arcobacteraceae"    ,    "Acetobacteraceae"   ,   
                                    "Erwiniaceae"          ,  "Pectobacteriaceae"  ,    "Marinilabiliaceae"  ,    "Rhodanobacteraceae" ,   
                                    "Nitrobacteraceae"     ,  "Paenibacillaceae"   ,    "Sulfuricellaceae"   ,    "Lactobacillaceae"   ,   
                                    "Fusobacteriaceae"     ,  "Hymenobacteraceae"  ,    "Devosiaceae"        ,    "Lachnospiraceae"    ,   
                                    "Peptostreptococcaceae",  "Methylocystaceae"   ,    "Idiomarinaceae"     ,    "Carnobacteriaceae"  ,   
                                    "Pasteurellaceae"      ,  "Oscillospiraceae"   ,    "Acidilutibacteraceae",   "Planococcaceae"     ,   
                                    "Methylococcaceae"     ,  "Phyllobacteriaceae" ,    "Selenomonadaceae"   ,    "Vibrionaceae"       ,   
                                    "Anaeromyxobacteraceae",  "Peptococcaceae"     ,    "Hyphomonadaceae"    ,    "Veillonellaceae"    ,   
                                    "Bacteroidaceae"       ,  "Erysipelotrichaceae",    "Labilitrichaceae"   ,    "Pseudoalteromonadaceae",
                                    "Unclassified")))
other_family <- ARG_Family[23:126,] %>% group_by(Sample_type) %>% summarise(mean = sum(mean))
other_family$Family <- "Unclassified/Others"
ARG_Family <- ARG_Family[1:22,]
ARG_Family <- rbind(ARG_Family, other_family)
## Plot
ARG_Family$Family <- factor(ARG_Family$Family, levels = c("Mycobacteriaceae"     ,    "Xanthomonadaceae" ,    "Weeksellaceae"      ,   
                                                          "Comamonadaceae"       ,  "Pseudomonadaceae"   ,    "Moraxellaceae"      ,    "Azonexaceae"        ,   
                                                          "Gordoniaceae"         ,  "Aeromonadaceae"     ,    "Bacillaceae"        ,    "Enterobacteriaceae" ,
                                                          "Unclassified/Others"))
ARG_Family$Sample_type <- factor(ARG_Family$Sample_type, 
                                levels = c("AT","ARP","ODP"))
p <- ggplot(ARG_Family, aes(x = Sample_type, y = mean, fill = Family)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)") +
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
        legend.background = element_rect(fill='transparent')) + #transparent legend bg)
        scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                             "#FDB462","#B3DE69","#FCCDE5", "#FFED6F","#BC80BD",
                             "#CCEBC5","#D9D9D9"))
print(p)
# ggsave("ARG_coverage_family.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG_coverage",
#        width = 3.72, height = 5,
#        units = "in", bg='transparent') # save to png format



############### Genus ################
ARG_Genus <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Genus) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Genus$Sample_type <- gsub("1|2|3|4|5","",ARG_Genus$Sample)
## Calculate mean
ARG_Genus <- ARG_Genus %>% ungroup() %>% group_by(Genus, Sample_type) %>% mutate(mean = sum(coverage)/5)
## Filter necesssary column for figure
ARG_Genus <- ARG_Genus %>% select(Genus,Sample_type,mean) %>% unique()
## Expand color
nb.cols <- 178
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
## Plot
ggplot(ARG_Genus, aes(x = Sample_type, y = mean, fill = Genus)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)") +
  scale_fill_manual(values = mycolors)
## Calculate others by summing  minimum arg
ARG_Genus <- ARG_Genus %>% arrange(desc(mean))
unique(ARG_Genus$Genus)
ARG_Genus <-ARG_Genus %>% 
  arrange(factor(Genus,  levels = c( "Mycolicibacterium", "Pseudomonas"     ,  "Xanthomonas"      , "Acinetobacter"   ,
                                     "Mycobacterium"    , "Gordonia"        ,  "Aeromonas"        , "Chryseobacterium",  "Dechloromonas"    ,
                                     "Acidovorax"       , "Lysobacter"      ,  "Lysinibacillus"   , "Nitrosococcus"   ,  "Ferribacterium"   ,
                                     "Empedobacter"     , "Riemerella"      ,  "Ottowia"          , "Stenotrophomonas",  "Klebsiella"       ,
                                     "Pigmentiphaga"    , "Nitrosomonas"    ,  "Variovorax"       , "Nocardioides"    ,  "Streptomyces"     ,
                                     "Macrococcus"      , "Cloacibacterium" ,  "Niveispirillum"   , "Hydrogenophaga"  ,  "Actinacidiphila"  ,
                                     "Streptococcus"    , "Chromobacterium" ,  "Thiothrix"        , "Staphylococcus"  ,  "Desulfovibrio"    ,
                                     "Casimicrobium"    , "Enterococcus"    ,  "Paenarthrobacter" , "Gemmobacter"     ,  "Corynebacterium"  ,
                                     "Brucella"         , "Burkholderia"    ,  "Comamonas"        , "Trichlorobacter" ,  "Bacillus"         ,
                                     "Defluviicoccus"   , "Achromobacter"   ,  "Shinella"         , "Chlorobaculum"   ,  "Aerococcus"       ,
                                     "Flavobacterium"   , "Massilia"        ,  "Sulfitobacter"    , "Vagococcus"      ,  "Morganella"       ,
                                     "Rubrivivax"       , "Pulveribacter"   ,  "Kaistella"        , "Denitrificimonas",  "Schlegelella"     ,
                                     "Enterobacter"     , "Glutamicibacter" ,  "Rhodococcus"      , "Escherichia"     ,  "Fluviibacter"     ,
                                     "Pseudonocardia"   , "Microbacterium"  ,  "Stutzerimonas"    , "Bordetella"      ,  "Hyphomicrobium"   ,
                                     "Janibacter"       , "Paraclostridium" ,  "Trueperella"      , "Aliarcobacter"   ,  "Rhodoferax"       ,
                                     "Raoultella"       , "Elizabethkingia" ,  "Aromatoleum"      , "Myroides"        ,  "Roseomonas"       ,
                                     "Citricoccus"      , "Salmonella"      ,  "Pectobacterium"   , "Melaminivora"    ,  "Protaetiibacter"  ,
                                     "Alkalitalea"      , "Arenimonas"      ,  "Thermomonas"      , "Bradyrhizobium"  ,  "Massilia group"   ,
                                     "Paenibacillus"    , "Diaphorobacter"  ,  "Rhizorhabdus"     , "Sulfurimicrobium",  "Paludibacterium"  ,
                                     "Sphingomonas"     , "Moraxella"       ,  "Fusobacterium"    , "Erwinia"         ,  "Pontibacter"      ,
                                     "Cupriavidus"      , "Cedecea"         ,  "Paradevosia"      , "Microcella"      ,  "Lacrimispora"     ,
                                     "Tahibacter"       , "Romboutsia"      ,  "Metabacillus"     , "Mesobacillus"    ,  "Magnetospirillum" ,
                                     "Methylosinus"     , "Paraburkholderia",  "Roseburia"        , "Idiomarina"      ,  "Jeotgalibaca"     ,
                                     "Actinobacillus"   , "Sphingopyxis"    ,  "Pasteurella"      , "Acidilutibacter" ,  "Rathayibacter"    ,
                                     "Ralstonia"        , "Amniculibacterium", "Planococcus"      , "Collimonas"      ,  "Methylomonas"     ,
                                     "Celeribacter"     , "Peribacillus"    ,  "Mesorhizobium"    , "Nordella"        ,  "Megamonas"        ,
                                     "Vibrio"           , "Pseudoxanthomonas", "Maribacter"       , "Anaeromyxobacter",  "Pantoea"          ,
                                     "Clostridium"      , "Quatrionicoccus" ,  "Luteibacter"      , "Acetivibrio"     ,  "Pseudoduganella"  ,
                                     "Thioflavicoccus"  , "Phycicoccus"     ,  "Desulfofarcimen"  , "Haematobacter"   ,  "Hyphomonas"       ,
                                     "Megasphaera"      , "Xylophilus"      ,  "Bacteroides"      , "Erysipelothrix"  ,  "Kitasatospora"    ,
                                     "Novosphingobium"  , "Nitrobacter"     ,  "Flavonifractor"   , "Paracoccus"      ,  "Azospira"         ,
                                     "Labilithrix"      , "Pseudoalteromonas", "Unclassified"   )))
other_genus <- ARG_Genus[23:204,] %>% group_by(Sample_type) %>% summarise(mean = sum(mean))
other_genus$Genus <- "Unclassified/Others"
ARG_Genus <- ARG_Genus[1:22,]
ARG_Genus <- rbind(ARG_Genus, other_genus)
# Plot
ARG_Genus$Genus <- factor(ARG_Genus$Genus, levels = c("Mycolicibacterium", "Pseudomonas"     ,  "Xanthomonas"      , "Acinetobacter"   ,
                                                      "Mycobacterium"    , "Gordonia"        ,  "Aeromonas"        , "Chryseobacterium",  "Dechloromonas"    ,
                                                      "Acidovorax"       , "Lysobacter"      ,  "Unclassified/Others"))
ARG_Genus$Sample_type <- factor(ARG_Genus$Sample_type, 
                                 levels = c("AT","ARP","ODP"))
p <- ggplot(ARG_Genus, aes(x = Sample_type, y = mean, fill = Genus)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)") +
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
        legend.background = element_rect(fill='transparent')) + #transparent legend bg)
        scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                             "#FDB462","#B3DE69","#FCCDE5", "#FFED6F","#BC80BD",
                             "#CCEBC5","#D9D9D9"))
print(p)
# ggsave("ARG_coverage_genus.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG_coverage",
#        width = 3.64, height = 5,
#        units = "in", bg='transparent') # save to png format


######### Species ##########
ARG_Species <- final_ARG_coverage  %>% ungroup() %>% group_by(Sample, Species) %>% 
  summarise(coverage = sum(coverage)) %>% filter(!(coverage == 0))
ARG_Species$Sample_type <- gsub("1|2|3|4|5","",ARG_Species$Sample)
## Calculate mean
ARG_Species <- ARG_Species %>% ungroup() %>% group_by(Species, Sample_type) %>% mutate(mean = sum(coverage)/5)
## Filter necessary column for figure
ARG_Species <- ARG_Species %>% select(Species,Sample_type,mean) %>% unique()
## Expand color
nb.cols <- 434
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
## Plot
ggplot(ARG_Species, aes(x = Sample_type, y = mean, fill = Species)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)") +
  scale_fill_manual(values = mycolors)
## Calculate others by summing  minimum arg
ARG_Species <- ARG_Species %>% arrange(desc(mean))
unique(ARG_Species$Species)
ARG_Species <-ARG_Species %>% 
  arrange(factor(Species, levels = c( "Mycolicibacterium insubricum"                ,  
                                      "Pseudomonas aeruginosa"                      , "Mycolicibacterium confluentis"             , 
                                      "Mycolicibacterium neoaurum"                  , "Gordonia amarae"                           ,  
                                      "Xanthomonas euroxanthea"                     , "Aeromonas sp. ASNIH1"                      ,  
                                      "Xanthomonas translucens pv. Phleipratensis"  , "Mycobacterium sp. DL592"                   ,  
                                      "Chryseobacterium sp. KACC 21268"             , "Nitrosococcus wardiae"                     ,  
                                      "Ferribacterium limneticum"                   , "Riemerella anatipestifer"                  ,  
                                      "Mycolicibacterium vaccae 95051"              , "Pigmentiphaga aceris"                      ,  
                                      "Lysinibacillus capsici"                      , "Gordonia insulae"                          ,  
                                      "Acinetobacter schindleri"                    , "Dechloromonas denitrificans"               ,  
                                      "Mycolicibacterium anyangense"                , "Dechloromonas sp. HYN0024"                 ,  
                                      "Stenotrophomonas maltophilia"                , "Acidovorax sp. YS12"                       ,  
                                      "Nitrosomonas eutropha C91"                   , "Empedobacter brevis"                       , 
                                      "Klebsiella quasipneumoniae"                  , "Mycobacterium sp. MS1601"                  ,  
                                      "Mycolicibacterium phocaicum"                 , "Chryseobacterium sp. POL2"                 ,  
                                      "Niveispirillum cyanobacteriorum"             , "Mycolicibacterium fallax"                  ,  
                                      "Actinacidiphila bryophytorum"                , "Ottowia testudinis"                        , 
                                      "Acinetobacter variabilis"                    , "Mycolicibacterium brumae"                  ,  
                                      "Mycolicibacterium chitae"                    , "Thiothrix litoralis"                       ,  
                                      "Desulfovibrio vulgaris"                      , "Casimicrobium huifangae"                   ,  
                                      "Paenarthrobacter ureafaciens"                , "Gemmobacter aquarius"                      ,  
                                      "Lysinibacillus sp. 2017"                     , "Hydrogenophaga sp. RAC07"                  ,  
                                      "Chromobacterium rhizoryzae"                  , "Cloacibacterium normanense"                ,  
                                      "Ottowia oryzae"                              , "Acinetobacter sp. TR3"                     ,  
                                      "Brucella anthropi"                           , "Acinetobacter baumannii"                   ,  
                                      "Chryseobacterium sp. Y16C"                   , "Mycolicibacterium rhodesiae NBB3"          ,  
                                      "Acinetobacter sp. C32I"                      , "Burkholderia cepacia complex"              ,  
                                      "Trichlorobacter lovleyi SZ"                  , "Chryseobacterium manosquense"              ,  
                                      "Comamonas sp. NLF-7-7"                       , "Defluviicoccus vanus"                      ,  
                                      "Empedobacter stercoris"                      , "Mycolicibacterium helvum"                  ,  
                                      "Chryseobacterium sp."                        , "Mycobacterium sp. IDR2000157661"           ,  
                                      "Shinella zoogloeoides"                       , "Chlorobaculum parvum NCIB 8327"            ,  
                                      "Aerococcus urinaeequi"                       , "Nocardioides daphniae"                     ,  
                                      "Flavobacterium branchiophilum FL-15"         , "Mycolicibacterium sarraceniae"             ,  
                                      "Sulfitobacter sp. BSw21498"                  , "Xanthomonas arboricola"                    ,  
                                      "Klebsiella pneumoniae"                       , "Mycolicibacterium litorale"                , 
                                      "Morganella morganii"                         , "Rubrivivax gelatinosus IL144"              , 
                                      "Pulveribacter suum"                          , "Dechloromonas sp. TW-R-39-2"               , 
                                      "Acidovorax radicis"                          , "Riemerella anatipestifer Yb2"              ,  
                                      "Chryseobacterium taklimakanense"             , "Kaistella sp. BT6-1-3"                     ,  
                                      "Denitrificimonas caeni"                      , "Staphylococcus aureus"                     ,  
                                      "Schlegelella aquatica"                       , "Mycolicibacterium fluoranthenivorans"      ,  
                                      "Acidovorax sp. JMULE5"                       , "Acinetobacter towneri"                     ,  
                                      "Bacillus anthracis"                          , "Enterobacter hormaechei"                   ,  
                                      "Mycolicibacterium poriferae"                 , "Aeromonas veronii"                         ,  
                                      "Rhodococcus coprophilus"                     , "Escherichia coli"                          ,  
                                      "Cloacibacterium caeni"                       , "Fluviibacter phosphoraccumulans"           ,  
                                      "Lysobacter sp. CW239"                        , "Pseudomonas chlororaphis subsp. Aureofaciens",
                                      "Pseudonocardia sp. HH130630-07"              , "Acidovorax sp. HDW3"                       ,  
                                      "Hyphomicrobium denitrificans ATCC 51888"     , "Stenotrophomonas acidaminiphila"           ,  
                                      "Nocardioides sp. JS614"                      , "Staphylococcus pseudoxylosus"              ,  
                                      "Lysinibacillus sp. CD3-6"                    , "Lysinibacillus sp. G01H"                   ,  
                                      "Paraclostridium bifermentans"                , "Stenotrophomonas sp. 169"                  ,  
                                      "Trueperella pyogenes"                        , "Pseudomonas sp. GCEP-101"                  ,  
                                      "Corynebacterium diphtheriae"                 , "Variovorax sp. PAMC28562"                  ,  
                                      "Bordetella trematum"                         , "Hydrogenophaga sp. NH-16"                  ,  
                                      "Rhodoferax sp. BAB1"                         , "Enterococcus faecalis"                     , 
                                      "Variovorax paradoxus"                        , "Nitrosomonas ureae"                        ,  
                                      "Raoultella ornithinolytica"                  , "Glutamicibacter nicotianae"                ,  
                                      "Elizabethkingia anophelis"                   , "Mycobacterium heraklionense"               ,  
                                      "Achromobacter xylosoxidans"                  , "Aromatoleum petrolei"                      , 
                                      "Nocardioides sp. LMS-CY"                     , "Acidovorax sp. RAC01"                      , 
                                      "Myroides odoratimimus"                       , "Acinetobacter sp. NEB 394"                 ,  
                                      "Citricoccus sp. NR2"                         , "Chromobacterium vaccinii"                  ,  
                                      "Nocardioides ungokensis"                     , "Variovorax sp. HW608"                      ,  
                                      "Macrococcus armenti"                         , "Stutzerimonas stutzeri group"              ,  
                                      "Thiothrix winogradskyi"                      , "Melaminivora suipulveris"                  ,  
                                      "Hydrogenophaga sp. PBL-H3"                   , "Protaetiibacter intestinalis"              ,  
                                      "Mycolicibacterium tokaiense"                 , "Massilia sp. H6"                           ,  
                                      "Corynebacterium freneyi"                     , "Microbacterium esteraromaticum"            ,  
                                      "Alkalitalea saponilacus"                     , "Arenimonas daejeonensis"                   ,  
                                      "Vagococcus lutrae"                           , "Enterococcus faecium"                      ,  
                                      "Mycolicibacterium psychrotolerans"           , "Streptomyces genisteinicus"                ,  
                                      "Lysobacter enzymogenes"                      , "Variovorax sp. RKNM96"                     , 
                                      "Mycobacterium grossiae"                      , "Paenibacillus urinalis"                    ,  
                                      "Massilia sp. YMA4"                           , "Acidovorax sp. 1608163"                    ,  
                                      "Rhizorhabdus wittichii"                      , "Sulfurimicrobium lacus"                    ,  
                                      "Paludibacterium sp. B53371"                  , "Mycolicibacterium gilvum Spyr1"            ,  
                                      "Mycolicibacterium sediminis"                 , "Lysobacter sp. S4-A87"                     ,  
                                      "Moraxella osloensis"                         , "Fusobacterium polymorphum"                 ,  
                                      "Burkholderia sp. FERM BP-3421"               , "Streptomyces actuosus"                     ,  
                                      "Pseudomonas sp. L5B5"                        , "Streptococcus pneumoniae"                  ,  
                                      "Erwinia sp. J780"                            , "Aliarcobacter cryaerophilus"               ,  
                                      "Diaphorobacter sp. JS3050"                   , "Mycolicibacterium madagascariense"         ,  
                                      "Acidovorax avenae"                           , "Pontibacter pudoricolor"                   ,  
                                      "Cupriavidus pauculus"                        , "Cedecea neteri"                            ,  
                                      "Paradevosia shaoguanensis"                   , "Diaphorobacter sp. ED-3"                   ,  
                                      "Thermomonas brevis"                          , "Variovorax sp. PAMC26660"                  ,  
                                      "Chryseobacterium suipulveris"                , "Gordonia pseudamarae"                      ,  
                                      "Bradyrhizobium diazoefficiens"               , "Acinetobacter indicus"                     ,  
                                      "Tahibacter sp. W38"                          , "Romboutsia sp. CE17"                       ,  
                                      "Vagococcus carniphilus"                      , "Microbacterium foliorum"                   ,  
                                      "Metabacillus sp. B2-18"                      ,"Corynebacterium jeikeium"                   , 
                                      "Sphingomonas sp. CL5.1"                      ,"Vagococcus fluvialis"                       , 
                                      "Pseudomonas pohangensis"                     , "Hydrogenophaga taeniospiralis"             ,  
                                      "Mesobacillus jeotgali"                       , "Rhodococcus sp. PBTS 1"                    ,  
                                      "Magnetospirillum magneticum AMB-1"           , "Gordonia sp. PP30"                         ,  
                                      "Lysobacter soli"                             , "Methylosinus trichosporium OB3b"           ,  
                                      "Paraburkholderia xenovorans LB400"           , "Glutamicibacter creatinolyticus"           ,  
                                      "Klebsiella pneumoniae MGH 39"                , "Roseburia intestinalis XB6B4"              ,  
                                      "Lysobacter sp. TY2-98"                       , "Klebsiella variicola subsp. Variicola"     ,  
                                      "Jeotgalibaca porci"                          , "Achromobacter xylosoxidans A8"             ,  
                                      "Actinobacillus pleuropneumoniae"             , "Pseudomonas sp. WJP1"                      ,  
                                      "Streptomyces albus"                          , "Lysobacter antibioticus"                   ,  
                                      "Stutzerimonas stutzeri"                      , "Sphingopyxis macrogoltabida"               ,  
                                      "Acinetobacter gyllenbergii"                  , "Pasteurella multocida"                     ,  
                                      "Streptomyces sp. NHF165"                     , "Aliarcobacter cryaerophilus D2610"         ,  
                                      "Enterococcus faecium DO"                     , "Enterococcus cecorum"                      ,  
                                      "Corynebacterium vitaeruminis DSM 20294"      , "Pseudomonas sp. J380"                      ,  
                                      "Mycolicibacterium celeriflavum"              , "Acidovorax sp. KKS102"                     ,  
                                      "Streptomyces sp. SirexAA-E"                  , "Mycolicibacterium fortuitum"               ,  
                                      "Mycolicibacterium mucogenicum DSM 44124"     , "Mycolicibacterium arabiense"               ,  
                                      "Mycolicibacterium moriokaense"               , "Thermomonas aquatica"                      ,  
                                      "Corynebacterium xerosis"                     , "Clostridium sphenoides JCM 1415"           ,  
                                      "Acidilutibacter cellobiosedens"              , "Aerococcus viridans"                       ,    
                                      "Nocardioides cynanchi"                       , "Streptomyces tubbatahanensis"              ,  
                                      "Rathayibacter sp. VKM Ac-2759"               , "Ralstonia mannitolilytica"                 ,    
                                      "Streptomyces durmitorensis"                  , "Amniculibacterium sp. G2-70"               ,      
                                      "Lysobacter oculi"                            , "Streptomyces sp. INR7"                     ,     
                                      "Mycobacterium haemophilum DSM 44634"         , "Achromobacter sp. MFA1 R4"                 ,  
                                      "Enterococcus sp. FDAARGOS_375"               , "Thermomonas sp. XSG"                       , 
                                      "Empedobacter falsenii"                       , "Planococcus sp. 107-1"                     ,  
                                      "Burkholderia sp. Bp7605"                     , "Acidovorax sp. NCPPB 2350"                 ,  
                                      "Collimonas arenae"                           , "Methylomonas koyamae"                      ,  
                                      "Celeribacter indicus"                        , "Lysobacter capsici"                        ,  
                                      "Nocardioides sambongensis"                   , "Staphylococcus saprophyticus"              ,  
                                      "Peribacillus psychrosaccharolyticus"         , "Nordella sp. HKS 07"                       ,  
                                      "Megamonas funiformis"                        , "Acinetobacter sp. SCLZS86"                 ,  
                                      "Mycobacterium noviomagense"                  , "Chryseobacterium sp. 6424"                 ,  
                                      "Pseudoxanthomonas sp. SL93"                  , "Pseudomonas sp. Q1-7"                      ,  
                                      "Maribacter sp. Hal144"                       , "Chryseobacterium sp. MEBOG07"              ,  
                                      "Bordetella genomosp. 13"                     , "Streptomyces bathyalis"                    ,  
                                      "Pseudomonas fluorescens R124"                , "Anaeromyxobacter oryzae"                   ,  
                                      "Pantoea rwandensis"                          , "Acinetobacter johnsonii"                   ,  
                                      "Clostridium perfringens"                     , "Burkholderia multivorans"                  ,  
                                      "Quatrionicoccus australiensis"               , "Corynebacterium variabile DSM 44702"       ,  
                                      "Luteibacter anthropi"                        , "Acetivibrio thermocellus ATCC 27405"       ,  
                                      "Pseudoduganella flava"                       , "Comamonas sp. NLF-1-9"                     ,  
                                      "Achromobacter sp. AONIH1"                    , "Thioflavicoccus mobilis 8321"              ,  
                                      "Phycicoccus endophyticus"                    , "Thermomonas sp. IMCC34681"                 ,  
                                      "Mycolicibacterium austroafricanum"           , "Aeromonas sp. FDAARGOS 1410"               ,  
                                      "Mycolicibacterium mengxianglii"              , "Desulfofarcimen acetoxidans DSM 771"       ,  
                                      "Haematobacter massiliensis"                  , "Bradyrhizobium sp. 1(2017)"                ,  
                                      "Acidovorax monticola"                        , "Lacrimispora saccharolytica"               ,  
                                      "Burkholderia plantarii"                      , "Roseomonas sp. FDAARGOS_362"               ,  
                                      "Streptococcus suis"                          , "Hyphomonas atlantica"                      ,  
                                      "Massilia forsythiae"                         , "Megasphaera stantonii"                     ,  
                                      "Xylophilus rhododendri"                      , "Comamonas aquatica"                        ,  
                                      "Bacteroides faecium"                         , "Staphylococcus xylosus"                    ,  
                                      "Pseudomonas siliginis"                       , "Erysipelothrix rhusiopathiae"              ,  
                                      "Nocardioides sp. S5"                         , "Kitasatospora sp. MMS16-BH015"             ,  
                                      "Novosphingobium sp. KACC 22771"              , "Massilia sp. NP310"                        ,  
                                      "Nitrobacter hamburgensis X14"                , "Burkholderia gladioli"                     ,  
                                      "Pseudomonas muyukensis"                      , "Flavonifractor plautii"                    ,  
                                      "Paracoccus aminophilus JCM 7686"             , "Mycolicibacterium rutilum"                 ,  
                                      "Streptococcus oriscaviae"                    , "Diaphorobacter sp. HDW4B"                  ,  
                                      "Sphingomonas sp. S2-65"                      , "Azospira restricta"                        , 
                                      "Labilithrix luteola"                         , "Sphingomonas profundi"                     , 
                                      "Bradyrhizobium sp. CCBAU 051011"             , "Pseudomonas orientalis"                    ,  
                                      "Pseudoalteromonas luteoviolacea"             , "Cupriavidus necator"                       ,
                                      "Unclassified"                                )))
other_species <- ARG_Species[17:364,] %>% group_by(Sample_type) %>% summarise(mean = sum(mean))
other_species$Species <- "Unclassified/Others"
ARG_Species <- ARG_Species[1:16,]
ARG_Species <- rbind(ARG_Species, other_species)
# Plot
ARG_Species$Species <- factor(ARG_Species$Species, levels = c("Mycolicibacterium insubricum"                ,  
                                                              "Pseudomonas aeruginosa"                      , "Mycolicibacterium confluentis"             , 
                                                              "Mycolicibacterium neoaurum"                  , "Gordonia amarae"                           ,  
                                                              "Xanthomonas euroxanthea"                     , "Aeromonas sp. ASNIH1"                      ,  
                                                              "Xanthomonas translucens pv. Phleipratensis"  , "Mycobacterium sp. DL592"                   ,  
                                                              "Chryseobacterium sp. KACC 21268"             , "Nitrosococcus wardiae"                     ,
                                                              "Unclassified/Others"))
ARG_Species$Sample_type <- factor(ARG_Species$Sample_type, 
                                levels = c("AT","ARP","ODP"))
p <- ggplot(ARG_Species, aes(x = Sample_type, y = mean, fill = Species)) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)") +
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
        legend.background = element_rect(fill='transparent')) + #transparent legend bg)
        scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                             "#FDB462","#B3DE69","#FCCDE5", "#FFED6F","#BC80BD",
                             "#CCEBC5","#D9D9D9"))
print(p)
# ggsave("ARG_coverage_species.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG_coverage",
#        width = 5.32, height = 5,
#        units = "in", bg='transparent') # save to png format

