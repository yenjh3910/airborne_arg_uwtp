# bin_abundance.R

# Import
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
bin_abundance <- read.table("../../airborne_arg_uwtp_result/metawrap_bin/bin_quant/bin_abundance_table.tab",
                          header = TRUE, sep="\t")
colnames(bin_abundance)[1] <- "MAG"
source("./bin_gtdbtk.R")
colnames(sep_bin_tax)[1] <- "MAG"
source("./bins_ARG_diamond.R")
source("./bins_MGE_diamond.R")
source("./bins_VF_diamond.R")

# Filter ARG carrying bin
SARG_MAG <- bin_abundance %>% filter(MAG %in% SARG$BinID)
# Join bin with gtdbtk tax name
SARG_MAG <- left_join(SARG_MAG, sep_bin_tax)
SARG_MAG_sp <- SARG_MAG %>% filter(!(Species=='s__')) %>% 
                            select(!(Domain)) %>% select(!(Phylum)) %>%
                            select(!(Class)) %>% select(!(Order)) %>%
                            select(!(Family)) %>% select(!(Genus))
SARG_MAG_g <- SARG_MAG %>% filter(Species=='s__') %>% 
                           filter(!(Genus=='g__')) %>% 
                           select(!(Domain)) %>% select(!(Phylum)) %>%
                           select(!(Class)) %>% select(!(Order)) %>%
                           select(!(Family)) %>% select(!(Species))
SARG_MAG_fam <- SARG_MAG %>% filter(Species=='s__') %>% 
                             filter(Genus=='g__') %>% 
                             filter(!(Family=='f__')) %>% 
                             select(!(Domain)) %>% select(!(Phylum)) %>%
                             select(!(Class)) %>% select(!(Order)) %>%
                             select(!(Genus)) %>% select(!(Species))
colnames(SARG_MAG_sp)[12] <- "Taxonomy"
colnames(SARG_MAG_g)[12] <- "Taxonomy"
colnames(SARG_MAG_fam)[12] <- "Taxonomy"
SARG_MAG <- rbind(SARG_MAG_sp,SARG_MAG_g)
SARG_MAG <- rbind(SARG_MAG,SARG_MAG_fam)
colnames(SARG)[12] <- "MAG" # View later
SARG <- left_join(SARG,SARG_MAG) # View later
SARG <- SARG %>% select(MAG,subtype,type,Taxonomy) # View later
SARG$Taxonomy <- paste(SARG$MAG, SARG$Taxonomy, sep="_") # View later
SARG_MAG$Taxonomy <- paste(SARG_MAG$MAG, SARG_MAG$Taxonomy, sep="_")

# Annotation row
# Tax
tax_annontate <- SARG_MAG %>% select(Taxonomy)
colnames(tax_annontate) <- "tmp"
species_annotate <- tax_annontate %>% filter(grepl('s__',tax_annontate$tmp)) %>% 
                                      mutate(Taxonomy = 'Species')
genus_annotate  <- tax_annontate %>% filter(grepl('g__',tax_annontate$tmp)) %>% 
                                    mutate(Taxonomy = 'Genus')
family_annotate <- tax_annontate %>% filter(grepl('f__',tax_annontate$tmp)) %>% 
                                     mutate(Taxonomy = 'Family')
tax_annontate <- rbind(species_annotate,genus_annotate)
tax_annontate <- rbind(tax_annontate,family_annotate)
# Gene 
# ARG+MGE+VF
ARG_MGE_VF_MAG <- bin_abundance %>% filter(MAG %in% SARG$MAG) %>% 
  filter(MAG %in% MGE$BinID) %>% 
  filter(MAG %in% VF$BinID) %>% select(MAG)
# ARG+VF
ARG_VF_MAG <- bin_abundance %>% filter(MAG %in% SARG$MAG) %>% 
  filter(MAG %in% VF$BinID) %>% 
  filter(!(MAG %in% MGE$BinID)) %>% select(MAG)
# ARG+MGE
ARG_MGE_MAG <- bin_abundance %>% filter(MAG %in% SARG$MAG) %>% 
  filter(MAG %in% MGE$BinID) %>% 
  filter(!(MAG %in% VF$BinID)) %>% select(MAG)
## Join
ARG_MGE_VF_MAG <- left_join(ARG_MGE_VF_MAG,SARG_MAG) %>% select(MAG,Taxonomy)
row.names(ARG_MGE_VF_MAG) <- ARG_MGE_VF_MAG$Taxonomy
ARG_MGE_VF_MAG <- ARG_MGE_VF_MAG %>% mutate(Genes = "ARGs+MGEs+VFs")
ARG_MGE_VF_MAG <- ARG_MGE_VF_MAG %>% select(Genes)

ARG_VF_MAG <- left_join(ARG_VF_MAG,SARG_MAG) %>% select(MAG,Taxonomy)
row.names(ARG_VF_MAG) <- ARG_VF_MAG$Taxonomy
ARG_VF_MAG <- ARG_VF_MAG %>% mutate(Genes = "ARGs+VFs")
ARG_VF_MAG <- ARG_VF_MAG %>% select(Genes)

ARG_MGE_MAG <- left_join(ARG_MGE_MAG,SARG_MAG) %>% select(MAG,Taxonomy)
row.names(ARG_MGE_MAG) <- ARG_MGE_MAG$Taxonomy
ARG_MGE_MAG <- ARG_MGE_MAG %>% mutate(Genes = "ARGs+MGEs")
ARG_MGE_MAG <- ARG_MGE_MAG %>% select(Genes)

annotation_row <- rbind(ARG_MGE_VF_MAG,ARG_VF_MAG)
annotation_row <- rbind(annotation_row,ARG_MGE_MAG)
# Join tax and gene
annotation_row$tmp <- row.names(annotation_row)
annotation_row <- full_join(tax_annontate,annotation_row)
row.names(annotation_row) <- annotation_row[,1]
annotation_row <- annotation_row %>% select(!(tmp))

# taxa as column name
SARG_MAG <- SARG_MAG[,-1]
row.names(SARG_MAG) <- SARG_MAG[,11]
SARG_MAG <- SARG_MAG[,-11]
# Log scale
log_SARG_MAG <- log10(SARG_MAG)
log_SARG_MAG[log_SARG_MAG == -Inf] <- NA

# Annotation col
annotation_col = data.frame(Sample = factor(rep(c("AT", "ARP"), 
                                                c(5, 5))))
rownames(annotation_col) = colnames(log_SARG_MAG)
# Select color
display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(12, "Set3")
ann_colors = list(Sample = c(AT = "#fab0ab", ARP = "#76d9dd"),
                  Taxonomy = c(Family = "#FFED6F", Genus = "#FCCDE5", Species ="#CCEBC5"),
                  Genes = c('ARGs+MGEs+VFs' ="#8DD3C7",'ARGs+MGEs' ="#FDB462",'ARGs+VFs' ="#BEBADA"))
# Remove taxa_prefix
row.names(annotation_row) <- gsub(".__","",row.names(annotation_row))
row.names(annotation_row) <- gsub("bin.","Bin ",row.names(annotation_row))
row.names(log_SARG_MAG) <- gsub(".__","",row.names(log_SARG_MAG))
row.names(log_SARG_MAG) <- gsub("bin.","Bin ",row.names(log_SARG_MAG))
# Plot
p <-pheatmap(log_SARG_MAG, 
         annotation_col = annotation_col, annotation_row = annotation_row,
         annotation_colors = ann_colors, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         fontsize = 12, fontsize_row = 12, fontsize_col = 10,
         cellwidth = 12, cellheight = 17, bg = "transparent")
print(p)

# ggsave("MAG_heatmap.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/binning",
#        width = 8.5, height = 7,
#        units = "in", bg='transparent') # save to png format


# View carrying ARG
MAG_name <- SARG %>% select(MAG,Taxonomy) %>% unique()
colnames(MGE)[12] <- "MAG"
MGE <- MGE %>% select(MAG,subtype,type)
MGE <- left_join(MGE,MAG_name)
MGE <- MGE[complete.cases(MGE), ]
arg_mge_bin <- rbind(SARG,MGE)

# GGGene
# Import manual csv
library(openxlsx)
library(gggenes)
bin_gene <- read.xlsx("../../airborne_arg_uwtp_result/bins_diamond/mag_gene_position.xlsx")

ggplot(bin_gene, aes(xmin = adjust_start, xmax = adjust_end, y = MAG, fill = Type,
                         forward = Direction, label = Subtype)) +
  geom_gene_arrow() +
  geom_gene_label(align = "left") +
  facet_wrap(~ MAG, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

ggplot(bin_gene, aes(xmin = adjust_start, xmax = adjust_end, y = MAG, fill = Type,
                     forward = Direction)) +
  geom_gene_arrow() +
  facet_wrap(~ MAG, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  theme(legend.position = "top") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))



# See VF+ARG bin
VF_SARG <- VF %>% filter(BinID %in% SARG$MAG)
# VF number
VF_bin <- VF_SARG$BinID %>% unique() 
for (bin in VF_bin){ 
tmp <- VF_SARG %>% filter(BinID == bin) %>% select(gene_abbr) %>% unique() %>% nrow()
print(paste(bin,'(number):',tmp))
}



# Calculate mean value
bin_gather <- gather(bin_abundance,key='sample',value='abundance',AT4:ARP2)
## Add sample type column
bin_gather$sample_type <- bin_gather$sample
bin_gather$sample_type <- gsub("1|2|3|4|5","",bin_gather$sample_type)
## mean & sd
mean_sd_bin <- bin_gather %>% group_by(MAG,sample_type) %>% 
  mutate(mean = mean(abundance)) %>% 
  mutate(sd = sd(abundance)) %>% 
  select(MAG, sample_type, mean, sd) %>% 
  unique()
## Convert into spread
bin_spread <- mean_sd_bin %>% select(!(sd)) %>% 
  spread(mean_sd_bin, key='sample_type',value='mean')
bin_spread <- data.frame(bin_spread$MAG,c(unlist(bin_spread$ARP)),c(unlist(bin_spread$AT)))
#write.table(bin_spread, "../../airborne_arg_uwtp_result/metawrap_bin/bin_refinement/mag_mean.txt",
#          ,sep=' ',quote = FALSE, row.names = FALSE)