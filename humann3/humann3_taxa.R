# humann3_taxa.R

# Import
library(tidyverse)
library(ontologyIndex)
library(pheatmap)
library(RColorBrewer)
library(FSA)
humann3_go <- read.table("../../airborne_arg_uwtp_result/humann3/GO/genefamilies_rename_go_cpm.tsv",
                         sep = "\t", header = TRUE, quote = "", comment.char = "")
# Rename column
colnames(humann3_go) <- c("Gene",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "AT1","AT2","AT3","AT4","AT5",
                          "ODP1","ODP2","ODP3","ODP4","ODP5")
# Separate taxa
humann3_go <- humann3_go %>% separate(col = 'Gene' ,sep = '\\|',c("Gene","Taxa"))
humann3_go <- humann3_go %>% drop_na()
humann3_go <- humann3_go %>% separate(col = 'Taxa' ,sep = '.s__',c("Genus","Species"))
humann3_go$Genus <- gsub("g__", "", humann3_go$Genus)
humann3_go$Species[is.na(humann3_go$Species)] <- 'unclassified'

# Transform to gather format
gather_go <- gather(humann3_go, key = "sample", value = "cpm", 
                    ARP1:ODP5)
# Add sample type
gather_go$sample_type <- gather_go$sample
gather_go$sample_type <- gsub("1|2|3|4|5","",gather_go$sample)

osmotic <- gather_go %>% filter(grepl(pattern = "^GO:0071470", x = Gene))
osmotic <- osmotic %>% group_by(Genus,sample_type) %>% summarise(mean_cpm = mean(cpm))
osmotic <- osmotic %>% filter(!(Genus == 'unclassified'))
osmotic <- osmotic %>% filter(!(sample_type == 'ODP'))

salt <- gather_go %>% filter(grepl(pattern = "^GO:0009651", x = Gene))
salt <- salt %>% group_by(Genus,sample_type) %>% summarise(mean_cpm = mean(cpm))
salt <- salt %>% filter(!(Genus == 'unclassified'))
salt <- salt %>% filter(!(sample_type == 'ODP'))

water <- gather_go %>% filter(grepl(pattern = "^GO:0009414", x = Gene))
water <- water %>% group_by(Genus,sample_type) %>% summarise(mean_cpm = mean(cpm))
water <- water %>% filter(!(Genus == 'unclassified'))
water <- water %>% filter(!(sample_type == 'ODP'))

starvation <- gather_go %>% filter(grepl(pattern = "^GO:0009267", x = Gene))
starvation <- starvation %>% group_by(Genus,sample_type) %>% summarise(mean_cpm = mean(cpm))
starvation <- starvation %>% filter(!(Genus == 'unclassified'))
starvation <- starvation %>% filter(!(sample_type == 'ODP'))

hyperosmotic <- gather_go %>% filter(grepl(pattern = "^GO:0006072", x = Gene))
hyperosmotic <- hyperosmotic %>% group_by(Genus,sample_type) %>% summarise(mean_cpm = mean(cpm))
hyperosmotic <- hyperosmotic %>% filter(!(Genus == 'unclassified'))
hyperosmotic <- hyperosmotic %>% filter(!(sample_type == 'ODP'))

DNA_break <- gather_go %>% filter(grepl(pattern = "^GO:0006281", x = Gene))
DNA_break <- DNA_break %>% group_by(Genus,sample_type) %>% summarise(mean_cpm = mean(cpm))
DNA_break <- DNA_break %>% filter(!(Genus == 'unclassified'))
DNA_break <- DNA_break %>% filter(!(sample_type == 'ODP'))

# ## set the levels in order we want
# tmp <- tmp %>% ungroup() %>% group_by(Genus) %>% mutate(total = sum(mean_cpm))
# tmp$Genus <- factor(tmp$Genus,                                    # Factor levels in decreasing order
#                     levels = tmp$Genus[order(tmp$total, decreasing = TRUE)])

# Transform to spread format
tmp <- spread(tmp, key = "sample_type", value = "mean_cpm")
# Remove AT = 0
tmp <- tmp %>% filter(!(AT == 0))
tmp <- tmp %>% group_by(Genus) %>% mutate(ratio = ARP/AT)

ggplot(tmp, aes(x=Genus, y=mean_cpm, fill=sample_type)) + 
  geom_bar(stat="identity", position=position_dodge(), alpha=0.6) +
  theme_bw()