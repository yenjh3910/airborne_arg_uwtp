# humann3_heatmap.R

# Import
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
humann3_go <- read.table("../../airborne_arg_uwtp_result/humann3/GO/genefamilies_rename_go_cpm.tsv",
                         sep = "\t", header = TRUE, quote = "", comment.char = "")
# Rename column
colnames(humann3_go) <- c("Gene",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "AT1","AT2","AT3","AT4","AT5",
                          "ODP1","ODP2","ODP3","ODP4","ODP5")
# Filter out taxonomy row
go <- humann3_go %>% filter(!(grepl(pattern = "\\|", x = Gene)))
# Transform to gather format
gather_go <- gather(go, key = "sample", value = "cpm", 
                          ARP1:ODP5)
# Calculate z-score 
gather_go <- gather_go %>% group_by(Gene) %>% 
  mutate(z_score = (cpm-mean(cpm))/sd(cpm))
# View mean z-score
gather_go$sample_type <- gather_go$sample
gather_go$sample_type <- gsub("1|2|3|4|5","",gather_go$sample)
gather_go <- gather_go %>% group_by(sample_type) %>% mutate(mean_z = mean(z_score)) 
gather_go <- gather_go %>% arrange(desc(mean_z))
# Remove unnecessary column
go_zscore <- gather_go %>% ungroup() %>% select(Gene, sample, z_score)
# Transform to spread format
go_zscore <- spread(go_zscore, key = "sample", value = "z_score")
# Reorder gene by dominant gene in ARP
go_zscore <- go_zscore %>% group_by(Gene) %>% mutate(ARP_mean_z = (ARP1+ARP2+ARP3+ARP4+ARP5)/5)
go_zscore <- go_zscore %>% arrange(desc(ARP_mean_z))
go_zscore <- go_zscore %>% ungroup() %>% select(!(ARP_mean_z))

select_gene <- "GO:0034599|GO:0000302|
               |GO:0019430|GO:0051409|
               |GO:0047484|GO:0071470|
               |GO:0009414|GO:0009269|
               |GO:0009411|GO:0080183|
               |GO:0030435|GO:0043937|
               |GO:0030436|GO:0009847|
               |GO:0042173|GO:0042174|
               |GO:0043934|GO:0045881|
               |GO:0042244|GO:0070590|
               |GO:0000746|GO:0009291|
               |GO:0046999|GO:0009294|
               |GO:0030420|GO:0006303|
               |GO:0000724|GO:0009405"

tmp <- go_zscore %>% filter(grepl(pattern = select_gene, x = Gene))


# tmp <- go %>% filter((grepl(pattern = "GO:0009651|GO:1901002", x = Gene)))

tmp_mat <- as.matrix(tmp[,-1])
row.names(tmp_mat) <- tmp$Gene

pheatmap(tmp_mat, cluster_cols = FALSE)















humann3_kegg <- read.table("../../airborne_arg_uwtp_result/humann3/KEGG/genefamilies_rename_kegg-orthology_cpm.tsv",
                    sep = "\t", header = TRUE, quote = "", comment.char = "")
colnames(humann3_kegg) <- c("Gene",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "AT1","AT2","AT3","AT4","AT5",
                          "ODP1","ODP2","ODP3","ODP4","ODP5")
kegg <- humann3_kegg %>% filter(!(grepl(pattern = "\\|", x = Gene)))



humann3_pfam <- read.table("../../airborne_arg_uwtp_result/humann3/pfam/genefamilies_rename_pfam_cpm.tsv",
                           sep = "\t", header = TRUE, quote = "", comment.char = "")
colnames(humann3_pfam) <- c("Gene",
                            "ARP1","ARP2","ARP3","ARP4","ARP5",
                            "AT1","AT2","AT3","AT4","AT5",
                            "ODP1","ODP2","ODP3","ODP4","ODP5")
pfam <- humann3_pfam %>% filter(!(grepl(pattern = "\\|", x = Gene)))





path <- read.table("../../airborne_arg_uwtp_result/humann3/path_abundance/pathabundance_cpm.tsv",
                           sep = "\t", header = TRUE, quote = "", comment.char = "")
colnames(path) <- c("Pathway",
                    "ARP1","ARP2","ARP3","ARP4","ARP5",
                    "AT1","AT2","AT3","AT4","AT5",
                    "ODP1","ODP2","ODP3","ODP4","ODP5")
path <- path %>% filter(!(grepl(pattern = "\\|", x = Pathway)))
