# humann3_heatmap.R

# Import
library(tidyverse)
library(ontologyIndex)
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
# Caluclate p=value & view mean z-score
gather_go$sample_type <- gather_go$sample
gather_go$sample_type <- gsub("1|2|3|4|5","",gather_go$sample)

# # Filter out the 0 value >= 2
# gather_go %>% group_by(Gene, sample_type) %>% filter(cpm == 0)
# t <- gather_go %>% group_by(Gene, sample_type) %>% filter(cpm == 0) %>% count("cpm")
# t %>% count("Gene")

## P value
### Calculate p.value among each sample
res_AT_ARP <- gather_go %>% filter(!(sample_type == "ODP")) %>% group_by(Gene) %>% 
  do(w = wilcox.test(z_score~sample_type, data=., paired=FALSE, exact = FALSE)) %>% 
  summarise(Gene, AT_ARP_Wilcox = w$p.value)
res_AT_ODP <- gather_go %>% filter(!(sample_type == "ARP")) %>% group_by(Gene) %>% 
  do(w = wilcox.test(z_score~sample_type, data=., paired=FALSE, exact = FALSE)) %>% 
  summarise(Gene, AT_ODP_Wilcox = w$p.value)
res_ARP_ODP <- gather_go %>% filter(!(sample_type == "AT")) %>% group_by(Gene) %>% 
  do(w = wilcox.test(z_score~sample_type, data=., paired=FALSE, exact = FALSE)) %>% 
  summarise(Gene, ARP_ODP_Wilcox = w$p.value)
### Join p.value dataframe
res <- full_join(res_AT_ARP, res_AT_ODP, "Gene")
res <- full_join(res, res_ARP_ODP, "Gene")
### Filter out the p.value that is bigger than 0.05 in all sample type
res <- res %>% filter(!((AT_ARP_Wilcox > 0.05)&(AT_ODP_Wilcox > 0.05)&(ARP_ODP_Wilcox > 0.05)))
## mean z-score
gather_go <- gather_go %>% group_by(sample_type) %>% mutate(mean_z = mean(z_score)) 
gather_go <- gather_go %>% arrange(desc(mean_z))
# Remove unnecessary column
go_zscore <- gather_go %>% ungroup() %>% select(Gene, sample, z_score)
# Transform to spread format
go_zscore <- spread(go_zscore, key = "sample", value = "z_score")
# Filter p value by res
go_zscore <- go_zscore %>% filter(Gene %in% res$Gene)
# Reorder gene by dominant gene in ARP
go_zscore <- go_zscore %>% group_by(Gene) %>% mutate(ARP_mean_z = (ARP1+ARP2+ARP3+ARP4+ARP5)/5)
go_zscore <- go_zscore %>% arrange(desc(ARP_mean_z))
go_zscore <- go_zscore %>% ungroup() %>% select(!(ARP_mean_z))
# Import GO database & select gene
GO_db <- get_OBO("../../airborne_arg_uwtp_result/humann3/GO/go-basic.obo")
GO_db$name[["GO:0000747"]] # View GO name
get_term_property(ontology=GO_db , property="ancestors", 
                  term="GO:0000747", as_names=TRUE) # View GO ancestors
get_term_property(ontology=GO_db , property="children", 
                  term="GO:0006950", as_names=TRUE) # View GO children
##### response to stress #####
stress_subGene <- GO_db$id[grep(GO_db$ancestors, pattern="GO:0006950")]
# Fiter select gene from GO database
tmp <- go_zscore %>% filter(grepl(pattern = stress_subGene[1], x = Gene))
for (i in 2:length(stress_subGene)){
  tmp2 <- go_zscore %>% filter(grepl(pattern = stress_subGene[i], x = Gene))
  tmp <- rbind(tmp,tmp2)
}
stress_go <- tmp
# Transform to matrix 
stress_go_mat <- as.matrix(stress_go[,-1])
row.names(stress_go_mat) <- stress_go$Gene
# Reorder column
stress_go_mat <- stress_go_mat[, c(6,7,8,9,10,1,2,3,4,5,11,12,13,14,15)]
# Edit row.name 
row.names(stress_go_mat) <- gsub('.\\[BP\\]','',row.names(stress_go_mat))
# Annotation row
annotation_row = data.frame(Sample = factor(rep(c("AT", "ARP", "ODP"), 
                                                c(5, 5, 5))))
rownames(annotation_row) = colnames(stress_go_mat)
# Color
display.brewer.all()
brewer.pal(7, "Set3")
ann_colors = list(Sample = c(AT = "#FB8072", ARP = "#80B1D3", ODP = "#B3DE69"))
# Plot
p <-pheatmap(stress_go_mat, cluster_cols = FALSE, 
         clustering_distance_rows = "euclidean",
         annotation_col = annotation_row, annotation_colors = ann_colors,
         fontsize = 10, fontsize_row = 6, fontsize_col = 7,
         cellwidth = 9, cellheight = 6, bg = "transparent")

# ggsave("GO_stress.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/humann3",
#        width = 7, height = 7,
#        units = "in", bg='transparent') # save to png format

pt <-pheatmap(t(stress_go_mat), cluster_rows = FALSE, 
             clustering_distance_rows = "euclidean",
             annotation_row = annotation_row, annotation_colors = ann_colors,
             fontsize = 10, fontsize_row = 6, fontsize_col = 7,
             cellwidth = 7, cellheight = 6, bg = "transparent")



################ Try
stress_subGene <- GO_db$id[grep(GO_db$ancestors, pattern="GO:0009292")] # Horizontal gene transfer
stress_subGene <- GO_db$id[grep(GO_db$ancestors, pattern="GO:0043934")] # Sporulation
stress_subGene <- GO_db$id[grep(GO_db$ancestors, pattern="GO:0006970")]
# Fiter select gene from GO database
tmp <- go_zscore %>% filter(grepl(pattern = stress_subGene[1], x = Gene))
for (i in 2:length(stress_subGene)){
  tmp2 <- go_zscore %>% filter(grepl(pattern = stress_subGene[i], x = Gene))
  tmp <- rbind(tmp,tmp2)
}
stress_go <- tmp
# Transform to matrix 
stress_go_mat <- as.matrix(stress_go[,-1])
row.names(stress_go_mat) <- stress_go$Gene
# Reorder column
stress_go_mat <- stress_go_mat[, c(6,7,8,9,10,1,2,3,4,5,11,12,13,14,15)]
# Annotation row
annotation_row = data.frame(Sample = factor(rep(c("AT", "ARP", "ODP"), 
                                                c(5, 5, 5))))
rownames(annotation_row) = colnames(stress_go_mat)
# Plot
pheatmap(stress_go_mat, cluster_cols = FALSE, 
         clustering_distance_rows = "euclidean",
         annotation_col = annotation_row)










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

tmp_gene <- "GO:0047484|GO:0006970|GO:0006974|GO:0009607|GO:0071260|GO:0009612|
|GO:0016048|GO:0071476|GO:0071474|GO:0042594|GO:0009267|GO:0006855"
tmp <- go_zscore %>% filter(grepl(pattern = tmp_gene, x = Gene))

starvation_gene <- "starvation"
tmp <- go_zscore %>% filter(grepl(pattern = starvation_gene, x = Gene))

stress_gene <- "stress"
tmp <- go_zscore %>% filter(grepl(pattern = stress_gene, x = Gene))

death_gene <- "death"
tmp <- go_zscore %>% filter(grepl(pattern = death_gene, x = Gene))

HGT_gene <- "transformation|conjugation"
tmp <- go_zscore %>% filter(grepl(pattern = HGT_gene, x = Gene))

efflux_gene <- "efflux"
tmp <- go_zscore %>% filter(grepl(pattern = efflux_gene, x = Gene))

membrane_gene <- "membrane"
tmp <- go_zscore %>% filter(grepl(pattern = membrane_gene, x = Gene))

transmembrane_gene <- "transmembrane"
tmp <- go_zscore %>% filter(grepl(pattern = transmembrane_gene, x = Gene))

secretion_gene <- "secretion"
tmp <- go_zscore %>% filter(grepl(pattern = secretion_gene, x = Gene))











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
