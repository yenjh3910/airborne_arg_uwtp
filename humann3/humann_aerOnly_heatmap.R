# humann3_aerOnly_heatmap.R

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
# Remove ODP sample
humann3_go <- humann3_go[,-(12:16)]
# Filter out taxonomy row
go <- humann3_go %>% filter(!(grepl(pattern = "\\|", x = Gene)))
# Transform to gather format
gather_go <- gather(go, key = "sample", value = "cpm", 
                    ARP1:AT5)
# Calculate z-score 
gather_go <- gather_go %>% group_by(Gene) %>% 
  mutate(z_score = (cpm-mean(cpm))/sd(cpm))
# Caluclate p=value & view mean z-score
gather_go$sample_type <- gather_go$sample
gather_go$sample_type <- gsub("1|2|3|4|5","",gather_go$sample)

## P value (FDR adjust by "holm")
### Calculate p.value between AT & ARP
res_AT_ARP <- gather_go %>% group_by(Gene) %>% 
  do(w = wilcox.test(cpm~sample_type, data=., p.adjust.method = "holm",
                     paired=FALSE, exact = FALSE)) %>%
  summarise(Gene, AT_ARP_Wilcox = w$p.value)
### Join p.value dataframe
res <- res_AT_ARP
### Filter out the p.value that is bigger than 0.05 in all sample type
res <- res %>% filter((AT_ARP_Wilcox < 0.05))
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
GO_db$name[["GO:0006950"]] # View GO name
get_term_property(ontology=GO_db , property="ancestors", 
                  term="GO:0006950", as_names=TRUE) # View GO ancestors
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
stress_go_mat <- stress_go_mat[, c(6,7,8,9,10,1,2,3,4,5)]
# Edit row.name 
row.names(stress_go_mat) <- gsub('.\\[BP\\]','',row.names(stress_go_mat))
# Annotation row
annotation_row = data.frame(Sample = factor(rep(c("AT", "ARP"), 
                                                c(5, 5))))
rownames(annotation_row) = colnames(stress_go_mat)
# Color
display.brewer.all()
brewer.pal(7, "Set3")
ann_colors = list(Sample = c(AT = "#FB8072", ARP = "#80B1D3"))
# Plot
p <-pheatmap(stress_go_mat, cluster_cols = TRUE, 
             clustering_distance_rows = "euclidean",
             annotation_col = annotation_row, annotation_colors = ann_colors,
             fontsize = 10, fontsize_row = 7, fontsize_col = 8,
             cellwidth = 10, cellheight = 8, bg = "transparent")

# ggsave("GO_aer_stress.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/humann3",
#        width = 7, height = 7,
#        units = "in", bg='transparent') # save to png format


# Statistic
gat_stress_go <- stress_go %>% gather(key="sample", value = "z", ARP1:AT5)
gat_stress_go$sample_type <- gat_stress_go$sample
gat_stress_go$sample_type <- gsub("1|2|3|4|5","",gat_stress_go$sample)
z_mean_sd <- gat_stress_go %>% group_by(sample_type,Gene) %>%
  mutate(mean = mean(z)) %>%
  mutate(sd = sd(z)) %>%
  select(Gene,sample_type,mean,sd) %>%
  unique()