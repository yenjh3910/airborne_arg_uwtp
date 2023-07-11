# humann3_heatmap_without_stress.R

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
## P value (FDR adjust by "holm")
### Calculate p.value among each sample
res_AT_ARP <- gather_go %>% filter(!(sample_type == "ODP")) %>% group_by(Gene) %>% 
  do(w = wilcox.test(z_score~sample_type, data=., p.adjust.method = "holm",
                     paired=FALSE, exact = FALSE)) %>% 
  summarise(Gene, AT_ARP_Wilcox = w$p.value)
res_AT_ODP <- gather_go %>% filter(!(sample_type == "ARP")) %>% group_by(Gene) %>% 
  do(w = wilcox.test(z_score~sample_type, data=., p.adjust.method = "holm",
                     paired=FALSE, exact = FALSE)) %>% 
  summarise(Gene, AT_ODP_Wilcox = w$p.value)
res_ARP_ODP <- gather_go %>% filter(!(sample_type == "AT")) %>% group_by(Gene) %>% 
  do(w = wilcox.test(z_score~sample_type, data=., p.adjust.method = "holm",
                     paired=FALSE, exact = FALSE)) %>% 
  summarise(Gene, ARP_ODP_Wilcox = w$p.value)
### Join p.value dataframe
res <- full_join(res_AT_ARP, res_AT_ODP, "Gene")
res <- full_join(res, res_ARP_ODP, "Gene")
### Filter out the p.value that is bigger than 0.05 in all sample type
res <- res %>% filter((AT_ARP_Wilcox < 0.05)|(AT_ODP_Wilcox < 0.05)|(ARP_ODP_Wilcox < 0.05))
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
###### First look #######
GO_db$name[["GO:0043934"]] # View GO name
get_term_property(ontology=GO_db , property="ancestors", 
                  term="GO:0043934", as_names=TRUE) # View GO ancestors
get_term_property(ontology=GO_db , property="children", 
                  term="GO:0043934", as_names=TRUE) # View GO children
#Select gene
subGene <- GO_db$id[grep(GO_db$ancestors, pattern="GO:0043934")]
# Filter select gene from GO database
tmp <- go_zscore %>% filter(grepl(pattern = subGene[1], x = Gene))
for (i in 2:length(subGene)){
  tmp2 <- go_zscore %>% filter(grepl(pattern = subGene[i], x = Gene))
  tmp <- rbind(tmp,tmp2)
}
tmp$Trait <- GO_db$name[["GO:0043934"]]
##########################
# See sporulation species
tmp3 <- humann3_go %>% filter(grepl(pattern = subGene[1], x = Gene))
for (i in 2:length(subGene)){
  tmp4 <- humann3_go %>% filter(grepl(pattern = subGene[i], x = Gene))
  tmp3 <- rbind(tmp3,tmp4)
}
spore_species <- tmp3
spore_species <- spore_species %>% separate(Gene, c("Gene","Taxa"), sep = "\\|")
spore_species <- spore_species %>% separate(Taxa, c("Genus","Species"), sep = ".s__")
spore_species$Genus <- gsub("g__", "", spore_species$Genus)
##########################
###### Second look #######
GO_db$name[["GO:0043937"]] # View GO name
get_term_property(ontology=GO_db , property="ancestors", 
                  term="GO:0043937", as_names=TRUE) # View GO ancestors
get_term_property(ontology=GO_db , property="children", 
                  term="GO:0043937", as_names=TRUE) # View GO children
#Select gene
subGene <- GO_db$id[grep(GO_db$ancestors, pattern="GO:0043937")]
# Fiter select gene from GO database
for (i in 1:length(subGene)){
  tmp2 <- go_zscore %>% filter(grepl(pattern = subGene[i], x = Gene))
  tmp2$Trait <- GO_db$name[["GO:0043937"]]
  tmp <- rbind(tmp,tmp2)
}
###### Third look #######
GO_db$name[["GO:0009292"]] # View GO name
get_term_property(ontology=GO_db , property="ancestors", 
                  term="GO:0009292", as_names=TRUE) # View GO ancestors
get_term_property(ontology=GO_db , property="children", 
                  term="GO:0009292", as_names=TRUE) # View GO children
#Select gene
subGene <- GO_db$id[grep(GO_db$ancestors, pattern="GO:0009292")]
# Fiter select gene from GO database
for (i in 1:length(subGene)){
  tmp2 <- go_zscore %>% filter(grepl(pattern = subGene[i], x = Gene))
  tmp2$Trait <- GO_db$name[["GO:0009292"]]
  tmp <- rbind(tmp,tmp2)
}
###### Fourth look #######
GO_db$name[["GO:0009416"]] # View GO name
get_term_property(ontology=GO_db , property="ancestors", 
                  term="GO:0009416", as_names=TRUE) # View GO ancestors
get_term_property(ontology=GO_db , property="children", 
                  term="GO:0009416", as_names=TRUE) # View GO children
# Select gene
subGene <- GO_db$id[grep(GO_db$ancestors, pattern="GO:0009416")]
# Filter select gene from GO database
for (i in 1:length(subGene)){
  tmp2 <- go_zscore %>% filter(grepl(pattern = subGene[i], x = Gene))
  tmp2$Trait <- GO_db$name[["GO:0009416"]]
  tmp <- rbind(tmp,tmp2)
}

# Create trait dataframe
gene_trait <- tmp %>% select(Trait) %>% as.data.frame()
row.names(gene_trait) <- tmp$Gene
gene_trait$Trait <- str_to_sentence(gene_trait$Trait) # Covert first letter to uppercase
gene_trait$Trait[1:11] <- 'Sporulation / Regulation of sporulation'
# Order trait
trait_order <- c("Horizontal gene transfer",
                 "Sporulation / Regulation of sporulation",
                 "Response to light stimulus")
gene_trait <- gene_trait  %>%
  mutate(Trait =  factor(Trait, levels = trait_order)) %>%
  arrange(Trait) 
# gene_trait$Trait <- gsub('Regulation of ',
#                          '',
#                          gene_trait$Trait)
# Order tmp
tmp <- tmp %>%
  mutate(Gene =  factor(Gene, levels = row.names(gene_trait))) %>%
  arrange(Gene) 
select_go <- tmp
select_go <- select_go[,-17]

# Remove non bacteria host
# Unclassified, Fungi, Plant:
# GO:0030437, GO:0043935, GO:0048315, GO:0009704,
# GO:0009416, GO:0009585, GO:0009637, GO:0009638,
# GO:0009639, GO:0009640, GO:0009642, GO:0010205,
# GO:0010224, GO:0048573, GO:0048574
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0030437", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0043935", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0048315", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0009704", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0009416", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0009585", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0009637", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0009638", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0009639", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0009640", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0009642", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0010205", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0010224", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0048573", x = Gene)))
select_go <- select_go %>% filter(!(grepl(pattern = "^GO:0048574", x = Gene)))

# Transform to matrix 
select_go_mat <- as.matrix(select_go[,-1])
row.names(select_go_mat) <- select_go$Gene
# Reorder column
select_go_mat <- select_go_mat[, c(6,7,8,9,10,1,2,3,4,5,11,12,13,14,15)]
# Edit row.name 
row.names(select_go_mat) <- gsub('.\\[BP\\]','',row.names(select_go_mat))
row.names(gene_trait) <- gsub('.\\[BP\\]','',row.names(gene_trait))
# Annotation row
annotation_col = data.frame(Sample = factor(rep(c("AT", "ARP", "ODP"), 
                                                c(5, 5, 5))))
rownames(annotation_col) = colnames(select_go_mat)
# Color
display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")
ann_colors = list(Sample = c(AT = "#FB8072", ARP = "#80B1D3", ODP = "#B3DE69"),
                  Trait = c('Horizontal gene transfer' = "#BEBADA",
                            'Sporulation / Regulation of sporulation' = "#FFFFB3",
                            'Response to light stimulus' = "#FCCDE5"))

row.names(select_go_mat)

# Plot
p <-pheatmap(select_go_mat, cluster_cols = FALSE, cluster_rows = FALSE,
             annotation_row = gene_trait, annotation_col = annotation_col,
             clustering_distance_rows = "euclidean", annotation_colors = ann_colors,
             fontsize = 10, fontsize_row = 9, fontsize_col = 9,
             cellwidth = 10, cellheight = 10, bg = "transparent")

print(p)

# ggsave("GO_heatmap.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/humann3",
#        width = 11, height = 4,
#        units = "in", bg='transparent') # save to png format


# statistic
gat_select_go <- select_go %>% gather(key="sample", value = "z", ARP1:ODP5)
gat_select_go$sample_type <- gat_select_go$sample
gat_select_go$sample_type <- gsub("1|2|3|4|5","",gat_select_go$sample)
z_mean_sd <- gat_select_go %>% group_by(sample_type,Gene) %>% 
                               mutate(mean = mean(z)) %>% 
                               mutate(sd = sd(z)) %>% 
                               select(Gene,sample_type,mean,sd) %>% 
                               unique()