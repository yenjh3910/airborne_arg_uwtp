# bind_ARG_MGE_VF_coverage.R

library(ggpubr)
library(ggplot2)
library(gggenes)
# Run coverage calculation script
source("./ARG_coverage.R")
source("./MGE_coverage.R")
source("./VF_coverage.R")

### Three import data:
### final_ARG_coverage, MGE_coverage, VF_coverage

# Change ARG type to uppercase
final_ARG_coverage$type <- str_to_title(final_ARG_coverage$type) 
final_ARG_coverage$type[final_ARG_coverage$type == "Macrolide-Lincosamide-Streptogramin"] <- "MLS"
final_ARG_coverage$type[final_ARG_coverage$type == "Beta_lactam"] <- "Beta-lactam"
## Change specific MGE type
MGE_coverage$type[MGE_coverage$type == "insertion_element_IS91"] <- "IS91"
MGE_coverage$type[MGE_coverage$type == "transposase"] <- "Transposase"
MGE_coverage$type[MGE_coverage$type == "plasmid"] <- "Plasmid"
MGE_coverage$type[MGE_coverage$type == "integrase"] <- "Integrase"
# Create new column to index
final_ARG_coverage$Contigs_Sample <- paste(final_ARG_coverage$contigs,
                                           final_ARG_coverage$Sample,
                                           sep = "_")
MGE_coverage$Contigs_Sample <- paste(MGE_coverage$contigs,
                                     MGE_coverage$SampleID.x,
                                     sep = "_")
VF_coverage$Contigs_Sample <- paste(VF_coverage$contigs,
                                    VF_coverage$SampleID.x,
                                    sep = "_")
# Fiter by contigs only in ARG
ARG_contigs <- final_ARG_coverage$Contigs_Sample
MGE_coverage_ARG <- MGE_coverage %>% filter(Contigs_Sample %in% ARG_contigs)
VF_coverage_ARG <- VF_coverage %>% filter(Contigs_Sample %in% ARG_contigs)
# Create tax corresponding to contigs
tax_contigs <- final_ARG_coverage[,7:15]
tax_contigs <- tax_contigs[,-8]
tax_contigs <- unique(tax_contigs)
# Add tax tag
MGE_coverage_ARG <- left_join(MGE_coverage_ARG,tax_contigs,by="Contigs_Sample")
VF_coverage_ARG <- left_join(VF_coverage_ARG,tax_contigs,by="Contigs_Sample")
# Edit column before binding
final_ARG_coverage$subtype <- str_replace(final_ARG_coverage$subtype, ".+__", "")
final_ARG_coverage$Gene_type <- "ARG"
MGE_coverage_ARG <- MGE_coverage_ARG %>% select(!("Avg_fold")) %>% select(!("Length") ) %>%
                                         select(!("Std_Dev")) %>% select(!('taxonomy ID'))
colnames(MGE_coverage_ARG)[2] <- "Sample"
MGE_coverage_ARG$Sample_type <- gsub("1|2|3|4|5","",MGE_coverage_ARG$Sample)
MGE_coverage_ARG$Gene_type <- "MGE"
VF_coverage_ARG <- VF_coverage_ARG %>% select(!("Avg_fold")) %>% select(!("Length") ) %>%
                                       select(!("Std_Dev")) %>% select(!('taxonomy ID'))
colnames(VF_coverage_ARG)[2] <- "Sample"
VF_coverage_ARG$Sample_type <- gsub("1|2|3|4|5","",VF_coverage_ARG$Sample)
VF_coverage_ARG$Gene_type <- "VF"
VF_coverage_ARG$gene <- str_replace(VF_coverage_ARG$gene, "\\[.+\\]", "")
VF_coverage_ARG$gene <- str_replace(VF_coverage_ARG$gene, ".$", "")
VF_coverage_ARG <- VF_coverage_ARG %>%
                   separate(gene, sep="^\\S*\\K\\s+",into=c("subtype","type"))
VF_coverage_ARG$subtype <- str_replace(VF_coverage_ARG$subtype, "\\)", "")
VF_coverage_ARG$subtype <- str_replace(VF_coverage_ARG$subtype, "\\(", "")
# Bind dataframe
final_coverage <- rbind(final_ARG_coverage,MGE_coverage_ARG)
final_coverage <- rbind(final_coverage,VF_coverage_ARG)


# Correlation between ARG & MGE
## Prepare df
mean_ARG <- final_ARG_coverage %>% select("Sample_type","coverage","Contigs_Sample") %>% 
                              group_by(Contigs_Sample) %>% mutate(ARG=mean(coverage)) %>% 
                              select(!("coverage")) %>% unique()
mean_MGE <- MGE_coverage_ARG %>% select("coverage","Contigs_Sample")%>% 
                             group_by(Contigs_Sample) %>% mutate(MGE=mean(coverage)) %>% 
                             select(!("coverage")) %>% unique()
ARG_MGE_cooccurence <- left_join(mean_ARG,mean_MGE,by="Contigs_Sample")
ARG_MGE_cooccurence <- na.omit(ARG_MGE_cooccurence)
ARG_MGE_cooccurence <- ARG_MGE_cooccurence %>%  filter(!(ARG == 0)) %>% filter(!(MGE == 0))
ARG_MGE_cooccurence$Sample_type <- factor(ARG_MGE_cooccurence$Sample_type, levels = c("AT","ARP","ODP"))
## Prepare df of each sample type
noODP_cor <-ARG_MGE_cooccurence %>% filter(!(Sample_type == "ODP"))
ODP_cor <-ARG_MGE_cooccurence %>% filter((Sample_type == "ODP"))
AT_cor <-ARG_MGE_cooccurence %>% filter(Sample_type == "AT")
ARP_cor <-ARG_MGE_cooccurence %>% filter(Sample_type == "ARP")
## Transfer to log value
ARP_cor$ARG <- log10(ARP_cor$ARG)
ARP_cor$MGE <- log10(ARP_cor$MGE)
AT_cor$ARG <- log10(AT_cor$ARG)
AT_cor$MGE <- log10(AT_cor$MGE)
ODP_cor$ARG <- log10(ODP_cor$ARG)
ODP_cor$MGE <- log10(ODP_cor$MGE)
noODP_cor$ARG <- log10(noODP_cor$ARG)
noODP_cor$MGE <- log10(noODP_cor$MGE)
# Regression
AT_model <- lm(MGE ~ ARG, data = AT_cor)
ARP_model <- lm(MGE ~ ARG, data = ARP_cor)
summary(AT_model)
summary(ARP_model)
# Plot 
###### Fine tune figure in powerpoint (add regression formula R^2 and move legend position)
## AT & ARP together
p <- ggscatter(noODP_cor, x = "ARG", y = "MGE",
          color = 'Sample_type', shape = 'Sample_type', 
          add = "reg.line", conf.int = TRUE, 
          size = 2.5, alpha = 0.6, ellipse.alpha = 0.6) + 
  # stat_cor(aes(color= Sample_type), show.legend = FALSE) +
  labs(x = expression(ARGs~coverage~(log[10]*"(×/Gb)")), 
       y = expression(MGEs~coverage~(log[10]*"(×/Gb)"))) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size=10.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        aspect.ratio=1, 
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) +  #transparent legend bg
 scale_color_manual(name = "", 
                     labels = c(expression(Aeration~tank),
                                expression(Aeration~tank~PM[2.5])),
                     values = c("AT" = "#F8766D",
                                "ARP" = "#00BFC4")) + 
  scale_fill_manual(name = "", 
                    labels = c(expression(Aeration~tank),
                               expression(Aeration~tank~PM[2.5])),
                    values = c("AT" = "#F8766D",
                               "ARP" = "#00BFC4"))
print(p)
# ggsave("Contigs_ARG_MGE_correlation.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/correlation",
#        width = 6, height = 5)

## AT only
p <- ggscatter(AT_cor, x = "ARG", y = "MGE",
          color = 'Sample_type', shape = 'Sample_type', 
          add = "reg.line", conf.int = TRUE, 
          size = 2, alpha = 0.6, ellipse.alpha = 0.6) + 
  #stat_cor(aes(color= Sample_type), show.legend = FALSE) +
  labs(x = expression(ARGs~coverage~(log[10]*"(×/Gb)")), 
       y = expression(MGEs~coverage~(log[10]*"(×/Gb)"))) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size=10.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        aspect.ratio=1, 
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) +  #transparent legend bg
  scale_color_manual(name = "", 
                     labels = c(expression(Aeration~tank)),
                     values = c("AT" = "#F8766D")) + 
  scale_fill_manual(name = "", 
                    labels = c(expression(Aeration~tank)),
                    values = c("AT" = "#F8766D"))
print(p)
# ggsave("AT_Contigs_ARG_MGE_correlation.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/correlation",
#        width = 6, height = 5)

## ARP only
p <- ggscatter(ARP_cor, x = "ARG", y = "MGE",
          color = 'Sample_type', shape = 'Sample_type', 
          add = "reg.line", conf.int = TRUE, 
          size = 2, alpha = 0.6, ellipse.alpha = 0.6) + 
  #stat_cor(aes(color= Sample_type), show.legend = FALSE) +
  labs(x = expression(ARGs~coverage~(log[10]*"(×/Gb)")), 
       y = expression(MGEs~coverage~(log[10]*"(×/Gb)"))) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size=10.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        aspect.ratio=1, 
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) +  #transparent legend bg
  scale_color_manual(name = "", 
                     labels = c(expression(Aeration~tank~PM[2.5])),
                     values = c("ARP" = "#00BFC4")) + 
  scale_fill_manual(name = "", 
                    labels = c(expression(Aeration~tank~PM[2.5])),
                    values = c("ARP" = "#00BFC4"))
# ggsave("ARP_Contigs_ARG_MGE_correlation.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/correlation",
#        width = 6, height = 5)


# Intersection between ARG,MGE,& VF
## Get ARG occurance with MGE & VF
MGE_contigs <- MGE_coverage$Contigs_Sample
VF_contigs <- VF_coverage$Contigs_Sample
ARG_cooccur_MGE <- final_coverage %>% filter(Contigs_Sample %in% MGE_contigs)
ARG_cooccur_VF <- final_coverage %>% filter(Contigs_Sample %in% VF_contigs)
cooccur_ARG_MGE_VF <- final_coverage %>% filter(Contigs_Sample %in% VF_contigs) %>% 
                                         filter(Contigs_Sample %in% MGE_contigs)
######## Intersection between ARG,MGE,& VF = 0 !!!! ############


# gggenes
## Join with ORF position
ORF_position <- read.csv("../../airborne_arg_uwtp_result/contigs_ORF_position/contigs_ORF_position.csv")
ARG_cooccur_MGE$ORF_SampleID <- paste(ARG_cooccur_MGE$ORF,
                                      ARG_cooccur_MGE$Sample,
                                      sep = "_")
ARG_cooccur_MGE <- left_join(ARG_cooccur_MGE, ORF_position, by = "ORF_SampleID")
## Finetune format for figure
ARG_cooccur_MGE["direction"][ARG_cooccur_MGE["direction"] == -1] <- 0
colnames(ARG_cooccur_MGE)[5] <- "Type"
colnames(ARG_cooccur_MGE)[6] <- "Contigs"


### ARG & MGE intersection more than three in "ARP"
ARG_MGE_IntersectMore3 <- ARG_cooccur_MGE %>% filter(Sample_type == "ARP") %>% 
  group_by(Contigs_Sample) %>%
  summarise(count = n()) %>% 
  as.data.frame() %>% 
  filter(count >= 3)
ARG_MGE_cooccurARP <- ARG_cooccur_MGE %>% 
                      filter(Contigs_Sample %in% ARG_MGE_IntersectMore3$Contigs_Sample)
### Select gene manually
select_contigs_sample <- c("k141_243726_ARP5","k141_480594_ARP2",
                           "k141_489310_ARP4","k141_50739_ARP4")
plot_contigs <- ARG_MGE_cooccurARP %>% filter(Contigs_Sample %in% select_contigs_sample)
#plot_contigs <- plot_contigs %>% filter(!(Type == "qacEdelta"))
### Plot
ggplot(plot_contigs, aes(xmin = start, xmax = end, y = Contigs, fill = Type,
                         forward = direction, label = subtype)) +
  geom_gene_arrow() +
  geom_gene_label(align = "left") +
  facet_wrap(~ Contigs, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

### ARG & MGE intersection more than three in "AT"
ARG_MGE_IntersectMore3 <- ARG_cooccur_MGE %>% filter(Sample_type == "AT") %>% 
  group_by(Contigs_Sample) %>%
  summarise(count = n()) %>% 
  as.data.frame() %>% 
  filter(count >= 3)
ARG_MGE_cooccurAT <- ARG_cooccur_MGE %>% 
  filter(Contigs_Sample %in% ARG_MGE_IntersectMore3$Contigs_Sample)
### Select gene manually
select_contigs_sample <- c("k141_313559_AT2","k141_34906_AT4",
                           "k141_454684_AT1","k141_55118_AT5",
                           "k141_712198_AT1","k141_77012_AT3")
plot_contigs <- ARG_MGE_cooccurAT %>% filter(Contigs_Sample %in% select_contigs_sample)
#plot_contigs <- plot_contigs %>% filter(!(Type == "qacEdelta"))
### Plot
ggplot(plot_contigs, aes(xmin = start, xmax = end, y = Contigs, fill = Type,
                         forward = direction, label = subtype)) +
  geom_gene_arrow() +
  geom_gene_label(align = "left") +
  facet_wrap(~ Contigs, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()


#### Plot gggene in aeration tank area (AT + ARP)
ARG_MGE_IntersectMore3 <- ARG_cooccur_MGE %>% filter(!(Sample_type == "ODP")) %>% 
  group_by(Contigs_Sample) %>%
  summarise(count = n()) %>% 
  as.data.frame() %>% 
  filter(count >= 3)
ARG_MGE_cooccurATARP <- ARG_cooccur_MGE %>% 
                        filter(Contigs_Sample %in% ARG_MGE_IntersectMore3$Contigs_Sample)
### Select gene manually
select_contigs_sample <- c("k141_313559_AT2",
                           "k141_712198_AT1",
                           "k141_243726_ARP5")
plot_contigs <- ARG_MGE_cooccurATARP %>% filter(Contigs_Sample %in% select_contigs_sample)
plot_contigs <- plot_contigs[,-11] # Deal with qacEdelta
plot_contigs["subtype"][plot_contigs["subtype"] == "qacEdelta"] <- "qacEdelta1" # Deal with qacEdelta
plot_contigs["Type"][plot_contigs["Type"] == "qacEdelta"] <- "Multidrug" # Deal with qacEdelta
## Order gene type
plot_contigs$Type <- factor(plot_contigs$Type, 
                            levels = c("Aminoglycoside","Florfenicol",
                                       "Multidrug","Sulfonamide",
                                       "Tetracycline",
                                       "Integrase","IS91",
                                       "Transposase"))
### Plot
#### View with gene label, but we need to add label afterward manually
ggplot(plot_contigs, aes(xmin = start, xmax = end, y = Contigs, fill = Type,
                         forward = direction, label = subtype)) +
  geom_gene_arrow() +
  geom_gene_label(align = "left") +
  facet_wrap(~ Contigs, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()
#### Output without label
p <- ggplot(plot_contigs, aes(xmin = start, xmax = end, y = Contigs, fill = Type,
                         forward = direction)) +
  geom_gene_arrow() +
  facet_wrap(~ Contigs, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

print(p)
# ggsave("ARG_MGE_cooccurance.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/gggene",
#        width = 10, height = 5,
#        units = "in", bg='transparent') # save to png format
