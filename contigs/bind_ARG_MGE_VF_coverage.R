# bind_ARG_MGE_VF_coverage.R

# Run coverage calculation script
source("./ARG_coverage.R")
source("./MGE_coverage.R")
source("./VF_coverage.R")

### Three import data:
### final_ARG_coverage, MGE_coverage, VF_coverage

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



# # Edit column name
# colnames(final_ARG_coverage)[3] <- "ARG_coverage"
# colnames(final_ARG_coverage)[4] <- "ARG_subtype"
# colnames(final_ARG_coverage)[5] <- "ARG_type"
# colnames(MGE_coverage)[6] <- "MGE_coverage"
# colnames(MGE_coverage)[7] <- "MGE_subtype"
# colnames(MGE_coverage)[8] <- "MGE_type"
# colnames(VF_coverage)[6] <- "VF_coverage"
# colnames(VF_coverage)[7] <- "VF_gene"
# # Select necessary column
# MGE_coverage <- MGE_coverage %>% select("MGE_coverage","MGE_subtype","MGE_type","Contigs_Sample")
# VF_coverage <- VF_coverage %>% select("VF_coverage","VF_gene","Contigs_Sample")
# # Join dataframe
# test <- left_join(final_ARG_coverage,MGE_coverage,by="Contigs_Sample")
