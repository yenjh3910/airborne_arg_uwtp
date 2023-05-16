# taxa_alpha_dicersity.R

library(vegan)
library(tidyverse)
library(ggbreak)
library(openxlsx)
library(FSA)

# Merge bracken file
## Import ARP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ARP/ARP1.S.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ARP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ARP/ARP",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ARP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
ARP_taxa<- tmp
## Import AT
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/AT/AT1.S.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "AT1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/AT/AT",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("AT", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
AT_taxa<- tmp
## Import ODP
tmp <- read.table("../../airborne_arg_uwtp_result/kraken2/ODP/ODP1.S.bracken", 
                  header = TRUE ,sep = "\t")
tmp <- tmp %>% select(name, new_est_reads)
colnames(tmp)[2] <- "ODP1"
for (i in 2:5) {
  tmp2 <- read.table(paste("../../airborne_arg_uwtp_result/kraken2/ODP/ODP",i,".S.bracken", sep = ""), 
                     header = TRUE ,sep = "\t") # Import
  tmp2 <- tmp2 %>% select(name, new_est_reads) # Select name and reads column
  colnames(tmp2)[2] <- paste("ODP", i, sep = "") # Change column headers to sample name
  tmp <- full_join(tmp ,tmp2) # Merge each sample file
  tmp <- tmp %>% replace(is.na(.), 0) # Replace NA with 0
}
colnames(tmp)[1] <- "taxa"
ODP_taxa<- tmp
## Full join all sample type file
ARP_AT_taxa <- full_join(ARP_taxa, AT_taxa)
merge_taxa <- full_join(ARP_AT_taxa, ODP_taxa)
merge_taxa <- merge_taxa %>% replace(is.na(.), 0) # Replace NA with 0
# Filtering low abundance read
THRESHOLD <- 0.00005 # 0.005 % as filtering threshold
for (i in 2:16) {
  individual_threshold <- sum(merge_taxa[,i])*THRESHOLD
  merge_taxa[,i] <- replace(merge_taxa[,i], merge_taxa[,i] < individual_threshold, 0)
}
# Fitst column as row name
rownames(merge_taxa) <- merge_taxa$taxa
merge_taxa <-  merge_taxa[,-1]
# Remove all 0 row
merge_taxa <- merge_taxa[rowSums(merge_taxa[])>0,]
# Transform row name back to column one
merge_taxa$taxa <- rownames(merge_taxa)
# Gather transformation
gat_merge_taxa <- merge_taxa %>% gather(key = "sample", value = "abundance",ARP1:ODP5)
# Add sample type column
gat_merge_taxa$sample_type <- gat_merge_taxa$sample
gat_merge_taxa$sample_type <- gsub("1|2|3|4|5","",gat_merge_taxa$sample)
# Calculate aplha diversity
alpha_diversity <- gat_merge_taxa %>%
  group_by(sample) %>%
  summarize(sobs = specnumber(abundance),
            shannon = diversity(abundance, index="shannon"),
            simpson = diversity(abundance, index="simpson"),
            invsimpson = 1/simpson,
            n = sum(abundance))
# Add sample type column
alpha_diversity$sample_type <- alpha_diversity$sample
alpha_diversity$sample_type <- gsub("1|2|3|4|5","",alpha_diversity$sample_type)
# Calculate mean and sd
mean_alpha <- alpha_diversity %>%  
                    group_by(sample_type) %>% 
                    mutate(sobs_mean = mean(sobs)) %>%
                    mutate(sobs_sd = sd(sobs)) %>%
                    mutate(shannon_mean = mean(shannon)) %>% 
                    mutate(shannon_sd = sd(shannon)) %>% 
                    mutate(simpson_mean = mean(simpson)) %>% 
                    mutate(simpson_sd = sd(simpson)) %>% 
                    mutate(invsimpson_mean = mean(invsimpson)) %>% 
                    mutate(invsimpson_sd = sd(invsimpson)) %>%
                    mutate(n_mean = mean(n)) %>% 
                    mutate(n_sd = sd(n)) %>%
                    select(!(sample)) %>% 
                    select(!(sobs)) %>% 
                    select(!(shannon)) %>% 
                    select(!(simpson)) %>% 
                    select(!(invsimpson)) %>%
                    select(!(n)) %>%
                    unique()
# Order
mean_alpha$sample_type <- factor(mean_alpha$sample_type, 
                                 levels = c("AT",
                                            "ARP",
                                            "ODP"))

# Plot Richness
p <- ggplot(mean_alpha, aes(x=sample_type, y=sobs_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=sobs_mean-sobs_sd, ymax=sobs_mean+sobs_sd), width=.2) +
  theme_bw() +
  xlab("")+ ylab("Observed Richness") + 
  scale_fill_manual(values=c("#F8766D", "#00BFC4", "#7CAE00")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) + 
  #ylim(0, 4000) +
  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 15),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        # hideaxis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank()) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 4000))

print(p)
#Save
# ggsave("richness.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/taxonomy",
#        width = 2.735, height = 5.815, units = "in") # save to png format

######################################################################

# Plot Shannon
p <- ggplot(mean_alpha, aes(x=sample_type, y=shannon_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=shannon_mean-shannon_sd, ymax=shannon_mean+shannon_sd), width=.2) +
  theme_bw() +
  xlab("")+ ylab("Alpha Diversity Measure") + 
  scale_fill_manual(values=c("#F8766D", "#00BFC4", "#7CAE00")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) +
  scale_y_break(c(1, 5.0), scales = 10) + 
  ylim(0, 8) +
  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 15),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        # hideaxis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())

print(p)
# # Save
# ggsave("shannon_index.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/taxonomy",
#        width = 2.909, height = 6, units = "in") # save to png format

######################################################################

# Plot Simpson
p <- ggplot(mean_alpha, aes(x=sample_type, y=invsimpson_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=invsimpson_mean-invsimpson_sd, ymax=invsimpson_mean+invsimpson_sd), width=.2) +
  theme_bw() +
  xlab("")+ ylab("Alpha Diversity Measure") + 
  scale_fill_manual(values=c("#F8766D", "#00BFC4", "#7CAE00")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) + 
  scale_y_break(c(0.1,0.98), scales = 10) + 
  ylim(0, 1.04) +
  theme(axis.text.x = element_text(angle = 45, hjust=1,size = 15),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 18),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        # hideaxis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank())
print(p)
# #Save
# ggsave("simpson_index.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/taxonomy",
#        width = 3, height = 6, units = "in") # save to png format

#################################################################################

# Statistic
## Richness
kruskal.test(sobs ~ sample_type, data = alpha_diversity)
dunnTest(sobs ~ sample_type, data=alpha_diversity, method="holm")
pairwise.wilcox.test(alpha_diversity$sobs, alpha_diversity$sample_type,
                     p.adjust.method = "holm")

## Shannon
kruskal.test(shannon ~ sample_type, data = alpha_diversity)
dunnTest(shannon ~ sample_type, data=alpha_diversity, method="holm")
pairwise.wilcox.test(alpha_diversity$shannon, alpha_diversity$sample_type,
                     p.adjust.method = "holm")
## Invsimpson
kruskal.test(invsimpson ~ sample_type, data = alpha_diversity)
dunnTest(invsimpson ~ sample_type, data=alpha_diversity, method="holm")
pairwise.wilcox.test(alpha_diversity$invsimpson, alpha_diversity$sample_type,
                     p.adjust.method = "holm")
