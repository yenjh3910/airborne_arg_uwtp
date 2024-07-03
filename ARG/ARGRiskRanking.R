# ARGRiskRanking.R

# Import library
library(tidyverse)
library(openxlsx)
# Read ARG_type file
gene_cell <- read.table("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.gene.txt",
                       header = TRUE, sep = "")
# gather data
gather_gene <- gather(gene_cell, key = "sample", value = "copy_per_cell", 
                             ARP1:ODP5) # Transform to gather format
## Add sample type column
gather_gene$sample_type <- NA
for (i in 1:nrow(gather_gene)) {
  if (grepl("^AT",gather_gene$sample[i])) {
    gather_gene$sample_type[i] <- "AT"
  } 
  if (grepl("^ARP",gather_gene$sample[i])) {
    gather_gene$sample_type[i] <- "ARP"
  } 
  if (grepl("^ODP",gather_gene$sample[i])) {
    gather_gene$sample_type[i] <- "ODP"
  }
}
# Read Rank file
rank_db <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/ARG_rank_db.xlsx",
                         sheet = 1)
rank_db <- rank_db %>% select(variant,rank,Type,Subtype)
colnames(rank_db)[1] <- "gene"
# Merge two file
gene_cell <- left_join(gather_gene,rank_db, by='gene')
# Calculate gene sum
sum_gene_cell <- gene_cell %>% group_by(sample,rank)  %>% summarize(copy_per_cell=sum(copy_per_cell))
## Add sample type column
sum_gene_cell$sample_type <- NA
for (i in 1:nrow(sum_gene_cell)) {
  if (grepl("^AT",sum_gene_cell$sample[i])) {
    sum_gene_cell$sample_type[i] <- "AT"
  } 
  if (grepl("^ARP",sum_gene_cell$sample[i])) {
    sum_gene_cell$sample_type[i] <- "ARP"
  } 
  if (grepl("^ODP",sum_gene_cell$sample[i])) {
    sum_gene_cell$sample_type[i] <- "ODP"
  }
}
# Calculate gene mean
mean_gene_cell <- sum_gene_cell %>% group_by(sample_type,rank)  %>% summarize(copy_per_cell=mean(copy_per_cell)) # Level of ARG in each sample
mean_gene_cell <- mean_gene_cell %>% group_by(sample_type) %>% mutate(sum = sum(copy_per_cell))
mean_gene_cell <- mean_gene_cell %>% mutate(percentage = copy_per_cell/sum*100)
# Calculate gene mean and sd
sum_gene_cell <- sum_gene_cell %>% group_by(sample_type,rank)  %>% mutate(mean=mean(copy_per_cell))
sum_gene_cell <- sum_gene_cell %>% group_by(sample_type,rank) %>% mutate(sd = sd(copy_per_cell))
# Filter necesssary column for figure
mean_gene_cell <- sum_gene_cell %>% select(rank,sample_type,mean,sd) %>% unique()
# Add st_mean for figure
mean_gene_cell <- mean_gene_cell %>% group_by(sample_type) %>% 
  mutate(sd_mean = mean)
mean_gene_cell <- mean_gene_cell %>% 
  arrange(factor(rank, levels = c('I','II','III','IV','notassessed'))) %>% 
  arrange(sample_type)
## sum ARP st_mean
for (i in c(1:4)) {
  mean_gene_cell[{i},5] <-  mean_gene_cell[{i},5] + 
    sum(mean_gene_cell[{i+1}:5,5])}
## sum AT st_mean
for (i in c(6:9)) {
  mean_gene_cell[{i},5] <-  mean_gene_cell[{i},5] + 
    sum(mean_gene_cell[{i+1}:10,5])}
## sum ODP st_mean
for (i in c(11:14)) {
  mean_gene_cell[{i},5] <- mean_gene_cell[{i},5] + 
    sum(mean_gene_cell[{i+1}:15,5])}

# Plot
# Order
mean_gene_cell$sample_type <- factor(mean_gene_cell$sample_type, 
                                    levels = c("ODP",
                                               "ARP",
                                               "AT"))
mean_gene_cell$rank[mean_gene_cell$rank == 'notassessed'] <- 'Not assessed'
rank_col <- c( "#FB8072", "#80B1D3", "#FDB462","#B3DE69","#BEBADA",
               "#FB8072", "#80B1D3", "#FDB462","#B3DE69","#BEBADA",
               "#FB8072", "#80B1D3", "#FDB462","#B3DE69","#BEBADA")
p<-ggplot(mean_gene_cell, aes(x = sample_type, y = mean, fill = rank)) + 
  geom_bar(stat="identity")+
  geom_errorbar( aes(x=sample_type, ymin=sd_mean-sd, ymax=sd_mean+sd), 
                 width=0.3, color=rank_col, alpha=0.9, size=0.6)+
  coord_flip()+
  theme_bw()+
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5])))+
  xlab("")+ ylab("Relative abundance (ARGs/cell)")+
  scale_fill_manual(values=c("#FB8072", "#80B1D3", "#FDB462",
                             "#B3DE69","#BEBADA"))+
  guides(fill=guide_legend(title="Risk Rank",byrow=TRUE))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        # plot.margin=unit(c(5,1,1,1),"cm"),
        legend.position = 'top',
        legend.spacing.x = unit(0.5, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        legend.key.size = unit(0.8, 'cm'),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) #transparent legend bg)
  
print(p)
# ggsave("risk_rank_abundance.png", p, path = "../../airborne_arg_uwtp_result/Figure/ARG",
#         width = 11, height = 4, units = "in") # save to png format

# Risk 1: 1.314313
mean_gene_cell$copy_per_cell[1]/mean_gene_cell$copy_per_cell[6]
# Risk 2: 1.14296
mean_gene_cell$copy_per_cell[2]/mean_gene_cell$copy_per_cell[7]
# Risk 3: 1.307035
mean_gene_cell$copy_per_cell[3]/mean_gene_cell$copy_per_cell[8]
# Risk 4: 1.753797
mean_gene_cell$copy_per_cell[4]/mean_gene_cell$copy_per_cell[9]
