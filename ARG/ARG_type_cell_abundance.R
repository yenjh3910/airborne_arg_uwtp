# ARG_type_cell_abundance.R
# Calculate and visualize ARG type/cell in each sample

# Import library
library(tidyverse)
library(stringr)
library(tibble)
library(FSA)

# Read ARG_type file
arg_type <- read.table("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.type.txt",
                   header = TRUE, sep = "")
# Preview raw data
gather_arg_type <- gather(arg_type, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5) # Transform to gather format
## Add sample type column
gather_arg_type$sample_type <- NA
for (i in 1:nrow(gather_arg_type)) {
  if (grepl("^AT",gather_arg_type$sample[i])) {
    gather_arg_type$sample_type[i] <- "AT"
  } 
  if (grepl("^ARP",gather_arg_type$sample[i])) {
      gather_arg_type$sample_type[i] <- "ARP"
  } 
  if (grepl("^ODP",gather_arg_type$sample[i])) {
      gather_arg_type$sample_type[i] <- "ODP"
  }
}
gather_arg_type %>% group_by(sample) %>% summarize(sum(copy_per_cell)) # Level of ARG in each sample
type_percentage <- gather_arg_type # Use for calculation of type percentage & type mean,sd finally
# Statastic of sample type
arg_sum <- gather_arg_type %>% group_by(sample) %>% mutate(sum = sum(copy_per_cell))
arg_sum <- arg_sum %>% group_by(sample_type) %>% mutate(mean = mean(sum))
arg_sum <- arg_sum %>% group_by(sample_type)%>% mutate(sd = sd(sum)) %>% 
  select(sample_type, mean, sd) %>% unique()
## Visualization
gather_arg_type$sample_type <- factor(gather_arg_type$sample_type, levels = c("AT", "ARP", "ODP"), 
                  labels = c(expression(Aeration~tank), 
                             expression(Aeration~tank~PM[2.5]), 
                             expression(Outdoor~PM[2.5])
                             )) # Change facet title
ggplot(gather_arg_type, aes(x = sample, y = copy_per_cell, fill = type)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~sample_type, scales = "free_x", labeller = label_parsed) + 
  theme_bw() + 
  xlab("") + ylab("Relative abundance (ARGs/cell)")


# Classify minimum ARG type to others
average_arg <- gather_arg_type %>%
  group_by(sample_type, type) %>%
  summarise_at(vars(copy_per_cell), funs(mean)) # Calculate mean value of ARG type in sample type

arg_order <- unique(
  average_arg[
    order(average_arg$copy_per_cell,
                    decreasing = T),
    ]$type) # unique ARG type order

arg_type <- arg_type %>% 
  arrange(factor(type, levels = arg_order)) # Order ARG type in spread format

## Convert the first column into row name
rownames(arg_type) <- arg_type[,1]
arg_type[,1] <- NULL

## Calculate others by summing  minimum arg
other_arg <- colSums(arg_type[11:26, ])
arg_type <- rbind(arg_type, other_arg)
rownames(arg_type)[rownames(arg_type) == "27"] <- "others" 
arg_type <- arg_type[-(11:26),] # delete minimum arg type
rownames(arg_type) <- str_to_title(rownames(arg_type)) # Covert first letter to uppercase
arg_type <- tibble::rownames_to_column(arg_type, "type") # Convert row name back to first column


# Plot
## Transform to gather format
gather_arg_type <- gather(arg_type, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5)
## Add sample type column
gather_arg_type$sample_type <- NA
for (i in 1:nrow(gather_arg_type)) {
  if (grepl("^AT",gather_arg_type$sample[i])) {
    gather_arg_type$sample_type[i] <- "AT"
  } 
  if (grepl("^ARP",gather_arg_type$sample[i])) {
    gather_arg_type$sample_type[i] <- "ARP"
  } 
  if (grepl("^ODP",gather_arg_type$sample[i])) {
    gather_arg_type$sample_type[i] <- "ODP"
  }
}
## Change facet title
gather_arg_type$sample_type <- factor(gather_arg_type$sample_type, levels = c("AT", "ARP", "ODP"), 
                                      labels = c(expression(Aeration~tank), 
                                                 expression(Aeration~tank~PM[2.5]), 
                                                 expression(Outdoor~PM[2.5])
                                      ))
## Change specific ARG type
gather_arg_type$type[gather_arg_type$type == "Macrolide-Lincosamide-Streptogramin"] <- "MLS"
gather_arg_type$type[gather_arg_type$type == "Beta_lactam"] <- "Beta-lactam"
## Order ARG type
gather_arg_type$type <- factor(gather_arg_type$type, 
                               levels = c("Multidrug", "Sulfonamide",
                               "MLS","Aminoglycoside", 
                               "Tetracycline", "Rifamycin", 
                               "Beta-lactam", "Bacitracin",
                               "Polymyxin", "Chloramphenicol", 
                               "Others"))
## ggplot
p <- ggplot(gather_arg_type, aes(x = sample, y = copy_per_cell, fill = type)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~sample_type, scales = "free_x", labeller = label_parsed) + 
  theme_bw() + 
  xlab("") + ylab("Relative abundance (ARGs/cell)") + 
  guides(fill=guide_legend(title="ARG type")) + 
  # scale_y_continuous(expand = c(0, 0)) + # y start at 0
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust=1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size=10.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        # panel.grid.major = element_blank(), #remove major gridlines
        # panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent')) + #transparent legend bg
        scale_fill_manual(values=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3",
                                   "#FDB462","#B3DE69","#FCCDE5", "#FFED6F","#BC80BD",
                                   "#D9D9D9"))

print(p)

# ggsave("ARG_type_cell.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 7, height = 5,
#        units = "in", bg='transparent') # save to png format

# Statistic
gather_arg_type <- gather_arg_type %>% group_by(sample) %>% 
                                       mutate(sum = sum(copy_per_cell))
gather_arg_type <- gather_arg_type %>% select(sample,sample_type,sum) %>% unique()
# Kruskal-Wallis Test
kruskal.test(sum ~ sample_type, data = gather_arg_type)
# Post hoc of Kruskal-Wallis Test (DunnTest)
dunnTest(sum ~ sample_type, data=gather_arg_type, method="holm")
# # Mann Whitney U Test (Wilcoxon Rank Sum Test)
pairwise.wilcox.test(gather_arg_type$sum, gather_arg_type$sample_type,
                     p.adjust.method = "holm")
pairwise.wilcox.test(gather_arg_type$sum, gather_arg_type$sample_type,
                     p.adjust.method = "bonferroni")
pairwise.wilcox.test(gather_arg_type$sum, gather_arg_type$sample_type,
                     p.adjust.method = "BH")

# Calculate mean abundance between sample type
## Covert first letter to uppercase
type_percentage$type <- str_to_title(type_percentage$type)
## Change specific ARG type
type_percentage$type[type_percentage$type == "Macrolide-Lincosamide-Streptogramin"] <- "MLS"
type_percentage$type[type_percentage$type == "Beta_lactam"] <- "Beta-lactam"
type_percentage$type <- gsub("_"," ",type_percentage$type)
type_percentage$type[type_percentage$type == "Tetracenomycin c"] <- "Tetracenomycin C"
# Calculate mean and sd of ARG type in each sample type
final_type_mean_sd <- type_percentage %>% group_by(sample_type,type) %>%
                      mutate(mean = mean(copy_per_cell)) %>% 
                      mutate(sd = sd(copy_per_cell)) %>% 
                      select(type,sample_type,mean,sd) %>% 
                      unique()
final_type_mean_sd %>% filter(sample_type == "AT") %>% arrange(desc(mean)) # Print mean & sd of AT
final_type_mean_sd %>% filter(sample_type == "ARP") %>% arrange(desc(mean)) # Print mean & sd of ARP
final_type_mean_sd %>% filter(sample_type == "ODP") %>% arrange(desc(mean)) # Print mean & sd of ODP
# Calculation type percentage
final_pecent <- type_percentage %>% group_by(sample_type,type) %>% 
                    mutate(abundance = sum(copy_per_cell)) %>% 
                    select(!(sample)) %>% select(!(copy_per_cell)) %>% 
                    unique() %>% 
                    ungroup() %>% 
                    group_by(sample_type) %>% 
                    mutate(sample_type_sum = sum(abundance)) %>% 
                    mutate(percentage = (abundance/sample_type_sum)*100)
final_pecent_spread <- final_pecent %>% select(!(abundance)) %>% 
                                        select(!(sample_type_sum)) %>% 
                                        spread(sample_type,percentage)
final_pecent_spread <- final_pecent_spread[,c(1,3,2,4)]
# Export to csv
# write.csv(final_pecent_spread,
#           "../../airborne_arg_uwtp_result/args_oap/ARG/ARG_type_percentage.csv",
#           row.names = FALSE)



# ARG type barplot
type_percentage %>% ggplot(aes(x=sample_type, y=copy_per_cell))+
  geom_boxplot(aes(fill=sample_type)) + 
  facet_wrap(~ type, ncol=7,scales='free') + 
  theme_bw()+ 
  labs(x="") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())