# ARG_alpha_dicersity.R

# Import library
library(tidyverse)
library(openxlsx)

# Read ARG_subtype file (Something wrong with read.table, so read by )
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
# Gather transformation
gat_arg <- arg_subtype %>% gather(key = "sample", value = "abundance",ARP1:ODP5)
# Add sample type column
gat_arg$sample_type <- gat_arg$sample
gat_arg$sample_type <- gsub("1|2|3|4|5","",gat_arg$sample_type)
# Calculate alpha diversity
alpha_diversity <- gat_arg %>%
  group_by(sample) %>%
  summarize(sobs = specnumber(abundance),
            shannon = diversity(abundance, index="shannon"),
            simpson = diversity(abundance, index="simpson"),
            invsimpson = diversity(abundance, index="invsimpson"),
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
                                 levels = c("ODP",
                                            "ARP",
                                            "AT"))
# Plot Richness
p <- ggplot(mean_alpha, aes(x=sample_type, y=sobs_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=sobs_mean-sobs_sd, ymax=sobs_mean+sobs_sd), width=.2) +
  theme_bw() +
  xlab("")+ ylab("Observed Richness") + 
  scale_fill_manual(values=c("#7CAE00","#00BFC4","#F8766D")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust=1,size = 15),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        # hideaxis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.line.y.left = element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 700)) +
  coord_flip() +
  ggtitle("Observed Richness")

print(p)
# #Save
# ggsave("arg_richness.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 5.5, height = 2.9, units = "in") # save to png format

# Plot Shannon
p <- ggplot(mean_alpha, aes(x=sample_type, y=shannon_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=shannon_mean-shannon_sd, ymax=shannon_mean+shannon_sd), width=.2) +
  theme_bw() +
  xlab("")+
  ggtitle("Shannon diversity") +
  scale_fill_manual(values=c("#7CAE00","#00BFC4","#F8766D")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) +
  scale_y_break(c(1, 4.0), scales = 10) + 
  ylim(0, 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust=1,size = 15),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        # hideaxis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.line.y.left = element_blank()) +
  coord_flip()

print(p)
# #Save
# ggsave("arg_shannon.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 5, height = 2.9, units = "in") # save to png format

# Plot Simpson
p <- ggplot(mean_alpha, aes(x=sample_type, y=invsimpson_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=invsimpson_mean-invsimpson_sd, ymax=invsimpson_mean+invsimpson_sd), width=.2) +
  theme_bw() +
  xlab("")+
  ggtitle("Inverse Simpson diversity") +
  scale_fill_manual(values=c("#7CAE00","#00BFC4","#F8766D")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) + 
  #scale_y_break(c(0.1,1), scales = 10) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust=1,size = 15),
        axis.text.y = element_blank(),
        axis.title =element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        # hideaxis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.line.y.left = element_blank())+
  coord_flip()
print(p)
# # Save
# ggsave("arg_simpson.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 5, height = 2.6, units = "in") # save to png format

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

############################################################################3
# Evenness
library(chemodiv)
## Prepare input
t_arg_subtype <- t(arg_subtype)
colnames(t_arg_subtype) <- t_arg_subtype[1,]
t_arg_subtype <- t_arg_subtype[-1,]
t_arg_subtype <- as.data.frame(t_arg_subtype)
subtype_eve_input <- as.data.frame(lapply(t_arg_subtype, as.numeric))
row.names(subtype_eve_input) <- row.names(t_arg_subtype)
## Calculate
eve <- calcDiv(subtype_eve_input, compDisMat = NULL, type = "PielouEven", q = 1)
## Post-process
### Add sample column
eve <- eve %>% mutate(sample = row.names(t_arg_subtype))
### Add sample type column
eve$sample_type <- eve$sample
eve$sample_type <- gsub("1|2|3|4|5","",eve$sample_type)
# Calculate mean and sd
mean_eve <- eve %>%  
  group_by(sample_type) %>% 
  mutate(even_mean = mean(PielouEven)) %>%
  mutate(even_sd = sd(PielouEven)) %>% 
  select(!(sample)) %>% 
  select(!(PielouEven)) %>% 
  unique()
# Statistic
kruskal.test(PielouEven ~ sample_type, data = eve)
dunnTest(PielouEven ~ sample_type, data=eve, method="holm")
pairwise.wilcox.test(eve$PielouEven, eve$sample_type,
                     p.adjust.method = "holm")
# Order
mean_eve$sample_type <- factor(mean_eve$sample_type, 
                                 levels = c("ODP",
                                            "ARP",
                                            "AT"))
### Plot
p <- ggplot(mean_eve, aes(x=sample_type, y=even_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=even_mean-even_sd, ymax=even_mean+even_sd), width=.2) +
  theme_bw() +
  xlab("")+ ggtitle('Pielou evenness') +
  scale_fill_manual(values=c("#7CAE00","#00BFC4","#F8766D")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) + 
  scale_y_break(c(0.1,0.7), scales = 10) + 
  ylim(0, 0.9) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust=1,size = 13),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        # hideaxis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.line.y.left = element_blank())+
  coord_flip()
print(p)
# Save
# ggsave("arg_evenness.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 4.9, height = 2.9, units = "in") # save to png format



# Legend
mean_eve$sample_type <- factor(mean_eve$sample_type, 
                               levels = c("AT",
                                          "ARP",
                                          "ODP"))
p <- ggplot(mean_eve, aes(x=sample_type, y=even_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=even_mean-even_sd, ymax=even_mean+even_sd), width=.2) +
  theme_bw() +
  xlab("") +
  scale_fill_manual(values=c("#F8766D","#00BFC4","#7CAE00"),
                    labels = c(expression(Aeration~tank),
                               expression(Aeration~tank~PM[2.5]),
                               expression(Outdoor~PM[2.5]))) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) + 
  scale_y_break(c(0.1,0.7), scales = 10) + 
  ylim(0, 0.9) +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        axis.text.x = element_text(angle = 45, hjust=1,size = 13),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        # hideaxis
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.text.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.line.y.left = element_blank())+
  coord_flip()
print(p)
# Save
# ggsave("arg_alpha_legend.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 4.9, height = 2.9, units = "in") # save to png format




# Aeration area only
aer_alpha <- mean_alpha %>% filter(sample_type!='ODP')
# Plot Richness
p <- ggplot(aer_alpha, aes(x=sample_type, y=sobs_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=sobs_mean-sobs_sd, ymax=sobs_mean+sobs_sd), width=.2) +
  theme_bw() +
  xlab("")+ ylab("Observed Richness") + 
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5])))+ 
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
        axis.line.y.right = element_blank())+
scale_y_continuous(expand = c(0, 0), limits = c(0, 700))

print(p)


# Plot Shannon
p <- ggplot(aer_alpha, aes(x=sample_type, y=shannon_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=shannon_mean-shannon_sd, ymax=shannon_mean+shannon_sd), width=.2) +
  theme_bw() +
  xlab("")+ ylab("Alpha Diversity Measure") + 
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]))) +
  scale_y_break(c(1, 4.0), scales = 10) + 
  #ylim(0, 8) +
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
        axis.line.y.right = element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5))

print(p)

# Plot Simpson
p <- ggplot(aer_alpha, aes(x=sample_type, y=invsimpson_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=invsimpson_mean-invsimpson_sd, ymax=invsimpson_mean+invsimpson_sd), width=.2) +
  theme_bw() +
  xlab("")+ ylab("Alpha Diversity Measure") + 
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]))) + 
  scale_y_break(c(0.1,1), scales = 10) + 
  #ylim(0, 1.04) +
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
        axis.line.y.right = element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05))

print(p)

# Even
aer_eve <- mean_eve %>% filter(sample_type!='ODP')
p <- ggplot(aer_eve, aes(x=sample_type, y=even_mean, fill=sample_type)) + 
  geom_bar(stat="identity",alpha=0.6) +
  geom_errorbar(aes(ymin=even_mean-even_sd, ymax=even_mean+even_sd), width=.2) +
  theme_bw() +
  xlab("")+ ylab("Evenness Measure") + 
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]))) + 
  scale_y_break(c(0.1,0.7), scales = 10) + 
  #ylim(0, 1.04) +
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

## Statistic
aer_alpha_statistic <- alpha_diversity %>% filter(sample_type!='ODP')
wilcox.test(sobs ~ sample_type, data = aer_alpha_statistic,
                   exact = FALSE)
wilcox.test(shannon ~ sample_type, data = aer_alpha_statistic,
            exact = FALSE)
wilcox.test(invsimpson ~ sample_type, data = aer_alpha_statistic,
            exact = FALSE)

aer_eve_statistic <- eve %>% filter(sample_type!="ODP")
wilcox.test(PielouEven ~ sample_type, data = aer_eve_statistic,
                   exact = FALSE)



div <- calcDiv(subtype_eve_input, compDisMat = NULL, type = c("Shannon","Simpson","PielouEven","HillEven"), q = 1)
## Post-process
### Add sample column
div <- div %>% mutate(sample = row.names(t_arg_subtype))
### Add sample type column
div$sample_type <- div$sample
div$sample_type <- gsub("1|2|3|4|5","",div$sample_type)
