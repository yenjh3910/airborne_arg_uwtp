# pathogen_coverage.R

library(stringr)
source("./ARG_coverage.R")
# Import pathogen list
pathogen_list <- read.table("./pathogen_list.txt",quote="",sep = ",",col.names = "Pathogen")
# Remove unclassified species
final_ARG_coverage <- final_ARG_coverage %>% filter(!(Species=="Unclassified"))
# Edit ARG type
final_ARG_coverage$type <- str_to_title(final_ARG_coverage$type)
final_ARG_coverage$type[final_ARG_coverage$type == "Macrolide-Lincosamide-Streptogramin"] <- "MLS"
final_ARG_coverage$type[final_ARG_coverage$type == "Beta_lactam"] <- "Beta-lactam"
final_ARG_coverage$type <- gsub("_"," ",final_ARG_coverage$type)
final_ARG_coverage$type[final_ARG_coverage$type == "Tetracenomycin c"] <- "Tetracenomycin C"
# Select ARG carrying pathogen
tmp <- data_frame()
for (i in 1:length(pathogen_list[,1])){
  tmp2 <- final_ARG_coverage %>% filter(grepl(pathogen_list[i,1], Species))
  tmp <- rbind(tmp,tmp2)
}
pathogen_arg_coverage <- tmp

# # Calculate coverage mean and sd
# pathogen_arg_coverage %>% group_by(Sample_type) %>% 
#                           mutate(mean = sum(coverage)/5) %>% 
#                           mutate(sd = sqrt(sum((coverage-mean)^2/(5-1)))) %>% 
#                           select(Sample_type,mean,sd) %>% 
#                           unique()

# # Calculate mean coverage by ARG type
# CoverageArgType <- pathogen_arg_coverage %>% group_by(type,Sample_type) %>% 
#                               summarise(coverage_sum = sum(coverage)/5)

# Calculate sum coverage by ARG type
CoverageArgType <- pathogen_arg_coverage %>% group_by(type,Sample_type) %>%
                              summarise(coverage_sum = sum(coverage))

# Calculate coverage sum
pathogen_arg_coverage %>% group_by(Sample_type) %>%
                          summarise(sum = sum(coverage))

# Order sample_type
CoverageArgType$Sample_type <- factor(CoverageArgType$Sample_type, 
                                 levels = c("ODP","ARP","AT"))
# Order ARG type
CoverageArgType <- CoverageArgType %>% arrange(desc(coverage_sum))
unique(CoverageArgType$type)
CoverageArgType$type <- factor(CoverageArgType$type,
                            levels = c("Beta-lactam"  , "Rifamycin"             , "MLS"            ,        
                                       "Tetracycline" , "Multidrug"             , "Aminoglycoside" ,        
                                       "Sulfonamide"  , "Polymyxin"             , "Chloramphenicol",       
                                       "Bacitracin"   , "Pleuromutilin tiamulin", "Quinolone"   ))
# Plot barplot
p <- ggplot(CoverageArgType, aes(x = Sample_type, y = coverage_sum, fill = type)) + 
  geom_bar(stat="identity",width = 0.7) + 
  theme_bw() + 
  xlab("") + ylab("Coverage (x/GB)") + 
  scale_fill_brewer(palette="Set3") +
  coord_flip() + 
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]),
                            "ODP" = expression(Outdoor~PM[2.5]))) +
  guides(fill=guide_legend(ncol=2, bycol=TRUE,title = "ARG type")) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent')) #transparent legend bg)

print(p)

# ggsave("ARG_pathogen_coverage.png", p, path = "../../airborne_arg_uwtp_result/Figure/pathogen",
#         width = 7.7, height = 2, units = "in") # save to png format
 
# Bubble chart of ARG coverage in pathogen
pathogen_arg_coverage$Species %>% unique()  # Check Species
## Change Species name to brief version
for (i in 1:length(pathogen_arg_coverage[,1])){
  pathogen_arg_coverage$Species[i] <- sub(".*",
                                          str_extract(pathogen_arg_coverage$Species[i],
                                                      "[^ ]* [^ ]*"),
                                          pathogen_arg_coverage$Species[i])}
pathogen_arg_coverage$Species %>% unique() # Check Species

# ## Calculate mean coverage by pathogen
# CoveragePathogen <- pathogen_arg_coverage %>% group_by(Species,Sample_type) %>% 
#                           summarise(coverage_sum = sum(coverage)/5)

## Calculate sum coverage by pathogen
CoveragePathogen <- pathogen_arg_coverage %>% group_by(Species,Sample_type) %>% 
  summarise(coverage_sum = sum(coverage))

## Filter pathogen only occur in  both AT and ARP
only_AT <- CoveragePathogen %>% filter(Sample_type == "AT")
only_ARP <- CoveragePathogen %>% filter(Sample_type == "ARP")
AerCoveragePathogen <- CoveragePathogen %>% filter(Species %in% only_AT$Species) %>% 
                                            filter(Species %in% only_ARP$Species)
# Order pathogen
AerCoveragePathogen %>% filter(Sample_type=="ARP") %>% 
                        arrange(desc(coverage_sum)) %>% 
                        select(Species) %>% unique()
AerCoveragePathogen$Species <- factor(AerCoveragePathogen$Species,
                               levels = c("Gordonia amarae"             ,
                                          "Pseudomonas aeruginosa"      ,
                                          "Stenotrophomonas maltophilia",
                                          "Burkholderia cepacia"        ,
                                          "Aeromonas veronii"           ,
                                          "Klebsiella pneumoniae"       ,
                                          "Enterococcus faecalis"       ,
                                          "Morganella morganii"         ,
                                          "Streptococcus suis"          ))
## Plot
p <- ggplot(AerCoveragePathogen, aes(x = Sample_type, y = Species)) + 
  geom_point(aes(size = coverage_sum,fill = Sample_type), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(0.00001, 100), range = c(1,22), breaks = c(1,5,15,30)) +
  theme_bw() + 
  coord_flip() +
  theme_bw() + 
  labs(x="",y="",size = "ARGs coverage (x/GB)", fill = "") +
  scale_x_discrete(labels=c("AT" = expression(Aeration~tank), 
                            "ARP" = expression(Aeration~tank~PM[2.5]))) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"), guide = FALSE) +
  theme(axis.text.x = element_text(size = 10, angle = 70, hjust=1),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.key=element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent'))
print(p)
# ggsave("ARG_pathogen_bubble.png", p, path = "../../airborne_arg_uwtp_result/Figure/pathogen",
#         width = 6.78, height = 3.62, units = "in") # save to png format
