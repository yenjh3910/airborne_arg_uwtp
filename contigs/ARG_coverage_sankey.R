# ARG_coverage_sankey.R

source("./ARG_coverage.R")
### final_ARG_coverage is the element we want ###
library(networkD3)
library(RColorBrewer)
# Aeration tank (AT)
## Calculate species mean coverage according to sum coverage in each sample
AT_species_coverage <- final_ARG_coverage %>% filter(Sample_type == "AT") %>% 
                       filter(!(Species == "Unclassified")) %>% 
                       group_by(Species,Sample) %>% 
                       summarise(sum_coverage = sum(coverage))
AT_species_coverage <-  AT_species_coverage %>% ungroup() %>% group_by(Species) %>% 
                        summarise(mean_coverage = mean(sum_coverage))
## join AT_species_coverage with coverage
taxa <- final_ARG_coverage %>% select("Domain","Phylum","Class",
                                      "Order","Family","Genus","Species") %>% 
                               unique()
AT_species_coverage <- left_join(AT_species_coverage,taxa,by="Species")
## Select top 10 species
AT_species_coverage <- AT_species_coverage %>% arrange(desc(mean_coverage))
AT_species_coverage <- AT_species_coverage[1:10,]
# Make Sankey node
taxa_list <- c("Domain","Phylum","Class",
               "Order","Family","Genus","Species")
tmp <- AT_species_coverage %>% select(taxa_list[1]) %>% unique()
colnames(tmp) <- "name"
for (i in 2:length(taxa_list)) {
  tmp2 <- AT_species_coverage %>% select(taxa_list[i]) %>% unique()
  colnames(tmp2) <- "name"
  tmp <- rbind(tmp,tmp2)
}
nodes <- tmp
# Make Sankey links
## Domain-Phylum
tmp <- AT_species_coverage %>% select(taxa_list[1],taxa_list[2]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- AT_species_coverage %>% group_by(Phylum) %>% mutate(value = sum(mean_coverage)) %>% 
                           select(Phylum,value) %>% unique()
colnames(tmp2)[1] <- "target"
links <- full_join(tmp,tmp2)
## Phylum-Class
tmp <- AT_species_coverage %>% select(taxa_list[2],taxa_list[3]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- AT_species_coverage %>% group_by(Class) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Class,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
## Class-Order
tmp <- AT_species_coverage %>% select(taxa_list[3],taxa_list[4]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- AT_species_coverage %>% group_by(Order) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Order,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
## Order-Family
tmp <- AT_species_coverage %>% select(taxa_list[4],taxa_list[5]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- AT_species_coverage %>% group_by(Family) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Family,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
## Family-Genus
tmp <- AT_species_coverage %>% select(taxa_list[5],taxa_list[6]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- AT_species_coverage %>% group_by(Genus) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Genus,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
## Genus-Species
tmp <- AT_species_coverage %>% select(taxa_list[6],taxa_list[7]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- AT_species_coverage %>% group_by(Species) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Species,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
### With networkD3, connection must be provided using id ###
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
# Sankey Diagram
## Add a 'group' column to the nodes data frame:
nodes$group <- as.factor(c(rep("d",length(unique(AT_species_coverage$Domain))),
                           rep("p",length(unique(AT_species_coverage$Phylum))),
                           rep("c",length(unique(AT_species_coverage$Class))),
                           rep("o",length(unique(AT_species_coverage$Order))),
                           rep("f",length(unique(AT_species_coverage$Family))),
                           rep("g",length(unique(AT_species_coverage$Genus))),
                           rep("s",length(unique(AT_species_coverage$Species)))
                           ))
# Add a 'group' column to each connection:
links$group <- as.factor(c(rep("p_link",length(unique(AT_species_coverage$Phylum))),
                           rep("c_link",length(unique(AT_species_coverage$Class))),
                           rep("o_link",length(unique(AT_species_coverage$Order))),
                           rep("f_link",length(unique(AT_species_coverage$Family))),
                           rep("g_link",length(unique(AT_species_coverage$Genus))),
                           rep("s_link",length(unique(AT_species_coverage$Species)))
                           ))
# Give a color for each group:
RColorBrewer::display.brewer.all()
display.brewer.pal(n=7,name="Set2")
brewer.pal(n=7,name="Set2")
display.brewer.pal(n=6,name="Pastel2")
brewer.pal(n=6,name="Pastel2")
my_color <- 'd3.scaleOrdinal() .domain(["d","p","c","o","f","g","s",
                                        "p_link","c_link","o_link","f_link","g_link","s_link"]) 
                                        .range(["#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
                                                "#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE"])'
## Plot
p<-sankeyNetwork(Links = links, Nodes= nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              colourScale=my_color, LinkGroup="group", NodeGroup="group",
              sinksRight=FALSE,fontSize = 13,fontFamily = "Arial")
print(p)
# save the widget
library(htmlwidgets)
#saveWidget(p, "../../airborne_arg_uwtp_result/Figure/Sankey/AT_top10_coverage_sankey.html")






# Aeration tank PM2.5 (ARP)
ARP_species_coverage <- final_ARG_coverage %>% filter(Sample_type == "ARP") %>% 
  filter(!(Species == "Unclassified")) %>% 
  group_by(Species,Sample) %>% 
  summarise(sum_coverage = sum(coverage))
ARP_species_coverage <-  ARP_species_coverage %>% ungroup() %>% group_by(Species) %>% 
  summarise(mean_coverage = mean(sum_coverage))
## join ARP_species_coverage with coverage
taxa <- final_ARG_coverage %>% select("Domain","Phylum","Class",
                                      "Order","Family","Genus","Species") %>% 
  unique()
ARP_species_coverage <- left_join(ARP_species_coverage,taxa,by="Species")
## Select top 10 species
ARP_species_coverage <- ARP_species_coverage %>% arrange(desc(mean_coverage))
ARP_species_coverage <- ARP_species_coverage[1:10,]
# Make Sankey node
taxa_list <- c("Domain","Phylum","Class",
               "Order","Family","Genus","Species")
tmp <- ARP_species_coverage %>% select(taxa_list[1]) %>% unique()
colnames(tmp) <- "name"
for (i in 2:length(taxa_list)) {
  tmp2 <- ARP_species_coverage %>% select(taxa_list[i]) %>% unique()
  colnames(tmp2) <- "name"
  tmp <- rbind(tmp,tmp2)
}
nodes <- tmp
# Make Sankey links
## Domain-Phylum
tmp <- ARP_species_coverage %>% select(taxa_list[1],taxa_list[2]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_species_coverage %>% group_by(Phylum) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Phylum,value) %>% unique()
colnames(tmp2)[1] <- "target"
links <- full_join(tmp,tmp2)
## Phylum-Class
tmp <- ARP_species_coverage %>% select(taxa_list[2],taxa_list[3]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_species_coverage %>% group_by(Class) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Class,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
## Class-Order
tmp <- ARP_species_coverage %>% select(taxa_list[3],taxa_list[4]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_species_coverage %>% group_by(Order) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Order,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
## Order-Family
tmp <- ARP_species_coverage %>% select(taxa_list[4],taxa_list[5]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_species_coverage %>% group_by(Family) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Family,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
## Family-Genus
tmp <- ARP_species_coverage %>% select(taxa_list[5],taxa_list[6]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_species_coverage %>% group_by(Genus) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Genus,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
## Genus-Species
tmp <- ARP_species_coverage %>% select(taxa_list[6],taxa_list[7]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_species_coverage %>% group_by(Species) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Species,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
### With networkD3, connection must be provided using id ###
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
# Sankey Diagram
## Add a 'group' column to the nodes dARPa frame:
nodes$group <- as.factor(c(rep("d",length(unique(ARP_species_coverage$Domain))),
                           rep("p",length(unique(ARP_species_coverage$Phylum))),
                           rep("c",length(unique(ARP_species_coverage$Class))),
                           rep("o",length(unique(ARP_species_coverage$Order))),
                           rep("f",length(unique(ARP_species_coverage$Family))),
                           rep("g",length(unique(ARP_species_coverage$Genus))),
                           rep("s",length(unique(ARP_species_coverage$Species)))
                           ))
# Add a 'group' column to each connection:
links$group <- as.factor(c(rep("p_link",length(unique(ARP_species_coverage$Phylum))),
                           rep("c_link",length(unique(ARP_species_coverage$Class))),
                           rep("o_link",length(unique(ARP_species_coverage$Order))),
                           rep("f_link",length(unique(ARP_species_coverage$Family))),
                           rep("g_link",length(unique(ARP_species_coverage$Genus))),
                           rep("s_link",length(unique(ARP_species_coverage$Species)))
                           ))
# Give a color for each group:
RColorBrewer::display.brewer.all()
display.brewer.pal(n=7,name="Set2")
brewer.pal(n=7,name="Set2")
display.brewer.pal(n=6,name="Pastel2")
brewer.pal(n=6,name="Pastel2")
my_color <- 'd3.scaleOrdinal() .domain(["d","p","c","o","f","g","s",
                                        "p_link","c_link","o_link","f_link","g_link","s_link"]) 
                                        .range(["#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
                                                "#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE"])'
## Plot
p <- sankeyNetwork(Links = links, Nodes= nodes,
                 Source = "IDsource", Target = "IDtarget",
                 Value = "value", NodeID = "name", 
                 colourScale=my_color, LinkGroup="group", NodeGroup="group",
                 sinksRight=FALSE,fontSize = 13,fontFamily = "Arial")
print(p)
# save the widget
#saveWidget(p, "../../airborne_arg_uwtp_result/Figure/Sankey/ARP_top10_coverage_sankey.html")


######## Barplot of sankey  which should be added in PPT manually ########
# AT
AT_bar <- final_ARG_coverage %>% filter(Species %in% AT_species_coverage$Species)
## Order species
AT_bar$Species <- factor(AT_bar$Species, 
                         levels = AT_species_coverage$Species[1:10])
# ARP
ARP_bar <- final_ARG_coverage %>% filter(Species %in% ARP_species_coverage$Species)
## Order species
ARP_bar$Species <- factor(ARP_bar$Species, 
                      levels = ARP_species_coverage$Species[1:10])
# Bind dataframe
taxa_bar <- rbind(AT_bar,ARP_bar)
# Edit ARG type
taxa_bar$type <- str_to_title(taxa_bar$type)
taxa_bar$type[taxa_bar$type == "Macrolide-Lincosamide-Streptogramin"] <- "MLS"
taxa_bar$type[taxa_bar$type == "Beta_lactam"] <- "Beta-lactam"
# Select color
library(RColorBrewer)
RColorBrewer::display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")
# Plot
p <- ggplot(taxa_bar, aes(x = Species, y = coverage, fill = type))+
  geom_bar(stat="identity", position = "fill",width = 0.5) +
  coord_flip()+
  scale_x_discrete(limits = rev(levels(taxa_bar$Species))) + 
  theme_bw()+
  #scale_fill_brewer(palette="Set3") +
  xlab("") + ylab("") +
  guides(fill=guide_legend(title="ARG type")) +
  scale_fill_manual(values=c("#8DD3C7","#FFED6F","#BEBADA","#FB8072","#80B1D3",
                             "#FDB462","#B3DE69","#FCCDE5", "#BC80BD","#D9D9D9")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) #transparent legend bg))
print(p) 

# Save
# ggsave("Sankey_barplot.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/Sankey",
#        width = 7, height = 7,
#        units = "in", bg='transparent') # save to png format



###############################################
################### Test ######################
# AT
## Class
AT_class_coverage <- final_ARG_coverage %>% filter(Sample_type == "AT") %>% 
  filter(!(Class == "Unclassified")) %>% 
  group_by(Class,Sample) %>% 
  summarise(sum_coverage = sum(coverage))
AT_class_coverage <-  AT_class_coverage %>% ungroup() %>% group_by(Class) %>% 
  summarise(mean_coverage = mean(sum_coverage))
AT_class_coverage <- AT_class_coverage %>%  filter(!(mean_coverage == 0))
## join AT_class_coverage with coverage
taxa <- final_ARG_coverage %>% select("Domain","Phylum","Class") %>% 
  unique()
AT_class_coverage <- left_join(AT_class_coverage,taxa,by="Class")
# Make Sankey node
taxa_list <- c("Domain","Phylum","Class")
tmp <- AT_class_coverage %>% select(taxa_list[1]) %>% unique()
colnames(tmp) <- "name"
for (i in 2:length(taxa_list)) {
  tmp2 <- AT_class_coverage %>% select(taxa_list[i]) %>% unique()
  colnames(tmp2) <- "name"
  tmp <- rbind(tmp,tmp2)
}
nodes <- tmp
# Make Sankey links
## Domain-Phylum
tmp <- AT_class_coverage %>% select(taxa_list[1],taxa_list[2]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- AT_class_coverage %>% group_by(Phylum) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Phylum,value) %>% unique()
colnames(tmp2)[1] <- "target"
links <- full_join(tmp,tmp2)
## Phylum-Class
tmp <- AT_class_coverage %>% select(taxa_list[2],taxa_list[3]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- AT_class_coverage %>% group_by(Class) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Class,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
### With networkD3, connection must be provided using id ###
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
# Sankey Diagram
## Plot
sankeyNetwork(Links = links, Nodes= nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              #colourScale=my_color, LinkGroup="group", NodeGroup="group",
              sinksRight=FALSE,fontSize = 15,fontFamily = "Arial")







# ARP
## Class
ARP_class_coverage <- final_ARG_coverage %>% filter(Sample_type == "ARP") %>% 
  filter(!(Class == "Unclassified")) %>% 
  group_by(Class,Sample) %>% 
  summarise(sum_coverage = sum(coverage))
ARP_class_coverage <-  ARP_class_coverage %>% ungroup() %>% group_by(Class) %>% 
  summarise(mean_coverage = mean(sum_coverage))
## join ARP_class_coverage with coverage
taxa <- final_ARG_coverage %>% select("Domain","Phylum","Class") %>% 
        unique()
ARP_class_coverage <- left_join(ARP_class_coverage,taxa,by="Class")
# Make Sankey node
taxa_list <- c("Domain","Phylum","Class")
tmp <- ARP_class_coverage %>% select(taxa_list[1]) %>% unique()
colnames(tmp) <- "name"
for (i in 2:length(taxa_list)) {
  tmp2 <- ARP_class_coverage %>% select(taxa_list[i]) %>% unique()
  colnames(tmp2) <- "name"
  tmp <- rbind(tmp,tmp2)
}
nodes <- tmp
# Make Sankey links
## Domain-Phylum
tmp <- ARP_class_coverage %>% select(taxa_list[1],taxa_list[2]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_class_coverage %>% group_by(Phylum) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Phylum,value) %>% unique()
colnames(tmp2)[1] <- "target"
links <- full_join(tmp,tmp2)
## Phylum-Class
tmp <- ARP_class_coverage %>% select(taxa_list[2],taxa_list[3]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_class_coverage %>% group_by(Class) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Class,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
### With networkD3, connection must be provided using id ###
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
# Sankey Diagram
## Plot
sankeyNetwork(Links = links, Nodes= nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   #colourScale=my_color, LinkGroup="group", NodeGroup="group",
                   sinksRight=FALSE,fontSize = 15,fontFamily = "Arial")

## Order
ARP_order_coverage <- final_ARG_coverage %>% filter(Sample_type == "ARP") %>% 
  filter(!(Order == "Unclassified")) %>% 
  group_by(Order,Sample) %>% 
  summarise(sum_coverage = sum(coverage))
ARP_order_coverage <-  ARP_order_coverage %>% ungroup() %>% group_by(Order) %>% 
  summarise(mean_coverage = mean(sum_coverage))
## join ARP_order_coverage with coverage
taxa <- final_ARG_coverage %>% select("Domain","Phylum","Class","Order") %>% 
  unique()
ARP_order_coverage <- left_join(ARP_order_coverage,taxa,by="Order")
# Make Sankey node
taxa_list <- c("Domain","Phylum","Class","Order")
tmp <- ARP_order_coverage %>% select(taxa_list[1]) %>% unique()
colnames(tmp) <- "name"
for (i in 2:length(taxa_list)) {
  tmp2 <- ARP_order_coverage %>% select(taxa_list[i]) %>% unique()
  colnames(tmp2) <- "name"
  tmp <- rbind(tmp,tmp2)
}
nodes <- tmp
# Make Sankey links
## Domain-Phylum
tmp <- ARP_order_coverage %>% select(taxa_list[1],taxa_list[2]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_order_coverage %>% group_by(Phylum) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Phylum,value) %>% unique()
colnames(tmp2)[1] <- "target"
links <- full_join(tmp,tmp2)
## Phylum-Class
tmp <- ARP_order_coverage %>% select(taxa_list[2],taxa_list[3]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_order_coverage %>% group_by(Class) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Class,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
## Class-Order
tmp <- ARP_order_coverage %>% select(taxa_list[2],taxa_list[3]) %>% unique()
colnames(tmp) <- c("source","target")
tmp2 <- ARP_order_coverage %>% group_by(Order) %>% mutate(value = sum(mean_coverage)) %>% 
  select(Order,value) %>% unique()
colnames(tmp2)[1] <- "target"
links_2 <- full_join(tmp,tmp2)
links <- rbind(links,links_2)
### With networkD3, connection must be provided using id ###
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
# Sankey Diagram
## Plot
sankeyNetwork(Links = links, Nodes= nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              #colourScale=my_color, LinkGroup="group", NodeGroup="group",
              sinksRight=FALSE,fontSize = 15,fontFamily = "Arial")
