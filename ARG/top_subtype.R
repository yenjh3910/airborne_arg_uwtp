# top_subtype.R

source("ARG_subtype_cell_abundance.R")
# Order by abundance
subtype_mean_sd_notlog <- subtype_mean_sd_notlog[order(subtype_mean_sd_notlog$mean, 
                                                     decreasing = T),]
# Adjust "_" in dataframe
subtype_mean_sd_notlog <- subtype_mean_sd_notlog %>% 
  separate(subtype, c("type","subtype"), 
           sep = "__") # Split subtype by "_"
# Replace subtype name 
subtype_mean_sd_notlog$subtype <- gsub("Bifidobacteria intrinsic ileS conferring resistance to mupirocin", 
                            "Bifidobacteria intrinsic ileS", subtype_mean_sd_notlog$subtype)
# Print top subtype
TOP_ARG_NUM <- 20
unique(subtype_mean_sd_notlog$subtype)[1:TOP_ARG_NUM]
# Select top subtype
top_subtype_list <- unique(subtype_mean_sd_notlog$subtype)[1:TOP_ARG_NUM]
# Assign top subtype dataframe
top_subtype <- filter(subtype_mean_sd_notlog, subtype %in% top_subtype_list)
# Order
top_subtype <- top_subtype %>% group_by(subtype) %>% 
  mutate(sum = sum(mean)) %>% 
  arrange(desc(sum))
top_subtype$subtype <- factor(top_subtype$subtype, levels=unique(top_subtype$subtype))
top_subtype$sample_type <- factor(top_subtype$sample_type, levels=c("AT","ARP","ODP"))

# Add st_mean for figure
top_subtype <- top_subtype %>% group_by(sample_type) %>% 
  mutate(sd_mean = mean)
top_subtype <- top_subtype %>% arrange(subtype,match(sample_type, c('ODP','ARP','AT')))
## sum st_mean
for (j in 0:(TOP_ARG_NUM-1)) {
  for (i in 2:(nrow(top_subtype)/TOP_ARG_NUM)) {
    k <- i+j*3
    top_subtype$sd_mean[k] <-  top_subtype$sd_mean[k-1] + top_subtype$sd_mean[k]
  }
}
## Change legend title
top_subtype$sample_type <- factor(top_subtype$sample_type, levels = c("AT", "ARP", "ODP"), 
                                      labels = c(expression(Aeration~tank), 
                                                 expression(Aeration~tank~PM[2.5]), 
                                                 expression(Outdoor~PM[2.5])
                                      ))
# Plot
p<-ggplot(top_subtype, aes(x = subtype, y = mean, fill = sample_type))+
  geom_bar(stat="identity",alpha=0.6)+
  theme_bw()+
  xlab("")+ ylab("Relative abundance (ARGs/cell)")+
  scale_fill_manual(values=c("#F8766D", "#00BFC4", "#7CAE00"),
                    name = 'Sample', 
                    labels = c(expression(Aeration~tank),
                               expression(Aeration~tank~PM[2.5]),
                               expression(Outdoor~PM[2.5])))+
  geom_errorbar(aes(x=subtype, ymin=sd_mean-sd, ymax=sd_mean+sd), 
                 width=0.4, colour="black", alpha=0.9, linewidth=0.2)+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust=1),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.text.align = 0,
        legend.position = c(0.75, 0.8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        #legend.key.size = unit(0.8, 'cm'),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill='transparent'))
print(p)
# ggsave("top_subtype.png", p,
#        path = "../../airborne_arg_uwtp_result/Figure/ARG",
#        width = 6, height = 5,
#        units = "in", bg='transparent') # save to png format
