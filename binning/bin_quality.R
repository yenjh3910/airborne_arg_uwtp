# bin_quality.R

# Import
library(tidyverse)
library(ggpubr)
bin_quality <- read.table("../../airborne_arg_uwtp_result/metawrap_bin/bin_refinement/metawrap_50_10_bins.stats",
                          header = TRUE)
# Quality category
highq_bin <- bin_quality %>% filter(completeness>=80 & contamination<=5) %>% 
                select(bin) %>% 
                mutate(Quality = "High")
medq_bin <- bin_quality %>% filter(completeness>=60 & contamination<=7) %>% 
                            select(bin) %>% 
                            mutate(Quality = "Medium")
medq_bin <- medq_bin[!(medq_bin$bin %in% highq_bin$bin),]

lowq_bin <- bin_quality[!(bin_quality$bin %in% highq_bin$bin),]
lowq_bin <- lowq_bin[!(lowq_bin$bin %in% medq_bin$bin),] %>% 
            select(bin) %>% 
            mutate(Quality = "Low")
bin_category <- rbind(highq_bin,medq_bin)
bin_category <- rbind(bin_category,lowq_bin)
bin_quality <- full_join(bin_quality,bin_category)
# Order sample type
bin_quality$Quality <- factor(bin_quality$Quality, levels = c("High","Medium","Low"))
# Select color
library(RColorBrewer)
RColorBrewer::display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")
# Plot
p <- ggscatter(bin_quality, x = "completeness", y = "contamination",
          color="Quality",alpha = 0.5) + 
  labs(x = "Completeness (%)", y = "Contamination (%)") + 
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        aspect.ratio=1, 
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) + #transparent legend bg
  scale_color_manual(values = c("High" = "#F8766D",
                                "Medium" = "#619CFF", 
                                "Low" = "#00BA38"))

print(p)
# ggsave("MAG_quality.png", p, path = "../../airborne_arg_uwtp_result/Figure/binning",
#        width = 4, height = 4, units = "in", bg='transparent') # save to png format
# Density
p <- ggplot(bin_quality, aes(x=completeness)) + 
  geom_density(color="#FCCDE5", fill="#FCCDE5",alpha=0.6) +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA), 
        legend.background = element_rect(fill='transparent'))
print(p)
# ggsave("completeness_density.png", p, path = "../../airborne_arg_uwtp_result/Figure/binning",
#        width = 2.63, height = 1, units = "in", bg='transparent') # save to png format
p <- ggplot(bin_quality, aes(x=contamination)) + 
  geom_density(color="#BEBADA", fill="#BEBADA",alpha=0.6) +
  coord_flip() +
  theme_classic() +
  theme(panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA), 
        legend.background = element_rect(fill='transparent'))
print(p)
# ggsave("contamination_density.png", p, path = "../../airborne_arg_uwtp_result/Figure/binning",
#        width = 1.2, height = 2.5, units = "in", bg='transparent') # save to png format
