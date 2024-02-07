# metacompare_plot.R

library(tidyverse)
library(scales)
library(scatterplot3d)
risk <- read.csv("../../airborne_arg_uwtp_result/metacompare/risk_score.csv")
select_risk <- risk %>% select(sample_type,
                               nARG.nContigs,
                               nARG.MGE.nContigs,
                               nARG.MGE.PAT.nContigs,)
# shape,color,tick
shapes = c(17, 16, 18) 
shapes <- shapes[as.factor(select_risk$sample_type)]
colors <- c("#80B1D3","#FB8072","#FDB462")
colors <- colors[as.factor(select_risk$sample_type)]
x_ticks <- c(0,sprintf("%.e",c(0.001,0.002,0.003,0.004,0.005)))
y_ticks <- c(0,sprintf("%.e", c(0.00005,0.00010,0.00015,0.00020,0.00025,0.00030,0.00035)))
z_ticks <- c(0,sprintf("%.e", c(0.00005,0.00010,0.00015,0.00020,0.00025,0.00030,0.00035)))
# Plot
s3d <-scatterplot3d(select_risk[,2:4], angle = 55,
              xlab = "Q(ARG)",
              ylab = "Q(ARG,MGE)",
              zlab = "Q(ARG,MGE,PATH)",
              pch = shapes,
              color = colors,
              type="h",
              x.ticklabs = x_ticks, y.ticklabs = y_ticks, z.ticklabs = z_ticks,
              cex.axis=0.7)
legend("top", legend = levels(as.factor(select_risk$sample_type)),
       col =  c("#80B1D3","#FB8072","#FDB462"), pch = c(17,16,18))
# Statistic
risk %>% select(sample_type,Risk.Score) %>% 
         group_by(sample_type) %>% 
         mutate(mean = mean(Risk.Score)) %>% 
         mutate(sd = sd(Risk.Score)) %>% 
         select(!(Risk.Score)) %>% 
         unique()

# Menuscript version
colors <- c("#00BFC4", "#F8766D", "#7CAE00")
colors <- colors[as.factor(select_risk$sample_type)]
s3d <-scatterplot3d(select_risk[,2:4], angle = 55,
                    xlab = "Q(ARG)",
                    ylab = "Q(ARG,MGE)",
                    zlab = "Q(ARG,MGE,PATH)",
                    pch = shapes,
                    color = colors,
                    type="h",
                    x.ticklabs = x_ticks, y.ticklabs = y_ticks, z.ticklabs = z_ticks,
                    cex.axis=0.7)
legend("top", legend = levels(as.factor(select_risk$sample_type)),
       col =  c("#00BFC4", "#F8766D", "#7CAE00"), pch = c(17,16,18))
