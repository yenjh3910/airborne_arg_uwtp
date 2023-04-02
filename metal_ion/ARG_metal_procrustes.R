# ARG_metal_procrustes.R
library(tidyverse)
library(openxlsx)
library(vegan)
# Import & prepare input dataframe
metal <- read.xlsx("../../airborne_arg_uwtp_result/metal_ion/UWTP_metal.xlsx", sheet = 1)
row.names(metal) <- metal[,1]
metal <- metal[,-1]
#metal <- as.data.frame(t(metal)) 
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx", sheet = 1)
row.names(arg_subtype) <- arg_subtype[,1]
arg_subtype <- arg_subtype[,-1]
arg_subtype <- as.data.frame(t(arg_subtype)) 
# Remove wastewater sample
metal <- metal[-(11:15),]
arg <-arg_subtype[-(6:10),]
# PCA for metal
metal.std=decostand(metal, method = "standardize")
pca.metal <- rda(metal.std)
plot(pca.metal, scaling =1, display = "sites", type = "text",
     main = "PCA for Metal Data")
# PCA for ARG
arg.hel=decostand(arg, method = "hellinger")
arg_bray<-vegdist(arg.hel, method="bray")
pcoa.arg = cmdscale(arg_bray, eig=TRUE)
pro.g.s<-procrustes(pca.metal,pcoa.arg,symmetric = T, scores = "sites", choice = c(1,2))
## Ckeck statistic result
summary(pro.g.s)
plot(pro.g.s, kind = 1,type="text")
plot(pro.g.s, kind = 2)
## Check statistic result by protest
prot <- protest(X = pca.metal, Y = pcoa.arg, permutations = 999)
prot
names(prot)
prot$signif  # p value
prot$ss  # M2

Y<-cbind(data.frame(pro.g.s$Yrot),data.frame(pro.g.s$X))
X<-data.frame(pro.g.s$rotation)
Y$sample_type<-rownames(Y)
Y$sample_type <- gsub("1|2|3|4|5","",Y$sample_type)
#color
library(RColorBrewer)
RColorBrewer::display.brewer.all()
display.brewer.pal(n=12,name="Set3")
brewer.pal(n=12,name="Set3")
color<-c("#FB8072","#80B1D3")


p <- ggplot(Y) +
  geom_segment(aes(x = X1, y = X2,
                   xend = PC1, yend = PC2, color=sample_type),
               # geom_segment 绘制两点间的直线
               size = 0.75,linetype="dashed",alpha=0.7) +
  geom_point(aes(X1, X2, color =sample_type),shape=16,size = 3,alpha=0.5) +
  geom_point(aes(PC1,PC2,color = sample_type),shape=17,size = 3,alpha=0.5) +
  scale_color_manual("Sample type", values = color, 
                     labels = c(c(expression(Aeration~tank~PM[2.5]), 
                                  expression(Outdoor~PM[2.5])))) +
  guides(color = guide_legend(label.hjust = 0)) + 
  theme_bw() + labs(title="Procrustes analysis") +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.3) +
  annotate('text', label = sprintf('M^2 == 0.4486 '),
           x = 0.30, y = 0.37, size = 5, parse = TRUE) +
  annotate('text', label = 'P==0.002',
           x = 0.30, y = 0.32, size = 5, parse = TRUE) +
  theme(title = element_text(size=14),
        axis.title = element_text(size=15),
        legend.title= element_text(size=13),
        legend.text = element_text(size=13),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.background = element_rect(fill='transparent')) #transparent legend bg)

print(p)

# # Save
# ggsave("ARG_metal_procrustes.png", p, 
#        path = "../../airborne_arg_uwtp_result/Figure/procrustes",
#        width = 7, height = 5, 
#        units = "in", bg='transparent') # save to png format