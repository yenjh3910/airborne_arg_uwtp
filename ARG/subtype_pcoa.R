# subtype_pcoa.R

library(tidyverse)
library(openxlsx)
theme_set(theme_test())
library(ape)
library(ade4)

# Import
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
# Adjust dataframe
subtype_data <- as.data.frame(t(arg_subtype))
colnames(subtype_data) <- subtype_data[1,]
subtype_data <- subtype_data[-1,]
subtype_data<- subtype_data %>% mutate(sample_type = row.names(subtype_data)) %>% 
                               relocate(sample_type)
row.names(subtype_data) <- 1:length(subtype_data[,1])
subtype_data$sample_type <- gsub("1|2|3|4|5","",subtype_data$sample_type)


pcoa.data <- iris %>% mutate(Sepal.Length = scale(Sepal.Length),
                             Sepal.Width = scale(Sepal.Width),
                             Petal.Length = scale(Petal.Length),
                             Petal.Width = scale(Petal.Width))

tab.dist <- dist(subtype_data[,-1])
pcoa <- pcoa(tab.dist)
pcoa_eig <- pcoa$values

sample_site <- pcoa$vectors%>%as.data.frame()

names(sample_site)[1:2] <- c("x","y")

sample_site <- data.frame(sample_site,subtype_data$sample_type)

sample_site%>%ggplot()+
  geom_point(aes(x,y,color=subtype_data.sample_type),size = 1.5)+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  stat_ellipse(aes(x,y,color = subtype_data.sample_type,fill = subtype_data.sample_type),
               geom ="polygon",level = 0.95,size = 0.5,alpha = 0.2)+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1,2], 2), '%'), 
       y = paste('PCoA2: ', round(100 * pcoa_eig[2,2], 2), '%'))


#=================================================================
library(tidyverse)
library(openxlsx)
library(vegan)
data("varespec")
df<-varespec
rownames(df)<-paste0("site",1:24)

arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
# Adjust dataframe
subtype_data <- as.data.frame(t(arg_subtype))
colnames(subtype_data) <- subtype_data[1,]
subtype_data <- subtype_data[-1,]
df <- subtype_data
as.numeric(df)

bray_dist<-vegdist(df,method = "bray")
library(ape)
df.pcoa<-pcoa(bray_dist,correction = "cailliez")
df.plot<-data.frame(df.pcoa$vectors)
head(df.plot)
library(ggplot2)
x_label<-round(df.pcoa$values$Rel_corr_eig[1]*100,2)
y_label<-round(df.pcoa$values$Rel_corr_eig[2]*100,2)
x_label
y_label
ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2))+
  geom_point()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),
       y=paste0("PCoA2 ",y_label,"%"))
df.plot$group<-ifelse(df.plot$Axis.1<0,"AAA","BBB")
ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2,
                        color=group,shape=group))+
  geom_point(size=5)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),
       y=paste0("PCoA2 ",y_label,"%"))+
  stat_ellipse(data=df.plot,
               geom = "polygon",
               aes(fill=group),
               alpha=0.3)+
  scale_fill_manual(values = c("#e31a1c","#1f78b4"))










## Import and join ARG subtype dataset
library(tidyverse)
library(openxlsx)
library(vegan)
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
download_fq <- read.table("../../airborne_arg_uwtp_result/download_fq/nineyear_args_oap/normalized_cell.subtype.txt",
                header = TRUE, sep = "\t", quote = "")
colnames(download_fq) <- c('subtype',paste0('AS',1:(length(download_fq)-1)))
arg_subtype <- full_join(arg_subtype, download_fq)

# ##commamox
# download_fq <- read.table("../../airborne_arg_uwtp_result/download_fq/commamox/normalized_cell.subtype.txt",
#                                                      header = TRUE, sep = "\t", quote = "")
# colnames(download_fq) <- c('subtype',paste0('commamox',1:(length(download_fq)-1)))
# arg_subtype <- full_join(arg_subtype, download_fq)

#PRJEB14051
download_fq <- read.table("../../airborne_arg_uwtp_result/download_fq/PRJEB14051/normalized_cell.subtype.txt",
                          header = TRUE, sep = "\t", quote = "")
colnames(download_fq) <- c('subtype',paste0('PRJEB14051',1:(length(download_fq)-1)))
arg_subtype <- full_join(arg_subtype, download_fq)

#PRJNA385831
download_fq <- read.table("../../airborne_arg_uwtp_result/download_fq/PRJNA385831/normalized_cell.subtype.txt",
                          header = TRUE, sep = "\t", quote = "")
colnames(download_fq) <- c('subtype',paste0('PRJNA385831',1:(length(download_fq)-1)))
arg_subtype <- full_join(arg_subtype, download_fq)

## Effluent
# download_fq <- read.table("../../airborne_arg_uwtp_result/download_fq/effluent_args_oap_1/normalized_cell.subtype.txt",
#                           header = TRUE, sep = "\t", quote = "")
# download_fq2 <- read.table("../../airborne_arg_uwtp_result/download_fq/effluent_args_oap_2/normalized_cell.subtype.txt",
#                           header = TRUE, sep = "\t", quote = "")
# download_fq3 <- read.table("../../airborne_arg_uwtp_result/download_fq/effluent_args_oap_3/normalized_cell.subtype.txt",
#                            header = TRUE, sep = "\t", quote = "")
# download_fq <- full_join(download_fq, download_fq2)
# download_fq <- full_join(download_fq, download_fq3)
# colnames(download_fq) <- c('subtype',paste0('EFF&INF',1:(length(download_fq)-1)))
# arg_subtype <- full_join(arg_subtype, download_fq)

arg_subtype[is.na(arg_subtype)] <- 0
row.names(arg_subtype) <- arg_subtype[,1]
arg_subtype <- arg_subtype[,-1]
## Transform dataframe
arg_subtype <-as.data.frame(t(arg_subtype))
## Normalization
arg_subtype <- decostand(arg_subtype, method = 'hellinger')
arg_bray<-vegdist(arg_subtype, method="bray")
library(ape)
df.pcoa<-pcoa(arg_bray,correction = "cailliez")
df.plot<-data.frame(df.pcoa$vectors)
head(df.plot)

df.plot$sample_type <- gsub("[0-9]","",row.names(df.plot))

ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2,color=sample_type))+
  geom_point()+
  theme_bw()
