# ARG_UpSet.R

# Import library
library(tidyverse)
library(UpSetR)
library(openxlsx)
library(RColorBrewer)
# Read ARG_subtype file (Something wrong with read.table, so read by )
arg_subtype <- read.xlsx("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.subtype.xlsx",
                         sheet = 1)
# First column to row name
row.names(arg_subtype) <- arg_subtype$subtype
arg_subtype <- arg_subtype[,-1]
# Replace no 0 value with ARG name
for (i in c(1:987)) {
  for (j in c(1:15)) {
    if (!(arg_subtype[i,j] == 0))
    {arg_subtype[i,j] <- row.names(arg_subtype[i,])}
  }
}
# Replace no 0 value with NA
arg_subtype[arg_subtype == 0] <- NA
# Select fill color
display.brewer.all()
brewer.pal(3, "Pastel1")
# Plot
## Save
pdf(file= "../../airborne_arg_uwtp_result/Figure/ARG/ARG_upset.pdf",
    width = 9,
    height = 7)
upset(fromList(arg_subtype),
      nsets = 15,
      nintersects = 52,
      order.by = "freq",
      mb.ratio = c(0.6,0.4),   # 左侧和上方条形图的比例关系
      text.scale = 1.5,
      mainbar.y.label = "Intersection ARG subtype numbers",
      sets.x.label = "ARG subtype numbers",
      main.bar.color = "#f4cccc", matrix.color = "#A6CEE3",
      sets.bar.color = "#A6CEE3",
      shade.color = "#A6CEE3",
      shade.alpha = 0.25, matrix.dot.alpha = 0.5,
      queries = list(
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP1","ODP2","ODP3","ODP4","ODP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP1","ODP2","ODP4","ODP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP2","ODP3","ODP4","ODP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP4"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP1","ODP3","ODP4","ODP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP3","ARP4","ARP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP2","ODP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP2","ODP3","ODP4"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP3","ARP4","ARP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP3"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP4","ODP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT2",
                          "ARP1"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT4",
                          "ARP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP5",
                          "ODP5"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP1"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP2"),
            color = "#FB8072",active = T),
          list(
            query = intersects, 
            params = list("AT1","AT2","AT3","AT4","AT5",
                          "ARP1","ARP2","ARP3","ARP4","ARP5",
                          "ODP5"),
            color = "#FB8072",active = T)
          )
      )
dev.off()
