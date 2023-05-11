# ARG_type_circos.R

# Import library
library(tidyverse)
library(circlize)
library(RColorBrewer)

# Read ARG_type file
arg_type <- read.table("../../airborne_arg_uwtp_result/args_oap/ARG/stage_two_output/normalized_cell.type.txt",
                       header = TRUE, sep = "")
# gather df
gather_arg_type <- gather(arg_type, key = "sample", value = "copy_per_cell", 
                          ARP1:ODP5) # Transform to gather format
# Add sample type column
gather_arg_type$sample_type <- gather_arg_type$sample
gather_arg_type$sample_type <- gsub("1|2|3|4|5","",gather_arg_type$sample_type)
# Calculate average abundance per sample type
gather_arg_type <- gather_arg_type %>% group_by(type, sample_type) %>% 
  mutate(mean = mean(copy_per_cell)) %>% 
  select(!(sample)) %>% 
  select(!(copy_per_cell)) %>% 
  unique()
# Calculate mean value and order of ARG type in all sample
average_arg <- gather_arg_type %>% group_by(type) %>% summarize(average =  mean(mean))
arg_order <- average_arg[order(average_arg$average, decreasing = T),]$type
# Spread df
arg_circos <- gather_arg_type %>% spread(key = sample_type, value = mean)
# Order ARG type
arg_circos <- arg_circos %>% 
  arrange(factor(type, levels = arg_order))
# Convert to dataframe
arg_circos <- as.data.frame(arg_circos)
# First column as row name
row.names(arg_circos) <- arg_circos$type
arg_circos <- arg_circos[,-1]
# Convert to matrix
arg_circos<- as.matrix(arg_circos)
# Calculate others ARG type
Others <- colSums(arg_circos[11:26, ])
arg_circos <- rbind(arg_circos, Others)
# Delete minimum arg type
arg_circos <- arg_circos[-(11:26),]
# Covert first letter to uppercase
rownames(arg_circos) <- str_to_title(rownames(arg_circos))
# Change specific ARG type
row.names(arg_circos)[3] <- "MLS"
row.names(arg_circos)[6] <- "Beta-lactam"
# Change column order
arg_circos <- arg_circos[, c("AT", "ARP", "ODP")] 
# Set color
display.brewer.all()
brewer.pal(n=6,name="Set2")
grid.col = c("#8DD3C7","#FFFFB3","#BEBADA",
             "#FB8072","#80B1D3","#FDB462",
             "#B3DE69","#FCCDE5","#FFED6F",
             "#BC80BD","#D9D9D9", # ARG type
             "#FBB4AE","#B3CDE3","#CCEBC5" # sample type 
            )


# Save
#pdf("../../airborne_arg_uwtp_result/Figure/ARG/ARG_circos.pdf")
graphics.off() # Reset figure window
circos.par(canvas.xlim = c(-0.5, 0.5), canvas.ylim = c(-1.5, 1))
chordDiagram(arg_circos, annotationTrack = "grid",
             preAllocateTracks = 1, grid.col = grid.col, 
             directional = 1, big.gap = 10, scale = 0)
circos.trackPlotRegion(track.index = 1, 
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         sector.name = get.cell.meta.data("sector.index")
                         circos.text(cex = 0.8,mean(xlim), ylim[1] + .5, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
                         circos.axis(h = "top", labels.cex = 0.4, major.tick.percentage = 1, sector.index = sector.name, track.index = 2)
                       }, bg.border = NA)
#dev.off()

#########################################################################

# Convert abundance to relative abundance
arg_circos <- apply(arg_circos,-1, function(x) x/sum(x))
# Save
graphics.off() # Reset figure window
# pdf("../../airborne_arg_uwtp_result/Figure/ARG/ARG_relative_circos.pdf")
circos.par(canvas.xlim = c(-0.5, 0.5), canvas.ylim = c(-1.5, 1))
chordDiagram(arg_circos, annotationTrack = "grid",
             preAllocateTracks = 1, grid.col = grid.col, 
             directional = 1, big.gap = 10, scale = 0)
circos.trackPlotRegion(track.index = 1, 
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         sector.name = get.cell.meta.data("sector.index")
                         circos.text(cex = 0.8,mean(xlim), ylim[1] + .5, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
                         circos.axis(h = "top", labels.cex = 0.4, major.tick.percentage = 1, sector.index = sector.name, track.index = 2)
                       }, bg.border = NA)
# dev.off()