# library
library(igraph)
library(networkD3)
library(readr)
library(readxl)
library(RColorBrewer)
library(ggplot2)

# Load data:
data <- read_xlsx(path = "Clinical_states.xlsx", sheet = 2, 
                  col_names = c("Source", "Target"))

# Turn it into igraph object
network <- graph_from_data_frame(d = data, directed = TRUE)

# Count the number of degree for each node:
deg <- degree(network, mode = "all")

# Set plot options:
par(mar = c(0, 0, 0, 0), bg = "white")

# Open jpeg file:
# pdf("Health_States.pdf")

# Plot network object:
set.seed(286)
myPlot <- plot.igraph(network, edge.arrow.size = .35, edge.width = 1.5, vertex.size = 14.5, 
     edge.curved = .062, edge.color = alpha("black", 1), vertex.label.color = "black", 
     vertex.label.font = 4, vertex.label.cex = .65, vertex.color = alpha("white", .8),
     vertex.frame.color = alpha("darkred", .75), layout = layout.davidson.harel,
     mark.groups = c(15), mark.col = alpha("forestgreen", .1), 
     mark.border = alpha("darkred", 1)
)
title("Human Papillomavirus", line = - 5, adj = .6, cex.main = 2)
text(x = .01, y = .4, "Note: Death state reached from\n every state in model. \n
        Omitted transition arrows \nfor ease of reading. \n 
     Adapted from Favato G., et al., Medical Care, 2012", 
     pos = 2, cex = .5, offset = -19, font = 2)
legend(x = .3, y = .8, legend = c("c = Cervical", 
                                  "ci = Cervical Intraepithelial Neoplasia",
                                  "w = Warts"), 
       cex = .5, box.lwd = 2)

# dev.off() 
myPlot
# End file ----------------------------------------------------------------