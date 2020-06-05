library(igraph)
library(tidyverse)
library(readxl)

Clinical_states <- read_excel("Clinical_states.xlsx", col_names = FALSE)

nodes <- Clinical_states
nodes

links <- c(from = nodes[, 1], to = nodes[, 2])
links

net <- graph_from_data_frame(d = links, directed = TRUE)

set.seed(25)
plot(net, layout = layout.davidson.harel, edge.arrow.size=.12, edge.arrow.width = 0.65,
     vertex.color = alpha("lightgrey", .25), vertex.label.color = "black", 
     vertex.label.font = 4, vertex.size = 27, vertex.label.cex = 0.75, edge.curved = .075,
     edge.color = "black", mark.groups = 1, mark.col = alpha("lightgreen", .2))
title(cex.main = 1.1, main = "Human Papillomavirus:
      a model of its health states and possible transitions")


# End file ----------------------------------------------------------------