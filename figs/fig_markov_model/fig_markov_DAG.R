library(igraph)
library(tidyverse)
library(readxl)

Clinical_states <- read_excel("figs/fig_markov_model/Clinical_states.xlsx", col_names = FALSE)

nodes <- Clinical_states
nodes

links <- c(from = nodes[, 1], to = nodes[, 2])
links

net <- graph_from_data_frame(d = links, directed = TRUE)

set.seed(1)
plot(net, layout = layout_as_tree, edge.arrow.size=.05, edge.arrow.width = 0.75,
     vertex.color = alpha("lightgrey", .95), vertex.label.color = "black", 
     vertex.frame.color = alpha("darkgrey", .35), vertex.label.font = 4, vertex.size = 13, 
     vertex.label.cex = 0.4, edge.curved = 0, edge.color = alpha("black", .75), 
     mark.groups = 1, mark.col = alpha("lightgreen", .025))
# title(main = "Human Papillomavirus: a Markov Model of its natural history")

# End file ----------------------------------------------------------------