library(igraph)
library(tidyverse)
library(readxl)

Clinical_states <- read_excel("Markov_DAG/Clinical_states.xlsx", col_names = FALSE)

nodes <- Clinical_states
nodes

links <- c(from = nodes[, 1], to = nodes[, 2])
links

net <- graph_from_data_frame(d = links, directed = TRUE)

# 57
# 58
# 59
# 10
set.seed(58)
plot(net, layout = layout.fruchterman.reingold, edge.arrow.size=.18, edge.arrow.width = 0.75,
     vertex.color = alpha("navyblue", .3), vertex.label.color = "black", 
     vertex.size = 19, vertex.label.cex = 0.75, edge.curved = .075,
     edge.color = "black")

# End file ----------------------------------------------------------------