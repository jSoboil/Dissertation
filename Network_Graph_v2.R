# library
library(igraph)
library(visNetwork)
library(networkD3)
library(readr)
library(readxl)
library(RColorBrewer)
library(ggplot2)

# Load data:
nodes <- data.frame(id = 1:17, label = c("He", "Expo", "Infec", "Imm",
                                         "De", "Geni", "Clear", "Reinfect",
                                         "LSIL", "HSIL-1", "HSIL-2", "Cancer", 
                                         "Year 1", "Year 2", "Year 3", "Clear", 
                                         "Reinfect"), shadow = TRUE)

edges <- data.frame(from = c(1, 1, 2, 2, 3, 2, 3, 4, 3, 6, 6, 7, 7, 8, 8, 3, 3, 9, 10,
                             9, 10, 11, 11, 9, 10, 11, 11, 12, 13, 14, 15, 16, 16,
                             16, 17, 17, 17),
                    to = c(1, 2, 2, 3, 2, 4, 4, 5, 6, 6, 7, 7, 8, 7, 6, 9, 10, 10, 9,
                           11, 11, 9, 10, 9, 10, 11, 12, 13, 14, 15, 16, 16, 16,
                           17, 16, 9, 10), smooth = FALSE,  length = 50, width = 2.5)

visNetwork(nodes, edges, shadow = TRUE, width = "120%") %>% 
 visNodes(size = 20, color = alpha("darkgrey", .9)) %>% 
 visEdges(shadow = TRUE, arrows = c("to"), smooth = FALSE,
          color = list(color = alpha("black", .85), 
                       highlight = "red")) %>%
 visIgraphLayout(randomSeed = 249, layout = "layout.davidson.harel") %>% 
 visClusteringByConnection(nodes = nodes)
# Seed options: 246 or 234.
