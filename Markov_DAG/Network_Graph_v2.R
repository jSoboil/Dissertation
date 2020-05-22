# library
library(igraph)
library(visNetwork)
library(networkD3)
library(readr)
library(readxl)
library(RColorBrewer)
library(ggplot2)

# Load data:
nodes <- data.frame(id = 1:12, label = c("Healthy", "Exposure", "Infection",
                                        "LSIL", "HSIL", "Cervical Cancer", "Year 1", 
                                        "Year 2", "Year 3", "Clearance", "Reinfection", 
                                        "Death"), 
                    shadow = FALSE,  font.color = "black", font.size = 30, 
                    font.bold = TRUE) 

edges <- data.frame(from = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 7, 8, 9, 4, 
                             5, 10, 11, 11, 11, 12),
                    to = c(1, 2, 2, 3, 4, 5, 4, 5, 4, 5, 6, 7, 8, 9, 10, 10, 
                           10, 11, 10, 4, 5, 12), 
                    smooth = FALSE,  length = 1, width = 2)

visNetwork(nodes, edges, shadow = FALSE, width = "100%", height = "350px") %>% 
 visNodes(size = 18, color = alpha("grey", .85), borderWidth = ) %>% 
 visEdges(selfReferenceSize = 24, shadow = TRUE, smooth = FALSE,
          color = list(color = alpha("black", 1)), 
          arrows = list(to = list(scaleFactor = .85))) %>%
 visIgraphLayout(layout = "layout_nicely", smooth = TRUE, randomSeed = 2)
# Seed options: 22, 50.

# End file ----------------------------------------------------------------