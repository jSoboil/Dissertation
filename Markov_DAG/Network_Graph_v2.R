# library
library(igraph)
library(visNetwork)
library(networkD3)
library(readr)
library(readxl)
library(RColorBrewer)
library(ggplot2)

# Load data:
nodes <- data.frame(id = 1:15, label = c("Healthy", "Exposure", "Infection", 
                                        "Condyloma", "Clearance", "Reinfection", 
                                        "LSIL", "HSIL", "Clearance", "Reinfection",
                                        "Cervical Cancer", "Year 1", "Year 2", "Year 3", 
                                        "Death"), 
                    shadow = TRUE,  font.color = "black", font.size = 25, 
                    font.bold = TRUE) 

edges <- data.frame(from = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 3, 7, 3, 7, 8, 8, 7, 8, 
                             9, 10, 10, 10, 8, 11, 12, 13, 14, 15),
                    to = c(1, 2, 2, 3, 2, 4, 4, 5, 5, 6, 5, 4, 7, 7, 8, 8, 7, 8, 9, 9, 
                           10, 9, 7, 8, 11, 12, 13, 14, 9, 15), 
                    smooth = TRUE,  length = 1, width = 1.5)

visNetwork(nodes, edges, shadow = FALSE, width = "100%", height = "350px") %>% 
 visNodes(size = 22, color = alpha("black", .85), borderWidth = 3) %>% 
 visEdges(selfReferenceSize = 24, shadow = FALSE, smooth = TRUE,
          color = list(color = alpha("darkgrey", 1)), 
          arrows = list(to = list(scaleFactor = 1.5))) %>%
 visIgraphLayout(layout = "layout_nicely", smooth = TRUE, randomSeed = 22)
# Seed options: 22, 50.

# End file ----------------------------------------------------------------