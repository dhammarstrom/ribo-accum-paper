##### Figure source ###############################



# Libraries for figures ####

library(tidyverse)
library(cowplot)

# Color scale 

multiple.color <- "#636363"
single.color <- "#bdbdbd"

# line sizes
line_size <- 0.3

# text sizes
label.size <- 10

text.size <- 8



sequential.colors <- c("#1b9e77",
                       "#d95f02",
                       "#7570b3",
                       "#e7298a",
                       "#66a61e",
                       "#e6ab02")

# Functions for plotting  #####

plot_theme <- function() {
  
  theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(), 
          axis.line = element_line(size = line_size), 
          axis.ticks = element_line(size = line_size), 
          axis.text = element_text(color = "black", size = 7), 
          axis.title = element_text(color = "black", size = 7),
          legend.title = element_blank(), 
          legend.background = element_rect(fill = "white"),
          legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
          legend.key = element_rect(fill = "white"),
          legend.position = c(0.85, 0.9)) 
  
  
  
}



# Themes, color, scales etc. #####








