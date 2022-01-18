##### Figure source ###############################



# Libraries for figures ####

library(tidyverse)
library(cowplot)
library(ggridges)
library(ggrepel)
library(ggtext)
# Color scale 
# color.scale <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
#                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Updated colors for Acta Physiologica
color.scale <- c("#ffffff", "#4f81bd", "#7f7f7f", "#b42a51")



multiple.color <- "#4f81bd"
single.color   <- "#7f7f7f"

# line sizes
line_size <- 0.3
error.size <- 0.25
# text sizes
label.size <- 10

text.size <- 8

sequential.colors <- c("#1b9e77",
                       "#d95f02",
                       "#7570b3",
                       "#e7298a",
                       "#66a61e",
                       "#e6ab02")


# Set Font to Arial 

windowsFonts()


# Functions for plotting  #####


plot_theme <- function() {
  
  theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_rect(fill = "#fefeda"),
          axis.line = element_line(size = line_size), 
          axis.ticks = element_line(size = line_size), 
          axis.text = element_text(color = "black", size = 7), 
          axis.title = element_text(color = "black", size = 7),
          legend.title = element_blank(), 
          legend.background = element_rect(fill = "white"),
       #   legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
          legend.key = element_rect(fill = "white", color = "white"),
          legend.position = c(0.85, 0.9)) 
  
  
}



# Themes, color, scales etc. #####








