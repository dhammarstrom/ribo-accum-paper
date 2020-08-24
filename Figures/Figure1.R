############# Figure 1 ##################################


## Panels
# A: Study design
# B: Performed volume per condition
# C/D: Strength and body composition changes. 


source("./R/figure-source.R")





# Colors for plot

control_color <- "white"

const_color <- "#56B4E9"
var_color <- "#009E73"
experiment_color <- "#4abdcb"



line_size <- 0.2
textsize <- 2.5


# design figure

biopsy_glyph <- "\U2193"
strength_glyph <-  "\U2299"
us_glyph <- "\U2295"

#### Design graph: EXP ####

design_exp <- expand.grid(cond = c("Variable volume", "Constant volume"),
                               session = seq(1:12), 
                               sets = rep(6)) %>%
  mutate(sets = if_else(cond == "Variable volume" & session %in% c(5, 6, 7, 8), 3, 
                        if_else(cond == "Variable volume" & session %in% c(9, 10, 11, 12), 9, sets))) %>%
  ggplot(aes(session, sets, fill = cond)) +
  
  annotate("text", x = -2, y = 11.0, label = c("Experimental group"), size = textsize + 0.5) +  
  
  geom_bar(stat = "identity", position = position_dodge(width = 0.55),
           width = 0.45, 
           color = "black", 
           size = 0) +
  scale_x_continuous(breaks = seq(1:12), limits = c(-4,15.5), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0,12), expand = c(0,0), 
                     labels = c("", "3", "6", "9", "")) +
  scale_fill_manual(values = c(const_color, var_color)) +

  annotate("text", x = c(-3.8, -3.8, -3.8), y = c(9.1, 8.2, 6), 
           label = c("Muscle biopsy", "Muscle thickness", "Strength tests"), 
           vjust = 0, hjust = 0, alpha = 0.8, size = textsize) + 
  
  # Biopsies
  annotate("text", x = c(0.5, 1.5, 4.5, 5.5, 8.5, 9.5, 12.5, 14.5), y = rep(9, 8), 
           label = rep(biopsy_glyph, 8), vjust = 0, family = "Cambria") +
  
  annotate("text", x = c(0.5, 12.5, 14.5), y = rep(8.2, 3), 
           label = rep(us_glyph, 3), vjust = 0, family = "Cambria") +
  # Strength tests
  annotate("text", x = c(-1.5, 13, 15), y = rep(6, 3), 
           label = rep(strength_glyph, 3), vjust = 0, family = "Cambria") +
  
  # time-frames indicators
  # Strength test prior to experiments
  annotate("segment", x = -1.5, xend = 0.5, y = 6.95, yend = 6.95, size = line_size) +
  annotate("segment", x = 0.5, xend = 0.5, y = 6.95, yend = 7.05, size = line_size) +
  annotate("segment", x = -1.5, xend = -1.5, y = 6.95, yend = 6.75, size = line_size) +
  annotate("segment", x = -0.5, xend = -0.5, y = 6.95, yend = 7.1, size = line_size) +
  annotate("text", x = -0.5, y = 7.55, label = c("> 7 days"), size = textsize) +
  
  
  # post 1st session biopsies
  annotate("segment", x = 1, xend = 1.5, y = 9.95, yend = 9.95, size = line_size) +
  annotate("segment", x = 1.5, xend = 1.5, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 1, xend = 1, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 1.25, xend = 1.25, y = 9.95, yend = 10.15, size = line_size) +
  annotate("text", x = 1.25, y = 10.5, label = c("48 h"), size = textsize) +
  
  # post last session biopsies
  annotate("segment", x = 12, xend = 12.5, y = 9.95, yend = 9.95, size = line_size) +
  annotate("segment", x = 12.5, xend = 12.5, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 12, xend = 12, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 12.25, xend = 12.25, y = 9.95, yend = 10.15, size = line_size) +
  annotate("text", x = 12.25, y = 10.5, label = c("48 h"), size = textsize) +
  
  
  annotate("segment", x = 12, xend = 14.5, y = 10.95, yend = 10.95, size = line_size) +
  annotate("segment", x = 14.5, xend = 14.5, y = 10.95, yend = 10.75, size = line_size) +
  annotate("segment", x = 12, xend = 12, y = 10.95, yend = 10.75, size = line_size) +
  annotate("segment", x = 13.5, xend = 13.5, y = 10.95, yend = 11.15, size = line_size) +
  annotate("text", x = 13.5, y = 11.75, label = c("8 days"), size = textsize) +
  
  # biopsies vs strength tests
  annotate("segment", x = 14.5, xend = 15, y = 9.95, yend = 9.95, size = line_size) +
  annotate("segment", x = 14.5, xend = 14.5, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 15, xend = 15, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 14.75, xend = 14.75, y = 9.95, yend = 10.15, size = line_size) +
  annotate("text", x = 14.75, y = 10.5, label = c("24 h"), size = textsize) +  
  
  
  xlab("Session") + 
  ylab("Number of sets per session") + 

  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank(), 
        legend.background = element_blank(),
        legend.spacing.x = unit(0.25, "lines"), 
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7))



### Design graph: Control group #####

design_ctrl <- expand.grid(cond = c("Variable volume", "Constant volume"),
                                    session = seq(1:12), 
                                    sets = rep(6)) %>%
  mutate(sets = if_else(cond == "Variable volume" & session %in% c(5, 6, 7, 8), 3, 
                        if_else(cond == "Variable volume" & session %in% c(9, 10, 11, 12), 9, sets))) %>%
  ggplot(aes(session, sets, fill = cond)) +
  
  annotate("text", x = -1.85, y = 11.2, label = c("Control group"), size = textsize + 0.5) +  
  
  annotate("text", x = c(-3.8, -3.8, -3.8), y = c(9.1, 8.2, 6), 
           label = c("Muscle biopsy", "Muscle thickness", "Strength tests"), 
           vjust = 0, hjust = 0, alpha = 0.8, size = textsize) + 
  
  # Biopsies
  annotate("text", x = c(0.5, 1.5, 12.5), y = rep(9, 3), 
           label = rep(biopsy_glyph, 3), vjust = 0, family = "Cambria") +
  
  annotate("text", x = c(0.5, 12.5), y = rep(8.2, 2), 
           label = rep(us_glyph, 2), vjust = 0, family = "Cambria") +
  # Strength tests
  annotate("text", x = c(-1.5, 13), y = rep(6, 2), 
           label = rep(strength_glyph, 2), vjust = 0, family = "Cambria") +
  scale_x_continuous(breaks = seq(1:12), limits = c(-4,15.5), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(5,11.5), expand = c(0,0), 
                     labels = c("", "3", "6", "9", "")) +
  
  
  # Strength vs biops
  annotate("segment", x = -1.5, xend = 0.5, y = 6.85, yend = 6.85, size = line_size) +
  annotate("segment", x = 0.5, xend = 0.5, y = 6.85, yend = 7.05, size = line_size) +
  annotate("segment", x = -1.5, xend = -1.5, y = 6.85, yend = 6.65, size = line_size) +
  annotate("segment", x = -0.5, xend = -0.5, y = 6.85, yend = 7.05, size = line_size) +
  annotate("text", x = -0.5, y = 7.35, label = c("> 7 days"), size = textsize) +
  
  
  # Biopsy 1
  annotate("segment", x = 0.5, xend = 1.5, y = 9.95, yend = 9.95, size = line_size) +
  annotate("segment", x = 1.5, xend = 1.5, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 0.5, xend = 0.5, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 1, xend = 1, y = 9.95, yend = 10.05, size = line_size) +
  annotate("text", x = 1, y = 10.25, label = c("48 h"), size = textsize) +
  
  # Biopsy 2
  annotate("segment", x = 0.5, xend = 12.5, y = 10.75, yend = 10.75, size = line_size) +
  annotate("segment", x = 12.5, xend = 12.5, y = 10.75, yend = 10.55, size = line_size) +
  annotate("segment", x = 0.5, xend = 0.5, y = 10.75, yend = 10.55, size = line_size) +
  annotate("segment", x = 6.25, xend = 6.25, y = 10.75, yend = 10.95, size = line_size) +
  annotate("text", x = 6.25, y = 11.25, label = c("2-4 weeks"), size = textsize) +
  
  # Second strength test
  annotate("segment", x = 12.5, xend = 13, y = 9.95, yend = 9.95, size = line_size) +
  annotate("segment", x = 12.5, xend = 12.5, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 13, xend = 13, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 12.75, xend = 12.75, y = 9.95, yend = 10.05, size = line_size) +
  annotate("text", x = 12.75, y = 10.25, label = c("24 h"), size = textsize) +  
  

  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank())





design_full_biopsies <- plot_grid(design_exp + ggtitle("Study design"), 
                                  design_ctrl, 
                                  ncol = 1, 
                                  rel_heights = c(1, 0.6)) 


# Panel B #############################








