######## ECSS conference presentation figures ####################

# Filename: ./figures/ecss2019_figures.R
# Author: Daniel H
# Notes:

# Purpose: Prepare plots for ecss presentation

# Output: Plots saved with ggsave

######## -----------------------####################

source("./R/libraries.R")


dpi <- 230

#### Theme for presenatation

# Changes plot background and text color

presentation_theme <- function() {
  
    theme(axis.text = element_text(size = 7, color = "black"),
          axis.title =  element_text(size = 7, color = "black"), 
          plot.background = element_rect(color = "white", fill = "white"))
  
}



### Background figure on Total RNA and protein synthesis

millward <- read_excel("./presentationData/litrev.xlsx", sheet = "millward1973") 

millward1973_fig1 <- millward %>%
  mutate(group = if_else(group == "A", "Protein fed", "Protein starved")) %>%
  filter(group != "Protein fed") %>%
  ggplot(aes(RNA, proteinSynthesis, fill = group, color = group)) +
  geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) + 
  geom_point(shape = 21, size = 3, alpha = 0.4) +
  scale_x_continuous(name = TeX("$RNA\\;(\\mu g \\; mg^{-1} \\, protein)$"), 
                     breaks = c(0, 4, 8, 12, 16),
                     minor_breaks = c(2, 6, 10, 14),
                     limits = c(0, 16), 
                     expand = c(0, 0)) +
  scale_y_continuous(name = TeX("$Protein\\,synthesis\\,rate\\;(days^{-1})"), 
                     breaks = c(0, 0.04, 0.08, 0.12, 0.16, 0.2),
                     minor_breaks = c(0.02, 0.06, 0.10, 0.14, 0.18, 0.22),
                     limits = c(0, 0.22), 
                     expand = c(0, 0)) + 
  publr_theme() +
  presentation_theme() +
  theme(legend.position = c(0.2, 0.8), 
        legend.title = element_blank(), 
        legend.background = element_rect(color = "white")) +
  scale_fill_manual(values = c("black")) +
  scale_color_manual(values = c("black")) +
  labs(title = "RNA concentrations determines protein synthesis",
       caption = "Data from Millward, D. J., et al. (1973). Nature 241: 204")

millward1973_fig1

millward1973 <- millward %>%
  mutate(group = if_else(group == "A", "Protein fed", "Protein starved")) %>%
  ggplot(aes(RNA, proteinSynthesis, fill = group, color = group)) +
  geom_smooth(method = "lm", se = FALSE, show.legend = FALSE) + 
  geom_point(shape = 21, size = 3, alpha = 0.4) +
  scale_x_continuous(name = TeX("$RNA\\;(\\mu g \\; mg^{-1} \\, protein)$"), 
                     breaks = c(0, 4, 8, 12, 16),
                     minor_breaks = c(2, 6, 10, 14),
                     limits = c(0, 16), 
                     expand = c(0, 0)) +
  scale_y_continuous(name = TeX("$Protein\\,synthesis\\,rate\\;(days^{-1})"), 
                     breaks = c(0, 0.04, 0.08, 0.12, 0.16, 0.2),
                     minor_breaks = c(0.02, 0.06, 0.10, 0.14, 0.18, 0.22),
                     limits = c(0, 0.22), 
                     expand = c(0, 0)) + 
  publr_theme() +
  presentation_theme() +
  theme(legend.position = c(0.2, 0.8), 
        legend.title = element_blank(), 
        legend.background = element_rect(color = "white")) +
  scale_fill_manual(values = c("gray50", "black")) +
  scale_color_manual(values = c("gray50", "black")) +
  labs(title = "RNA concentrations determines protein synthesis",
       caption = "Data from Millward, D. J., et al. (1973). Nature 241: 204")

millward1973


ggsave("figures/millward1973_fig1.tiff", device = "tiff", 
       plot = millward1973_fig1, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")



ggsave("figures/millward1973.tiff", device = "tiff", 
       plot = millward1973, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")




#### Background figure, total RNA content increase in response to training ####

total_rna <- read_excel("./presentationData/litrev.xlsx", sheet = 1)


total_rna_nonrepeated <- total_rna %>%
  filter(statistics == "summary") %>%
  mutate(study = factor(study, 
                        levels = c("RN1644",
                                   "RN1656",
                                   "RN1755",
                                   "RN1809",
                                   "RN1897",
                                   "RN2265"),
                        labels = c("Figueiredo (2015)",
                                   "Stec (2015)",
                                   "Stec (2016)",
                                   "Brook (2016)",
                                   "Reidy (2017)",
                                   "Hammarström (2019)"))) %>%
  group_by(study, subject, n_sessions) %>%
  summarise(fc = mean(totalRNA_foldchange)) %>%
  mutate(study_group = paste0(study, subject)) %>%
  # mutate(age = if_else(subject_age > 40, "old", "young")) %>%
  ggplot(aes(n_sessions, fc, fill = study)) + 
  
  
  geom_point(shape = 21, alpha = 0.6, size = 3) +
  
  scale_x_continuous(breaks = c(0, 1, 4, 9, 12, 16, 18, 31, 36), 
                     name = "Number of sessions") +
  scale_y_continuous(breaks = c(0.8, 1, 1.2, 1.4, 1.6),
                     limits = c(0.8, 1.6), 
                     expand = c(0, 0), 
                     name = "RNA abundance\nFold-change from pre-training") +
  publr_theme() +
  presentation_theme() +
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.background = element_rect(color = "white")) +
  labs(title = "Effect of the number of Resistance Training sessions\non RNA abundance",
       caption = "")


total_rna_repeated <- total_rna %>%
  filter(statistics == "summary") %>%
  mutate(study = factor(study, 
                        levels = c("RN1644",
                                   "RN1656",
                                   "RN1755",
                                   "RN1809",
                                   "RN1897",
                                   "RN2265"),
                        labels = c("Figueiredo (2015)",
                                   "Stec (2015)",
                                   "Stec (2016)",
                                   "Brook (2016)",
                                   "Reidy (2017)",
                                   "Hammarström (2019)"))) %>%
  group_by(study, subject, n_sessions) %>%
  summarise(fc = mean(totalRNA_foldchange)) %>%
  mutate(study_group = paste0(study, subject)) %>%
  # mutate(age = if_else(subject_age > 40, "old", "young")) %>%
  ggplot(aes(n_sessions, fc, fill = study)) + 
  
  geom_line(aes(group = study_group), color = "gray40") +
  
  geom_point(shape = 21, alpha = 0.6, size = 3) +

  scale_x_continuous(breaks = c(0, 1, 4, 9, 12, 16, 18, 31, 36), 
                     name = "Number of sessions") +
  scale_y_continuous(breaks = c(0.8, 1, 1.2, 1.4, 1.6),
                     limits = c(0.8, 1.6), 
                     expand = c(0, 0), 
                     name = "RNA abundance\nFold-change from pre-training") +
  publr_theme() +
  presentation_theme() +
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.background = element_rect(color = "white")) +
  labs(title = "Effect of the number of Resistance Training sessions\non RNA abundance",
       caption = "")

total_rna_nonrepeated

ggsave("figures/total_rna_background1.tiff", device = "tiff", 
       plot = total_rna_nonrepeated, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")

ggsave("figures/total_rna_background2.tiff", device = "tiff", 
       plot = total_rna_repeated, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")




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

#### Study design ####

design_biopsies <- expand.grid(cond = c("Variable volume", "Constant volume"),
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
  # Gray bars 
  #    annotate("rect",
  #                ymin =c(9.1, 6.1, 8.1),
  #                ymax =c(9.5, 6.5, 8.5),
  #                xmin =c(-3.8, -3.8, -3.8), 
  #                xmax =c(15.5, 15.5, 15.5),
  #            fill = "gray40", 
  #            alpha = 0.2) +
  
  # Muscle biopsies
  annotate("rect", 
           xmin = -4, 
           xmax = 15, 
           ymin = 8.95, 
           ymax = 9.65, 
           color = "white",
           fill = "red",
           alpha = 0.2) +  
  
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
  publr_theme() +
  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank(), 
        legend.background = element_blank(),
        legend.spacing.x = unit(0.25, "lines"), 
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7))

design_us <- expand.grid(cond = c("Variable volume", "Constant volume"),
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
  # Gray bars 
  #    annotate("rect",
  #                ymin =c(9.1, 6.1, 8.1),
  #                ymax =c(9.5, 6.5, 8.5),
  #                xmin =c(-3.8, -3.8, -3.8), 
  #                xmax =c(15.5, 15.5, 15.5),
  #            fill = "gray40", 
  #            alpha = 0.2) +
  
  # Muscle biopsies
  annotate("rect", 
           xmin = -4, 
           xmax = 15, 
           ymin = 8.05, 
           ymax = 8.85, 
           color = "white",
           fill = "red",
           alpha = 0.2) +  
  
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
  publr_theme() +
  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank(), 
        legend.background = element_blank(),
        legend.spacing.x = unit(0.25, "lines"), 
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7))


design_strength <- expand.grid(cond = c("Variable volume", "Constant volume"),
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
  # Gray bars 
  #    annotate("rect",
  #                ymin =c(9.1, 6.1, 8.1),
  #                ymax =c(9.5, 6.5, 8.5),
  #                xmin =c(-3.8, -3.8, -3.8), 
  #                xmax =c(15.5, 15.5, 15.5),
  #            fill = "gray40", 
  #            alpha = 0.2) +
  
  # Muscle biopsies
  annotate("rect", 
           xmin = -4, 
           xmax = 15.5, 
           ymin = 5.85, 
           ymax = 6.65, 
           color = "white",
           fill = "red",
           alpha = 0.2) +  
  
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
  publr_theme() +
  theme(legend.position = "none", 
        legend.direction = "vertical", 
        legend.title = element_blank(), 
        legend.background = element_blank(),
        legend.spacing.x = unit(0.25, "lines"), 
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7))


design_strength
design_biopsies
design_us
### Control group 


ctrl_design_biopsies <- expand.grid(cond = c("Variable volume", "Constant volume"),
                           session = seq(1:12), 
                           sets = rep(6)) %>%
  mutate(sets = if_else(cond == "Variable volume" & session %in% c(5, 6, 7, 8), 3, 
                        if_else(cond == "Variable volume" & session %in% c(9, 10, 11, 12), 9, sets))) %>%
  ggplot(aes(session, sets, fill = cond)) +
  
  # Muscle biopsies
  annotate("rect", 
           xmin = -4, 
           xmax = 14, 
           ymin = 8.95, 
           ymax = 9.50, 
           color = "white",
           fill = "red",
           alpha = 0.2) +
  
  
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
  
  
  # Strength vs biopsy
  
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
  
  
  publr_theme() +
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank())

ctrl_design_us <- expand.grid(cond = c("Variable volume", "Constant volume"),
                                    session = seq(1:12), 
                                    sets = rep(6)) %>%
  mutate(sets = if_else(cond == "Variable volume" & session %in% c(5, 6, 7, 8), 3, 
                        if_else(cond == "Variable volume" & session %in% c(9, 10, 11, 12), 9, sets))) %>%
  ggplot(aes(session, sets, fill = cond)) +
  
  # ultra sound
  annotate("rect", 
           xmin = -4, 
           xmax = 14, 
           ymin = 8.05, 
           ymax = 8.55, 
           color = "white",
           fill = "red",
           alpha = 0.2) +
  
  
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
  
  
  # Strength vs biopsy
  
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
  
  
  publr_theme() +
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank())

ctrl_design_strength <- expand.grid(cond = c("Variable volume", "Constant volume"),
                              session = seq(1:12), 
                              sets = rep(6)) %>%
  mutate(sets = if_else(cond == "Variable volume" & session %in% c(5, 6, 7, 8), 3, 
                        if_else(cond == "Variable volume" & session %in% c(9, 10, 11, 12), 9, sets))) %>%
  ggplot(aes(session, sets, fill = cond)) +
  
  # ultra sound
  annotate("rect", 
           xmin = -4, 
           xmax = 14, 
           ymin = 5.85, 
           ymax = 6.65, 
           color = "white",
           fill = "red",
           alpha = 0.2) +
  
  
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
  
  
  # Strength vs biopsy
  
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
  
  
  publr_theme() +
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank())




design_strength
design_biopsies
design_us

ctrl_design_strength
ctrl_design_biopsies
ctrl_design_us




design_full_biopsies <- plot_grid(design_biopsies + ggtitle("Study design"), 
                         ctrl_design_biopsies, 
                                   ncol = 1, 
                                   rel_heights = c(1, 0.6)) 

design_full_us <- plot_grid(design_us + ggtitle("Study design"), 
                             ctrl_design_us, 
                             ncol = 1, 
                             rel_heights = c(1, 0.6))

design_full_strength <- plot_grid(design_strength + ggtitle("Study design"), 
                             ctrl_design_strength, 
                             ncol = 1, 
                             rel_heights = c(1, 0.6))


ggsave("figures/design_biopsies.tiff", device = "tiff", 
       plot = design_full_biopsies, 
       width = 22.65/1.45, 
       height = 15.1/1.45,
       dpi = dpi,
       units = "cm")
ggsave("figures/design_us.tiff", device = "tiff", 
       plot = design_full_us, 
       width = 22.65/1.45, 
       height = 15.1/1.45,
       dpi = dpi,
       units = "cm")
ggsave("figures/design_strength.tiff", device = "tiff", 
       plot = design_full_strength, 
       width = 22.65/1.45, 
       height = 15.1/1.45,
       dpi = dpi,
       units = "cm")


### Training load and intensity plot ####

source("./R/libraries.R") 


### Exploratory graphs ### 
training_sum <-  read_excel("./data/tr010_training.xlsx", sheet = 1, na = "NA") %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  filter(exercise == "legext") %>%
  mutate(set.load = repetitions * load) %>%
  group_by(participant, leg, cond, session) %>%
  summarise(total.load = sum(set.load)) %>%
  group_by(cond, session) %>%
  summarise(total.load = mean(total.load))



training_load <- read_excel("./data/tr010_training.xlsx", sheet = 1, na = "NA") %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  filter(exercise == "legext") %>%
  mutate(set.load = repetitions * load) %>%
  group_by(participant, leg, cond, session) %>%
  summarise(total.load = sum(set.load)) %>%
  ggplot(aes(session, total.load, group = paste0(participant, leg), color = cond)) + 
  geom_line(alpha = 0.6) + 
  geom_line(data = training_sum, 
            aes(session, total.load, group = cond, color = cond), 
            size = 1.5) +
  scale_color_manual(values = c(const_color, var_color)) +
  geom_text(data = data.frame(session = c(3, 3),
                              total.load = c(5000, 5000), 
                              lab = c("Constant volume","Variable volume"),
                              cond = factor(c("const", "var"), 
                                            levels = c("const", "var"))), 
            aes(label = lab, group = NULL, color = NULL), size = textsize + 0.5) +
  
  scale_y_continuous(limits = c(0, 6000), breaks = c(0, 2000, 4000, 6000), expand = c(0, 0)) +
  scale_x_continuous(limits = c(1, 12), breaks = seq(1:12), expand = c(0.1, 0.1)) +
  ylab(TeX("Training load $\\,$ (kg $\\times$ repetitions)")) +
  xlab("Session") +
  facet_grid(. ~ cond) +
  publr_theme() +
  presentation_theme() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank(), 
        legend.position = "none")

load.stats <- read_excel("./data/tr010_training.xlsx", sheet = 1, na = "NA") %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  filter(exercise == "legext") %>%
  mutate(set.load = repetitions * load, 
         week = if_else(session %in% c(1:4), "W1", if_else(session %in% c(5:8), "W2", "W3"))) %>%
  group_by(participant, cond, session, week) %>%
  summarise(load = mean(load, na.rm = TRUE)) %>%
  mutate(week = factor(week, levels = c("W2", "W1", "W3"))) %>%
  print()

load.m1 <- lme(load ~ week * cond, random = list(participant = ~1), 
               data = load.stats)


average_intensity <- load.stats %>%
  group_by(cond, week) %>%
  summarise(m = mean(load, na.rm = TRUE), 
            s = sd(load, na.rm = TRUE)) %>%
  mutate(week = factor(week, levels = c("W1", "W2", "W3"), labels = c("Week 1", "Week 2", "Week 3"))) %>%
  ggplot(aes(week, m, fill = cond)) +
  geom_errorbar(aes(ymin = m - s, ymax = m + s), position = position_dodge(width = 0.6), width = 0.2) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5, color = "black", size = 0.2) +
  geom_segment(aes(x = 1, xend = 1.95, y = 50, yend = 50), size = line_size) +
  geom_segment(aes(x = 2.05, xend = 3, y = 56, yend = 56), size = line_size) +
  geom_segment(aes(x = 0.8, xend = 1.2, y = 42.5, yend = 42.5), size = line_size) +
  geom_segment(aes(x = 1.8, xend = 2.2, y = 47.5, yend = 47.5), size = line_size) +
  geom_segment(aes(x = 2.8, xend = 3.2, y = 53.2, yend = 53.2), size = line_size) +
  geom_segment(aes(x = 3, xend = 3, y = 56, yend = 53.2), size = line_size) +
  geom_segment(aes(x = 2.05, xend = 2.05, y = 56, yend = 47.5), size = line_size) +
  geom_segment(aes(x = 1.95, xend = 1.95, y = 50, yend = 47.5), size = line_size) +
  geom_segment(aes(x = 1, xend = 1, y = 50, yend = 42.5), size = line_size) +
  annotate(geom = "text", x = 1.5, y = 52, 
           label = pval(coef(summary(load.m1))[2, 5], plot = TRUE), 
           parse = TRUE, size = textsize) +
  annotate(geom = "text", x = 2.5, y = 58, 
           label = pval(coef(summary(load.m1))[3, 5], plot = TRUE), 
           parse = TRUE, size = textsize) +
 
   
  
  publr_theme() +
  scale_y_continuous(breaks = c(0, 20, 40, 60), expand = c(0, 0), limits = c(0, 60)) + 
  scale_fill_manual(values = c(const_color, var_color)) +
  ylab("Average 10RM (kg)") +
  theme(axis.title.x = element_blank(), 
        legend.position = "none") +
  presentation_theme() +
  labs(caption = "Mean + SD")

average_intensity

training_load_full <- plot_grid(training_load + ggtitle("Training load and intensity"), average_intensity, ncol = 1)


ggsave("figures/training_load.tiff", device = "tiff", 
       plot = training_load_full, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")

#### Training outcomes ##### 


### Muscle strength 


# Isokinetic strength

strength_absolute <- read_excel("./data/tr010_humac.xlsx", sheet = 1, na = "NA") %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  
  dplyr::select(participant, sex, leg, timepoint, cond, isokinetic_torque, isometric_torque) %>%
  gather(modality, torque, isokinetic_torque:isometric_torque) %>%
  mutate(modality = gsub("_torque", "", modality),
         timepoint = if_else(timepoint %in% c("post_ctrl","postctrl"), 
                             "post", timepoint)) %>%
  group_by(participant, sex, leg, timepoint, cond, modality) %>%
  summarise(torque = mean(torque, na.rm = TRUE)) %>%
  spread(timepoint, torque) %>%
  rowwise() %>%
  mutate(baseline = mean(c(B1, B2), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(participant:modality, B1:baseline) %>%
  mutate(post_1w = if_else(cond == "ctrl_leg", post, post_1w)) %>%
  gather(time, torque, post:baseline) %>%
  mutate(cond = factor(cond, levels = c("ctrl_leg", "const", "var")),
         grouping = paste(cond, time, sep = "_"), 
         test = paste(participant, leg, time)) %>%
  filter(modality == "isokinetic") %>%
  data.frame() %>%
  print()

# Notes on statistical tests:

# For comparison to control, post-values are used for comparison to both post time-points in experimental
# groups. The model is primarily used for plotting, inferense confirms data from ancova model.


isok_abs_post <- lme(torque ~ cond * time, random = list(participant = ~1, leg = ~ 1),
                     data = strength_absolute[strength_absolute$time != "post_1w", ], na.action = na.omit)  

isok_abs_post1w <- lme(torque ~ cond * time, random = list(participant = ~1, leg = ~ 1),
                       data = strength_absolute[strength_absolute$time != "post", ], na.action = na.omit) 


coef(summary(isok_abs_post))
coef(summary(isok_abs_post1w))


anova(isok_abs_post)


pos <- position_dodge(width = 0.2)

isokinetic_time_course <-  data.frame(emmeans(isok_abs_post, ~time|cond)) %>%
  filter(time != "baseline") %>%
  rbind(data.frame(emmeans(isok_abs_post1w, ~time|cond))) %>%
  filter(!(cond == "ctrl_leg" & time == "post")) %>%
  mutate(time = factor(time, levels = c("baseline", "post", "post_1w"), 
                       labels = c("Baseline", "Post\ntraining", "Post\nde-training"))) %>%
  ggplot(aes(time, emmean, fill = cond, group = cond)) + 
  geom_line(position = pos) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, position = pos) + 
  geom_point(shape = 21, size = 3, position = pos) + 
  scale_fill_manual(values = c(control_color, const_color, var_color)) +
  scale_y_continuous(name = "Mean isokinetic torque (Nm \U00B1 95% CI)", lim = c(100, 200), expand = c(0, 0), 
                     breaks = seq(from = 100, to = 200, by = 25)) +
  annotate(geom = "text", x = 2.05, y = 196, label = "\U2020", size = 4) +
  annotate(geom = "text", x = 3.07, y = 193, label = "\U2020", size = 4) +
  publr_theme() +
  presentation_theme() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  labs(title = "Isokinetic strength",  caption = "\U2020, increase from Baseline p < 0.05 vs. control")


isokinetic_time_course


# Isometric strength 

strength_absolute <- read_excel("./data/tr010_humac.xlsx", sheet = 1, na = "NA") %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  
  dplyr::select(participant, sex, leg, timepoint, cond, isokinetic_torque, isometric_torque) %>%
  gather(modality, torque, isokinetic_torque:isometric_torque) %>%
  mutate(modality = gsub("_torque", "", modality),
         timepoint = if_else(timepoint %in% c("post_ctrl","postctrl"), 
                             "post", timepoint)) %>%
  group_by(participant, sex, leg, timepoint, cond, modality) %>%
  summarise(torque = mean(torque, na.rm = TRUE)) %>%
  spread(timepoint, torque) %>%
  rowwise() %>%
  mutate(baseline = mean(c(B1, B2), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(participant:modality, B1:baseline) %>%
  mutate(post_1w = if_else(cond == "ctrl_leg", post, post_1w)) %>%
  gather(time, torque, post:baseline) %>%
  mutate(cond = factor(cond, levels = c("ctrl_leg", "const", "var")),
         grouping = paste(cond, time, sep = "_"), 
         test = paste(participant, leg, time)) %>%
  filter(modality == "isometric") %>%
  data.frame() %>%
  print()

# Notes on statistical tests:

# For comparison to control, post-values are used for comparison to both post time-points in experimental
# groups. The model is primarily used for plotting, inferense confirms data from ancova model.


isom_abs_post <- lme(torque ~ cond * time, random = list(participant = ~1, leg = ~ 1),
                     data = strength_absolute[strength_absolute$time != "post_1w", ], na.action = na.omit)  

isom_abs_post1w <- lme(torque ~ cond * time, random = list(participant = ~1, leg = ~ 1),
                       data = strength_absolute[strength_absolute$time != "post", ], na.action = na.omit) 

summary(isom_abs_post)
summary(isom_abs_post1w)



pos <- position_dodge(width = 0.2)

isometric_time_course <- data.frame(emmeans(isom_abs_post, ~time|cond)) %>%
  filter(time != "baseline") %>%
  rbind(data.frame(emmeans(isom_abs_post1w, ~time|cond))) %>%
  filter(!(cond == "ctrl_leg" & time == "post")) %>%
  mutate(time = factor(time, levels = c("baseline", "post", "post_1w"), 
                       labels = c("Baseline", "Post\ntraining", "Post\nde-training")),
         cond = factor(cond, levels = c("ctrl_leg", "const", "var"), 
                       labels = c("Control", "Constant", "Variable"))) %>%
  ggplot(aes(time, emmean, fill = cond, group = cond)) + 
  geom_line(position = pos) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, position = pos) + 
  geom_point(shape = 21, size = 3, position = pos) + 
  scale_fill_manual(values = c(control_color, const_color, var_color)) +
  scale_y_continuous(name = "Mean isometric torque (Nm \U00B1 95% CI)", 
                     lim = c(120, 280), 
                     expand = c(0, 0), 
                     breaks = seq(from = 120, to = 280, by = 40)) +

  publr_theme() +
  presentation_theme() +
  theme(legend.position = "left", 
        legend.direction = "vertical", 
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.background = element_rect(color = "white")) +
  labs(title = "Isometric strength", caption = " ")


strength <- plot_grid(isometric_time_course, isokinetic_time_course, ncol = 2, 
                      rel_widths = c(1, 0.7))




ggsave("figures/strength.tiff", device = "tiff", 
       plot = strength, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")

### Muscle size

us_long <- read_excel("./data/ultrasound/ultrasound_data.xlsx") %>%
  inner_join(read_csv("./data/ultrasound/ultrasound_codekey.csv")) %>%
  mutate(leg = gsub("VL", "", leg)) %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  dplyr::select(participant, time, leg, sex, cond, code, length) %>%
  mutate(group = if_else(participant %in% paste("P", 1:7, sep = ""), "experiment", "control")) %>%
  group_by(participant, time, leg, sex, cond, group) %>%
  summarise(thickness = mean(length, na.rm = TRUE)) %>%
  spread(time, thickness) %>%
  mutate(post1w = if_else(cond == "ctrl_leg", post, post1w)) %>%
  gather(time, thickness, post:pre) %>%
  ungroup() %>%
  mutate(cond = factor(cond, levels = c("ctrl_leg", "const", "var")),
         time = factor(time, levels = c("pre", "post", "post1w"))) %>%
  print()



us_abs_post <- lme(thickness ~ cond * time, random = list(participant = ~1, leg = ~ 1),
                   data = us_long[us_long$time != "post1w", ], na.action = na.omit)  

us_abs_post1w <- lme(thickness ~ cond * time, random = list(participant = ~1, leg = ~ 1),
                     data = us_long[us_long$time != "post", ], na.action = na.omit) 

summary(us_abs_post)

summary(us_abs_post1w)


us_time_course <- data.frame(emmeans(us_abs_post, ~time|cond)) %>%
  filter(time != "pre") %>%
  rbind(data.frame(emmeans(us_abs_post1w, ~time|cond))) %>%
  filter(!(cond == "ctrl_leg" & time == "post")) %>%
  mutate(time = factor(time, levels = c("pre", "post", "post1w"), 
                       labels = c("Baseline", "Post\ntraining", "Post\nde-training")),
         cond = factor(cond, levels = c("ctrl_leg", "const", "var"), 
                       labels = c("Control", "Constant", "Variable"))) %>%
  ggplot(aes(time, emmean, fill = cond, group = cond)) + 
  geom_line(position = pos) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0, position = pos) + 
  geom_point(shape = 21, size = 3, position = pos) + 
  scale_fill_manual(values = c(control_color, const_color, var_color)) +
  scale_y_continuous(name = "Vastus lateralis thickness (cm \U00B1 95% CI)", 
                     lim = c(15, 30), 
                     expand = c(0, 0), 
                     breaks = seq(from = 15, to = 30, by = 3)) +
  annotate(geom = "text", x = 2.05, y = 27, label = "\U2020", size = 4) +
  annotate(geom = "text", x = 1.95, y = 27.5, label = "\U2020", size = 4) +
  annotate(geom = "text", x = 3.06, y = 28, label = "\U2020", size = 4) +
  annotate(geom = "text", x = 3, y = 28, label = "\U2020", size = 4) +
  publr_theme() +
  presentation_theme() +
  theme(legend.direction = "vertical", 
        legend.title = element_blank(),
        axis.title.x = element_blank(), 
        legend.background = element_rect(color = "white")) +
  labs(title = "Muscle thickness (ultra sound)", caption = "\U2020, increase from Baseline p < 0.05 vs. control")


us_time_course <- plot_grid(us_time_course, rel_widths = c(1), nrow = 1)


ggsave("figures/ultra_sound.tiff", device = "tiff", 
       plot = us_time_course, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")


##### Total RNA plots ######


source("./R/libraries.R")

sample_setup <- read_excel("./data/tr010_mRNASamples.xlsx", na = "NA") %>%
  inner_join(read_excel("./data/tr010_tissue.xlsx", na = "NA") %>%
               filter(sample == "mRNA") %>%
               mutate(date = as.Date(as.numeric(date), origin = "1899-12-30")) %>%
               dplyr::select(participant, leg, time, samplenr, date)) %>%
  filter(IncludeTotalRNA == "YES") %>%
  dplyr::select(participant, leg, time, ExtractionNR, tissue_weight, elution, date) %>%
  print()



### Read files from Total RNA measurements
files <- list.files("./data/TotalRNA/")

results <- list()

for(i in 1:length(files)) {
  
  results[[i]] <- read_excel(paste0("./data/TotalRNA/", files[i], sep = ""), na = "NA")  
  
}

rna <- bind_rows(results) %>%
  filter(Sample != "Blank1") %>%
  separate(Sample, c("participant", "ExtractionNR")) %>%
  dplyr::select(participant, ExtractionNR, TotalRNA) %>%
  mutate(ExtractionNR = as.numeric(ExtractionNR)) %>%
  inner_join(sample_setup) %>%
  mutate(RNA = TotalRNA * elution,
         RNA.per.tissue = RNA/tissue_weight) %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  mutate(time.c = gsub("S", "", time), 
         detrain = if_else(time.c == "post1w", "detrain", "train"),
         time.c = if_else(time.c == "post1w", "12", time.c),
         time.c = if_else(time.c == "postctrl", "12", time.c),
         time.c = gsub("c", "", time.c), 
         time.c = as.numeric(time.c), 
         time = factor(time, levels = c("S0", "S1","S1c", "S4", "S5", "S8", "S9", "S12", "post1w", "postctrl"))) 


#### Remove outliers based on weight to rna relationship ###
# model for RNA concentration to weight relationship
l.mod <- lm(log(TotalRNA) ~ log(tissue_weight), data = rna) 

# Calculates standardized residuals
rna$resid<- resid(l.mod)/sd(resid(l.mod))
# store predicted values for diagnostic plotting
rna$pred<- predict(l.mod)

# plot predicted vs. standardized residuals
rna %>%
  mutate(outlier = if_else(resid < -2.5| resid > 2.5, "out", "in"),
         labels = paste(participant, time, cond)) %>%
  # filter(resid > 3) %>%
  ggplot(aes(log(tissue_weight), log(TotalRNA), color = outlier, label = labels)) + geom_point() 


# First question, does training (with variable and constant volume) 
# lead to RNA accumulation compared to control?

# Questions in modelling 

# 1. Compare a single session to control, is there evidence for rapid increase?
# 2. Compare accumulated training and brief de-training to no-training control
# 3. Estimate slope of increase per session
# 4. Does variable volume affect the slope?




# Data-set used in experiment vs- control comparisons

# Does a single session (of 6 sets) increase total-RNA at 48 h?

unique(rna$time)


dat1 <- rna %>%
  mutate(outlier = if_else(resid < -2.5 | resid > 2.5, "out", "in")) %>%
  # ggplot(aes(tissue_weight, TotalRNA, color = outlier)) + geom_point()
  filter(outlier == "in",
         time %in% c("S0", "S1", "S1c", "S12", "postctrl")) %>%
  mutate(sample = paste0(participant, leg, time),
         leg = paste0(leg),
         cond = factor(cond),
         timepoint = factor(if_else(as.character(time) == "S1c", "S1",
                                    if_else(as.character(time) == "postctrl", "S12", as.character(time))),
                            levels = c("S0", "S1", "S12"))) %>%
  print()



### Explorative graphs ###

# Between leg variations
dat1 %>%
  filter(timepoint == "S0") %>%
  group_by(participant, timepoint, leg) %>%
  summarise(RNA.per.tissue = mean(RNA.per.tissue)) %>%
  spread(leg, RNA.per.tissue) %>%
  ggplot(aes(L, R)) + geom_point()


dat1 %>%
  filter(timepoint == "S0") %>%
  group_by(participant, timepoint, leg) %>%
  summarise(RNA.per.tissue = mean(RNA.per.tissue)) %>%
  ggplot(aes(leg, RNA.per.tissue, group = participant)) + geom_line()



dat1 %>%
  group_by(participant, cond, leg) %>%
  mutate(n = n()) %>%
  filter(n > 4) %>%
  group_by(participant, timepoint, cond, leg) %>%
  summarise(RNA.per.tissue = mean(RNA.per.tissue)) %>%
  ungroup() %>%
  mutate(cond = factor(cond, levels = c("ctrl_leg", "const", "var")),
         group = factor(if_else(cond == "ctrl_leg", "control", "experimental"))) %>%
  ggplot(aes(timepoint, RNA.per.tissue, color = group)) + geom_boxplot() 






# Preliminary model for the control vs. experimental 

m1 <- dat1 %>%
  group_by(participant, cond, leg) %>%
  mutate(n = n()) %>%
  filter(n > 4) %>% # This removes non-repeated legs in the control group
  ungroup() %>%
  mutate(cond = factor(cond, levels = c("ctrl_leg", "const", "var")),
         group = factor(if_else(cond == "ctrl_leg", "control", "experimental"))) %>%
  group_by(participant, group, timepoint, leg) %>%
  summarise(RNA.per.tissue = log(mean(RNA.per.tissue))) %>%
  
  lme(RNA.per.tissue ~ timepoint * group, 
      random = list(participant = ~1, 
                    leg = ~1),
      data = ., 
      na.action = na.omit,
      control=list(msMaxIter = 120,
                   opt = "nloptwrap", msVerbose = FALSE), method = "ML")




plot(m1)



# Leg is kept as a random effect
summary(m1)


pos <- position_dodge(width = 0.2)

total_rna_exp_vs_cont<- data.frame(emmeans(m1, specs = ~ "timepoint|group")) %>%
  mutate(group = factor(group, levels = c("control", "experimental"), labels = c("Control group", "Experimental group\n(Variable + Constant)")),
         timepoint = factor(timepoint, levels = c("S0", "S1", "S12"), labels = c("Pre-\nexercise", 
                                                                                 "48-h Post\nSession 1",
                                                                                 "48-h post\nSession 12"))) %>%
  
  ggplot(aes(timepoint, exp(emmean), fill = group)) + 
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                position = pos, width = 0) +
  geom_point(position = pos, shape = 21, size = 3) +
  scale_y_continuous(breaks = seq(from = 250, to = 625, by = 75),
                     limits = c(250, 625), 
                     expand = c(0,0), 
                     name =   expression("Total RNA per tissue weight (ng mg"^"-1"~" \U00B1 95% CI)")) +
  scale_fill_manual(values = c(control_color, experiment_color)) + 
  annotate(geom = "text", x = 3.05, y = 618, label = "\U2020") +
  publr_theme() +
  presentation_theme() +
  theme(legend.title = element_blank(),
        legend.direction = "vertical",
        legend.background = element_blank(),
        legend.position = c(0.25, 0.9),
        axis.title.x = element_blank()) +
    labs(title = "Total RNA",  caption = "\U2020, increase from Baseline p < 0.05 vs. control")
  



total_rna_exp_vs_cont <- plot_grid(NULL, total_rna_exp_vs_cont, NULL, rel_widths = c(0.2, 1, 0.2), nrow = 1)


ggsave("figures/total_rna_exp_vs_cont.tiff", device = "tiff", 
       plot = total_rna_exp_vs_cont, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")

total_rna_exp_vs_cont


#### Total RNA time-course ### 



##### From ECSS abstract RNA ######
  
  sample_setup <- read_excel("./data/tr010_mRNASamples.xlsx", na = "NA") %>%
    filter(IncludeTotalRNA == "YES") %>%
    dplyr::select(participant, leg, time, ExtractionNR, tissue_weight, elution) 
  
  ### Read files from Total RNA measurements
  files <- list.files("./data/TotalRNA/")
  
  results <- list()
  
  for(i in 1:length(files)) {
    
    results[[i]] <- read_excel(paste0("./data/TotalRNA/", files[i], sep = ""), na = "NA")  
    
  }
  
  rna <- bind_rows(results) %>%
    filter(Sample != "Blank1") %>%
    separate(Sample, c("participant", "ExtractionNR")) %>%
    dplyr::select(participant, ExtractionNR, TotalRNA) %>%
    mutate(ExtractionNR = as.numeric(ExtractionNR)) %>%
    inner_join(sample_setup) %>%
    mutate(RNA = TotalRNA * elution,
           RNA.per.tissue = RNA/tissue_weight) %>%
    inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
    mutate(time.c = gsub("S", "", time), 
           detrain = if_else(time.c == "post1w", "detrain", "train"),
           time.c = if_else(time.c == "post1w", "12", time.c),
           time.c = if_else(time.c == "postctrl", "12", time.c),
           time.c = gsub("c", "", time.c), 
           time.c = as.numeric(time.c), 
           time = factor(time, levels = c("S0", "S1","S1c", "S4", "S5", "S8", "S9", "S12", "post1w", "postctrl")),
           sets = factor(if_else(cond == "var" & time %in% c("S5", "S8"), "sets3", 
                                 if_else(cond == "var" & time %in% c("S9", "S12"), "sets9", "sets6")),
                         levels = c("sets6", "sets3", "sets9"))) 
  
  
  #### Remove outliers based on weight to rna relationship ###
  # model for RNA concentration to weight relationship
  l.mod <- lm(TotalRNA ~ tissue_weight, data = rna) 
  
  # Calculates standardized residuals
  rna$resid<- resid(l.mod)/sd(resid(l.mod))
  # store predicted values for diagnostic plotting
  rna$pred<- predict(l.mod)
  
  # plot predicted vs. standardized residuals
  #rna %>%
  #  mutate(outlier = if_else(resid < -4 | resid > 4, "out", "in"),
  #         labels = paste(participant, time, cond)) %>%
  #  filter(resid > 4) %>%
  #  ggplot(aes(pred, resid, color = outlier, label = labels)) + geom_point() + geom_label()
  
  # creates a outlier variable and filters data
  mdat1 <- rna%>%
    mutate(outlier = if_else(resid < -2.5 | resid > 2.5, "out", "in")) %>%
    filter(outlier == "in",
           participant %in% c("P1", "P2", "P3", "P4", "P5", "P6", "P7"),
           detrain == "train") %>%
    #  group_by(participant, time, time.c, cond, leg) %>%
    #  summarise(RNA.per.tissue = mean(RNA.per.tissue, na.rm = TRUE)) %>%
    mutate(sample = paste0(participant, leg, time),
           leg = paste0(participant, leg),
           time1 = time.c,
           time2 = if_else(time %in% c("S0", "S1", "S4"), 0, time.c - 4),
           time3 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8)) 
  
 # To model Total-RNA increase, we first use all data with 6 sets per 
 # session to estimate the slopes.
  
  # A data frame for predictions 
  
  df <-  data.frame(time1 = rep(c(0, 1, 4, 5, 8, 9, 12), 2),
                    time2 = rep(c(0, 0, 0, 1, 4, 5, 8), 2),
                    time3 = rep(c(0, 0, 0, 0, 0, 1, 4), 2),
                    time.c= rep(c(0, 1, 4, 5, 8, 9, 12), 2),
                    cond = c(rep("var", 7), rep("const", 7)))  
  
  
  
  mod_linear <- mdat1 %>%
    filter(!(cond == "var" & time.c < 4)) %>%
    lme(log(RNA.per.tissue) ~ time.c,
        random  = list(participant = pdDiag( ~(time1+ time2+ time3)),
                       leg = ~ 1,
                       sample = ~1),
        data = ., 
        control=list(msMaxIter=120,
                     opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  df$predict <-  predict(mod_linear, level = 0, 
                         newdata = df)
  
  # Estimate from model
  estimate_linear <- paste0("Estimated increase ",
                            sprintf('%#.1f', data.frame((exp(intervals(mod_linear)$fixed) -1) * 100)[2, 2],1),
                            "% (95% CI: ",
                            round(data.frame((exp(intervals(mod_linear)$fixed) -1) * 100)[2, 1],1),
                            "-",
                            round(data.frame((exp(intervals(mod_linear)$fixed) -1) * 100)[2, 3],1),
                            "%)")
  
  
linear_trend_rna <- df %>%
    mutate(RNA.per.tissue = exp(predict)) %>%
    ggplot(aes(time.c, RNA.per.tissue)) + 
 
    geom_point(data = mdat1 %>%
                  filter(!(cond == "var" & time.c < 4)) %>%
                  group_by(participant, time.c) %>%
                  summarise(RNA.per.tissue = mean(RNA.per.tissue, na.rm = TRUE)), 
               aes(time.c, RNA.per.tissue), 
               alpha = 0.2, 
               size = 3) +
    geom_smooth(data = mdat1 %>%
                filter(!(cond == "var" & time.c < 4)) %>%
                group_by(participant, time.c) %>%
                summarise(RNA.per.tissue = mean(RNA.per.tissue, na.rm = TRUE)), 
              aes(time.c, RNA.per.tissue, group = participant), method = "lm",
                  se = FALSE, 
              alpha = 0.1, 
              color = "gray80") +
    geom_line(size = 2, color = const_color) +
    scale_x_continuous(breaks = c(0, 1, 4, 5, 8, 9, 12), name = "Session") +
    scale_y_continuous(breaks = seq(from = 150, 
                                    to = 675, 
                                    by = 75), 
                       expand = c(0, 0), 
                       limits = c(150, 675), 
                       name = expression("Total RNA per tissue weight (ng mg"^"-1"~")")) +
    labs(title = "Linear trend of Total RNA increase\n in response to constant volume", 
            subtitle = estimate_linear) +
    publr_theme() +
    presentation_theme()
  
  
linear_trend_rna <- plot_grid(NULL, linear_trend_rna, NULL, rel_widths = c(0.2, 1, 0.2), nrow = 1)


ggsave("figures/linear_trend_rna.tiff", device = "tiff", 
       plot = linear_trend_rna, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")
  

  
  
  ## Estimating the total RNA increase in piecewise regression 
  
  mod_picewise <- mdat1 %>%
    filter(!(cond == "var" & time.c < 4)) %>%
    lme(log(RNA.per.tissue) ~ (time1 + time2 + time3),
        random  = list(participant = pdDiag( ~(time1+ time2+ time3)),
                       leg = ~ 1,
                       sample = ~1),
        data = ., 
        control=list(msMaxIter=120,
                     opt = "nloptwrap",msVerbose=FALSE), method = "ML")

  
  
  df$predict <-  predict(mod_picewise, level = 0, 
                         newdata = df)
  
  
  
  
  estimate_piece1 <- paste0("Estimated increase ",
                            sprintf('%#.1f', data.frame((exp(intervals(mod_picewise, which = "fixed")$fixed) -1) * 100)[2, 2],1),
                            "%\n(95% CI: ",
                            sprintf('%#.1f', data.frame((exp(intervals(mod_picewise, which = "fixed")$fixed) -1) * 100)[2, 1],1),
                            ", ",
                            sprintf('%#.1f', data.frame((exp(intervals(mod_picewise, which = "fixed")$fixed) -1) * 100)[2, 3],1),
                            "%)")
  
  estimate_piece2 <- paste0("Estimated deviation ",
                            sprintf('%#.1f', data.frame((exp(intervals(mod_picewise, which = "fixed")$fixed) -1) * 100)[3, 2],1),
                            "%\n(95% CI: ",
                            sprintf('%#.1f', data.frame((exp(intervals(mod_picewise, which = "fixed")$fixed) -1) * 100)[3, 1],1),
                            ", ",
                            sprintf('%#.1f', data.frame((exp(intervals(mod_picewise, which = "fixed")$fixed) -1) * 100)[3, 3],1),
                            "%)")
  
  
piecewise_trend_rna <- df %>%
    mutate(RNA.per.tissue = exp(predict)) %>%
    ggplot(aes(time.c, RNA.per.tissue)) + 
    
    geom_point(data = mdat1 %>%
                 filter(!(cond == "var" & time.c < 4)) %>%
                 group_by(participant, time.c) %>%
                 summarise(RNA.per.tissue = mean(RNA.per.tissue, na.rm = TRUE)), 
               aes(time.c, RNA.per.tissue), 
               alpha = 0.2, 
               size = 3) +
  
    geom_line(size = 2, color = const_color) +
    scale_x_continuous(breaks = c(0, 1, 4, 5, 8, 9, 12), name = "Session") +
    scale_y_continuous(breaks = seq(from = 150, 
                                    to = 675, 
                                    by = 75), 
                       expand = c(0, 0), 
                       limits = c(150, 675), 
                       name = expression("Total RNA per tissue weight (ng mg"^"-1"~")")) +
    labs(title = "Piecewise Linear trend of Total RNA increase\nin response to constant volume", 
         subtitle = " ") +
    annotate("text", x = 3.2, y = 630, label = estimate_piece1) +
    annotate("text", x = 8.2, y = 250, label = estimate_piece2) +
    publr_theme() +
    presentation_theme()
  
  
piecewise_trend_rna <- plot_grid(NULL, piecewise_trend_rna, NULL, rel_widths = c(0.2, 1, 0.2), nrow = 1)
  
  
  ggsave("figures/piecewise_trend_rna.tiff", device = "tiff", 
         plot = piecewise_trend_rna, 
         width = 22.65/1.8, 
         height = 15.1/1.8,
         dpi = dpi,
         units = "cm")
  
  
  


  # Piecewise time trends in both conditions
  # To test if different loading schemes changes the slope of the RNA accumulation ###
  
  m1 <- lme(log(RNA.per.tissue) ~ (time1 + time2 + time3) * cond,
            random  = list(participant = pdDiag( ~(time1+ time2+ time3)),
                           leg = ~ 1,
                           sample = ~1),
            data = mdat1, 
            control=list(msMaxIter=120,
                         opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  summary(m1)
  
  ### Graphical representation 
  
  ### How to plot the model
  
 df <-  data.frame(time1 = rep(c(0, 1, 4, 5, 8, 9, 12), 2),
                   time2 = rep(c(0, 0, 0, 1, 4, 5, 8), 2),
                   time3 = rep(c(0, 0, 0, 0, 0, 1, 4), 2),
                   time.c= rep(c(0, 1, 4, 5, 8, 9, 12), 2),
             cond = c(rep("var", 7), rep("const", 7)))
  
  
 estimate_deviation <- paste0("Estimated deviation\nto session 8 ~ ",
                           round(data.frame((exp(intervals(m1, which = "fixed")$fixed) -1) * 100)[6, 2] +
                                     data.frame((exp(intervals(m1, which = "fixed")$fixed) -1) * 100)[7, 2], 1),
                           "%")
 
 summary(m1)
 
 
 df$predict <-  predict(m1, level = 0, 
          newdata = df)
 
 variable_volume_rna <- df %>%
   mutate(RNA.per.tissue = exp(predict), 
          cond = factor(cond, levels = c("const", "var"),
                        labels = c("Constant", "Variable"))) %>%
   ggplot(aes(time.c, RNA.per.tissue, color = cond)) + 
   
   geom_point(data = mdat1 %>%
                
                group_by(participant,cond, time.c) %>%
                summarise(RNA.per.tissue = mean(RNA.per.tissue, na.rm = TRUE)) %>%
                ungroup() %>%
                mutate(cond = factor(cond, levels = c("const", "var"),
                                     labels = c("Constant", "Variable"))), 
              aes(time.c, RNA.per.tissue), 
              alpha = 0.2, 
              size = 3) +
   
   geom_line(size = 2) +
   scale_x_continuous(breaks = c(0, 1, 4, 5, 8, 9, 12), name = "Session") +
   scale_y_continuous(breaks = seq(from = 150, 
                                   to = 675, 
                                   by = 75), 
                      expand = c(0, 0), 
                      limits = c(150, 675), 
                      name = expression("Total RNA per tissue weight (ng mg"^"-1"~")")) +
   labs(title = "Piecewise trends of Total RNA increase in\n response different training volumes", 
        subtitle = " ") +
   annotate("text", x = 8, y = 250, label = estimate_deviation) +
   scale_color_manual(values = c(const_color, var_color)) +
   publr_theme() +
   presentation_theme() +
   theme(legend.position = "right",
         legend.background = element_rect(color = "white"),
         legend.title = element_blank())
 
 variable_volume_rna <- plot_grid(NULL, variable_volume_rna, rel_widths = c(0.2, 1), nrow = 1)
 
 
 ggsave("figures/variable_volume_rna.tiff", device = "tiff", 
        plot = variable_volume_rna, 
        width = 22.65/1.8, 
        height = 15.1/1.8,
        dpi = dpi,
        units = "cm")
 
 ### Calculation CI for predictions Note: 
 # This code produces confidence bands on the prediction "conditional on 
 # the estimates of the random effect variances". This does not give 
 # much confidence in the confidence bands. Not advisable, see:
 # https://stackoverflow.com/questions/14358811/extract-prediction-band-from-lme-fit
 
 #create design matrix
 Designmat <- model.matrix(eval(eval(m1$call$fixed)[-2]), df[-ncol(df)])
 
 #compute standard error for predictions
 predvar <- diag(Designmat %*% m1$varFix %*% t(Designmat))
 df$SE <- sqrt(predvar) 
 df$SE2 <- sqrt(predvar+m1$sigma^2)
 
 
 
  
  df %>%
    filter(cond == "const") %>%
    ggplot(aes(time.c, predict, color = cond)) + 
    geom_line() +
    geom_ribbon(aes(ymin = predict-2 * SE2, ymax = predict + 2*SE2,
                    color = NULL,
                    group = cond),
                alpha=0.2, fill="red") 
  #  geom_ribbon(aes(ymin = predict-2 * SE,
  #                  ymax = predict + 2*SE, 
  #                  color = NULL,
  #                  group = cond),
  #              alpha=0.2, fill="blue")
  
  totrna.int <- round(data.frame(exp(intervals(m1)$fixed)), 2)
  
   plot(m1)
  
  mdat2 <- rna%>%
    mutate(outlier = if_else(resid < -4 | resid > 4, "out", "in")) %>%
    # group_by(participant, time, time.c, cond, leg, detrain) %>%
    # summarise(RNA.per.tissue = mean(RNA.per.tissue, na.rm = TRUE)) %>%
    filter(participant %in% c("P1", "P2", "P3", "P4", "P5", "P6", "P7"),
           time.c == 12,
           outlier == "in") %>%
    mutate(detrain = factor(detrain, levels = c("train", "detrain")),
           sample = paste0(participant, cond, detrain),
           leg = paste0(participant, leg)) 
  
  
  m2 <- lme(log(RNA.per.tissue) ~ detrain*cond,
            random  = list(participant = ~1, leg = ~1, sample = ~1),
            data = mdat2, 
            control=list(msMaxIter=120,
                         opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  
  pval(coef(summary(m2))[2,5], plot = TRUE)
  
  detrain_rna <- emmeans(ref_grid(m2), specs = ~"cond|detrain") %>%
    data.frame() %>%
    mutate(cond = factor(cond, levels = c("const", "var"), 
                         labels = c("Constant", "Variable")),
           expression = exp(emmean), 
           uci = exp(upper.CL)  ,
           lci = exp(lower.CL) , 
           detrain = factor(detrain, levels = c("train", "detrain"), 
                            labels = c("48-h\nPost-training",
                                       "8 days\n Post-training"))) %>%
    ggplot(aes(detrain, expression, fill = cond)) +
    geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, 
                  position = position_dodge(width = 0.2)) +
    geom_point(shape = 21,
               size = 3,
               position = position_dodge(width = 0.2)) +
    scale_fill_manual(values = c(const_color, var_color)) +
    scale_y_continuous(breaks = seq(from = 250, to = 750, by = 100), 
                       limits = c(250, 750), 
                       expand = c(0, 0), 
                       name = expression("Total RNA per tissue weight (ng mg"^"-1"~")")) +
    publr_theme() +
    presentation_theme() +
    theme(legend.position = c(0.85, 0.9), 
          legend.background = element_rect(color = "white"), 
          legend.title = element_blank(), 
          axis.title.x = element_blank()) +
    annotate("text", x = 1.95, y = 560, label = "*") +
    labs(title = "Total RNA in response to de-training", 
         caption = paste("*, ", pval(coef(summary(m2))[2,5], plot = FALSE),
                         " compared to 48-h"), 
                         parse = FALSE)
  
rna_detrain  <- plot_grid(NULL, detrain_rna, NULL, ncol = 3, rel_widths = c(0.2, 1, 0.2))
  
ggsave("figures/rna_detrain.tiff", device = "tiff", 
       plot = rna_detrain, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")
  

  
  ##### Western blot data #####
  
  western_data <- read_excel("./data/westernBlot_quant/tr010_western.xlsx", na = "NA") %>%
    filter(participant != "LADDER") %>%
    dplyr::select(participant, ExtractionNR, round, gel, well, mean.gray1:total.protein2) %>%
    inner_join(read_excel("./data/tr010_mRNASamples.xlsx", na = "NA") %>%
                 dplyr::select(participant, leg, time, ExtractionNR)) %>%
    inner_join(read_excel("./data/westernBlot_quant/tr010_westernECL.xlsx", na = "NA") %>%
                 dplyr::select(round, gel, well,target, signal)) %>%
    rowwise() %>%
    mutate(total.protein = mean(c(total.protein1, total.protein2)) - mean(c(mean.gray1, mean.gray2, mean.gray3))) %>%
    ungroup() %>%
    dplyr::select(participant, leg, time, gel, well, target, total.protein, signal) %>%
    group_by(gel, participant) %>%
    mutate(tp.factor = total.protein/max(total.protein)) %>%
    inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
    dplyr::select(participant, leg, cond, sex, time, gel, well, target, tp.factor, signal) %>%
    ungroup() %>%
    mutate(expression = signal / tp.factor) %>%
    mutate(time = factor(time, levels = c("S0", "S1","S1c", "S4", "S5", "S6", 
                                          "S8", "S9", "S12", "post1w", "postctrl"))) %>%
    mutate(time.c = gsub("S", "", time), 
           detrain = if_else(time.c == "post1w", "detrain", "train"),
           time.c = if_else(time.c == "post1w", "12", time.c),
           time.c = if_else(time.c == "postctrl", "12", time.c),
           time.c = gsub("c", "", time.c), 
           time.c = as.numeric(time.c), 
           time = factor(time, levels = c("S0", "S1","S1c", "S4", "S5", "S8", "S9", "S12", "post1w", "postctrl")))
  
  
  
  west_training1 <- western_data %>%
    filter(participant %in% c("P1", "P2", "P3", "P4", "P5", "P6", "P7"),
           detrain == "train") %>%
    filter(!(cond == "var" & time.c %in% c(5,8,9,12))) %>%
    #  group_by(participant, time, time.c, cond, leg) %>%
    #  summarise(RNA.per.tissue = mean(RNA.per.tissue, na.rm = TRUE)) %>%
    mutate(sample = paste0(participant, leg, time),
           leg = paste0(participant, leg),
           time1 = time.c,
           time2 = if_else(time %in% c("S0", "S1", "S4"), 0, time.c - 4),
           time3 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8),
           biopsy = if_else(time %in% c("S1", "S5", "S9"), "biopsy", "no.biopsy")) %>%
    data.frame()
  
  west_training <- western_data %>%
    filter(participant %in% c("P1", "P2", "P3", "P4", "P5", "P6", "P7"),
           detrain == "train") %>%
    #  group_by(participant, time, time.c, cond, leg) %>%
    #  summarise(RNA.per.tissue = mean(RNA.per.tissue, na.rm = TRUE)) %>%
    mutate(sample = paste0(participant, leg, time),
           leg = paste0(participant, leg),
           time1 = time.c,
           time2 = if_else(time %in% c("S0", "S1", "S4"), 0, time.c - 4),
           time3 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8),
           biopsy = if_else(time %in% c("S1", "S5", "S9"), "biopsy", "no.biopsy")) %>%
    data.frame()
  
  
  rpS6.m0 <- lme(log(expression) ~ time1,
                 random  = list(participant = pdDiag( ~time1),
                                sample = ~1),
                 na.action = na.exclude,
                 data = west_training1[west_training1$target == "t-s6", ], 
                 control=list(msMaxIter=120,
                              opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  rpS6.m1 <- lme(log(expression) ~ (time1 + time2 + time3),
                 random  = list(participant = pdDiag( ~time1+ time2+ time3),
                                sample = ~1),
                 na.action = na.exclude,
                 data = west_training1[west_training1$target == "t-s6", ], 
                 control=list(msMaxIter=120,
                              opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  
  
  
  anova(rpS6.m0, rpS6.m1)
  summary(rpS6.m0)
  
  
  df <- df %>%
    mutate(biopsy = if_else(time.c %in% c(1, 5, 9), "biopsy", "no.biopsy")) %>%
    print()
  
  
  summary(rpS6.m1)
  
  df$predict <-  predict(rpS6.m0, level = 0, 
                         newdata = df)
  
  # Estimate from model
  estimate_linear <- paste0("Estimated increase ",
                            sprintf('%#.1f', data.frame((exp(intervals(rpS6.m0)$fixed) -1) * 100)[2, 2],1),
                            "%\n(95% CI: ",
                            round(data.frame((exp(intervals(rpS6.m0)$fixed) -1) * 100)[2, 1],1),
                            "-",
                            round(data.frame((exp(intervals(rpS6.m0)$fixed) -1) * 100)[2, 3],1),
                            "%)")
  
  
  linear_trend_rps6 <- df %>%
    mutate(expression = exp(predict)) %>%
    ggplot(aes(time.c, expression / 10^6)) + 
    
    geom_point(data = west_training[west_training$target == "t-s6", ] %>%
                 filter(!(cond == "var" & time.c < 4)) %>%
                 group_by(participant, time.c) %>%
                 summarise(expression = mean(expression, na.rm = TRUE)), 
               aes(time.c, expression / 10^6), 
               alpha = 0.2, 
               size = 3) +
    geom_smooth(data = west_training[west_training$target == "t-s6", ] %>%
                  filter(!(cond == "var" & time.c < 4)) %>%
                  group_by(participant, time.c) %>%
                  summarise(expression = mean(expression, na.rm = TRUE)), 
                aes(time.c, expression / 10^6, group = participant), method = "lm",
                se = FALSE, 
                alpha = 0.1, 
                color = "gray80") +
    geom_line(size = 2, color = const_color) +
   scale_x_continuous(breaks = c(0, 1, 4, 5, 8, 9, 12), name = "Session") +
    scale_y_continuous(breaks = seq(from = 0, 
                                    to = 100, 
                                    by = 20), 
                       expand = c(0, 0), 
                       limits = c(0, 100), 
                       name = expression("RPS6 abundance (AU)")) +
    labs(title = "RPS6 protein", 
         subtitle = estimate_linear) +
    publr_theme() +
    presentation_theme()
  
  
  linear_trend_rps6
  
  
  rpS6.m2 <- lme(log(expression) ~ (time1 + time2 + time3) * cond,
                 random  = list(participant = pdDiag( ~time1+ time2+ time3),
                                sample = ~1),
                 na.action = na.exclude,
                 data = west_training[west_training$target == "t-s6", ], 
                 control=list(msMaxIter=120,
                              opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  df$predict <-  predict(rpS6.m2, level = 0, 
                         newdata = df)
  
  # Estimate from model
  estimate_linear <- paste0("Estimated increase ",
                            sprintf('%#.1f', data.frame((exp(intervals(rpS6.m1)$fixed) -1) * 100)[2, 2],1),
                            "% (95% CI: ",
                            round(data.frame((exp(intervals(rpS6.m1)$fixed) -1) * 100)[2, 1],1),
                            "-",
                            round(data.frame((exp(intervals(rpS6.m1)$fixed) -1) * 100)[2, 3],1),
                            "%)")
  
  
  linear_trend_rps6_volume <- df %>%
    mutate(expression = exp(predict)) %>%
    ggplot(aes(time.c, expression / 10^6, color = cond)) + 
    
    geom_point(data = west_training[west_training$target == "t-s6", ] %>%
                 filter(!(cond == "var" & time.c < 4)) %>%
                 group_by(participant, time.c, cond) %>%
                 summarise(expression = mean(expression, na.rm = TRUE)), 
               aes(time.c, expression / 10^6), 
               alpha = 0.2, 
               size = 3) +
   # geom_smooth(data = west_training[west_training$target == "t-s6", ] %>%
   #               filter(!(cond == "var" & time.c < 4)) %>%
   #               group_by(participant, time.c, cond) %>%
   #               summarise(expression = mean(expression, na.rm = TRUE)), 
   #             aes(time.c, expression / 10^6, group = paste(participant, cond)), method = "lm",
   #             se = FALSE, 
   #             alpha = 0.1, 
   #             color = "gray80") +
    geom_line(size = 2) +
    scale_color_manual(values = c(const_color, var_color)) +
    scale_x_continuous(breaks = c(0, 1, 4, 5, 8, 9, 12), name = "Session") +
    scale_y_continuous(breaks = seq(from = 0, 
                                    to = 100, 
                                    by = 20), 
                       expand = c(0, 0), 
                       limits = c(0, 100), 
                       name = expression("RPS6 abundance (AU)")) +
    labs(title = "RPS6 protein", 
         subtitle = " ") +
    publr_theme() +
    presentation_theme() +
    theme(legend.position = "none")
  
  linear_trend_rps6_volume

  summary(rpS6.m2)
  
  
  #### UBF #####
  ubf.m0 <- lme(log(expression) ~ (time1),
                random  = list(participant = pdDiag( ~time1),
                               sample = ~1),
                na.action = na.exclude,
                data = west_training1[west_training1$target == "t-UBF", ], 
                control=list(msMaxIter=120,
                             opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  
  ubf.m1 <- lme(log(expression) ~ (time1 + time2 + time3),
                 random  = list(participant = pdDiag( ~time1+ time2+ time3),
                                sample = ~1),
                 na.action = na.exclude,
                 data = west_training1[west_training1$target == "t-UBF", ], 
                 control=list(msMaxIter=120,
                              opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  
  anova(ubf.m0, ubf.m1)
  
  df <- df %>%
    mutate(biopsy = if_else(time.c %in% c(1, 5, 9), "biopsy", "no.biopsy")) %>%
    print()
  
  
  summary(ubf.m0)
  
  df$predict <-  predict(ubf.m0, level = 0, 
                         newdata = df)
  
  # Estimate from model
  estimate_linear <- paste0("Estimated increase ",
                            sprintf('%#.1f', data.frame((exp(intervals(ubf.m0)$fixed) -1) * 100)[2, 2],1),
                            "%\n(95% CI: ",
                            round(data.frame((exp(intervals(ubf.m0)$fixed) -1) * 100)[2, 1],1),
                            "-",
                            round(data.frame((exp(intervals(ubf.m0)$fixed) -1) * 100)[2, 3],1),
                            "%)")
  
  
  linear_trend_ubf <- df %>%
    mutate(expression = exp(predict)) %>%
    ggplot(aes(time.c, expression / 10^6)) + 
    
    geom_point(data = west_training[west_training$target == "t-UBF", ] %>%
                 filter(!(cond == "var" & time.c < 4)) %>%
                 group_by(participant, time.c) %>%
                 summarise(expression = mean(expression, na.rm = TRUE)), 
               aes(time.c, expression / 10^6), 
               alpha = 0.2, 
               size = 3) +
    geom_smooth(data = west_training[west_training$target == "t-UBF", ] %>%
                  filter(!(cond == "var" & time.c < 4)) %>%
                  group_by(participant, time.c) %>%
                  summarise(expression = mean(expression, na.rm = TRUE)), 
                aes(time.c, expression / 10^6, group = participant), method = "lm",
                se = FALSE, 
                alpha = 0.1, 
                color = "gray80") +
    geom_line(size = 2, color = const_color) +
    scale_x_continuous(breaks = c(0, 1, 4, 5, 8, 9, 12), name = "Session") +
    scale_y_continuous(breaks = seq(from = 0, 
                                    to = 100, 
                                    by = 20), 
                       expand = c(0, 0), 
                       limits = c(0, 100), 
                       name = expression("UBF abundance (AU)")) +
    labs(title = "UBF protein", 
         subtitle = estimate_linear) +
    publr_theme() +
    presentation_theme()
  
  
  linear_trend_rps6
  
  
  ubf.m2 <- lme(log(expression) ~ (time1 + time2 + time3) * cond,
                 random  = list(participant = pdDiag( ~time1+ time2+ time3),
                                sample = ~1),
                 na.action = na.exclude,
                 data = west_training[west_training$target == "t-UBF", ], 
                 control=list(msMaxIter=120,
                              opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  df$predict <-  predict(ubf.m2, level = 0, 
                         newdata = df)
  
  
  
  linear_trend_ubf_volume <- df %>%
    mutate(expression = exp(predict)) %>%
    ggplot(aes(time.c, expression / 10^6, color = cond)) + 
    
    geom_point(data = west_training[west_training$target == "t-UBF", ] %>%
                 group_by(participant, time.c, cond) %>%
                 summarise(expression = mean(expression, na.rm = TRUE)), 
               aes(time.c, expression / 10^6), 
               alpha = 0.2, 
               size = 3) +
   # geom_smooth(data = west_training[west_training$target == "t-UBF", ] %>%
   #               group_by(participant, time.c, cond) %>%
   #               summarise(expression = mean(expression, na.rm = TRUE)), 
   #             aes(time.c, expression / 10^6, group = paste(participant, cond)), 
   #                 method = "lm",
   #             se = FALSE, 
   #             alpha = 0.1, 
   #             color = "gray80") +
    geom_line(size = 2) +
    scale_color_manual(values = c(const_color, var_color)) +
    scale_x_continuous(breaks = c(0, 1, 4, 5, 8, 9, 12), name = "Session") +
    scale_y_continuous(breaks = seq(from = 0, 
                                    to = 100, 
                                    by = 20), 
                       expand = c(0, 0), 
                       limits = c(0, 100), 
                       name = expression("UBF abundance (AU)")) +
    labs(title = "UBF protein", 
         subtitle = " ") +
    publr_theme() +
    presentation_theme() +
    theme(legend.position = "none")
  
  linear_trend_ubf_volume
  
  summary(ubf.m2)
  
  
linear_trend_proteins <-  plot_grid(linear_trend_rps6, linear_trend_ubf)
  
piecewise_trend_protein <- plot_grid(linear_trend_rps6_volume, linear_trend_ubf_volume)
  
ggsave("figures/linear_trend_proteins.tiff", device = "tiff", 
       plot = linear_trend_proteins, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")

ggsave("figures/piecewise_trend_proteins.tiff", device = "tiff", 
       plot = piecewise_trend_protein, 
       width = 22.65/1.8, 
       height = 15.1/1.8,
       dpi = dpi,
       units = "cm")


  
  
  
  
  
 
  
  west_detrain <- western_data %>%
    filter(participant %in% c("P1", "P2", "P3", "P4", "P5", "P6", "P7"),
           time.c == 12) %>%
    mutate(detrain = factor(detrain, levels = c("train", "detrain")),
           sample = paste0(participant, cond, detrain),
           leg = paste0(participant, leg)) 
  
  
  
  
  s6d <- lme(log(expression) ~ detrain*cond,
             random  = list(participant = ~1, leg = ~1, sample = ~1),
             data = west_detrain[west_detrain$target == "t-s6",], 
             na.action = na.exclude,
             control=list(msMaxIter=120,
                          opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  summary(s6d)
  
  ubfd <- lme(log(expression) ~ detrain*cond,
              random  = list(participant = ~1, sample = ~1),
              data = west_detrain[west_detrain$target == "t-UBF",], 
              na.action = na.exclude,
              control=list(msMaxIter=120,
                           opt = "nloptwrap",msVerbose=FALSE), method = "ML")
  
  summary(ubfd)
  
  
  detrainS6 <- emmeans(ref_grid(s6d), specs = ~"cond|detrain") %>%
    data.frame() %>%
    mutate(cond = factor(cond, levels = c("const", "var"), 
                         labels = c("Constant", "Variable")),
           expression = exp(emmean) / 10^6, 
           uci = exp(upper.CL) / 10^6 ,
           lci = exp(lower.CL) / 10^6, 
           detrain = factor(detrain, levels = c("train", "detrain"), 
                            labels = c("48-h\nPost-training",
                                       "8 days\n Post-training"))) %>%
    ggplot(aes(detrain, expression, fill = cond)) +
    geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, 
                  position = position_dodge(width = 0.2)) +
    geom_point(shape = 21,
               size = 3,
               position = position_dodge(width = 0.2)) +
    scale_fill_manual(values = c(const_color, var_color)) +
    scale_y_continuous(breaks = seq(from = 20, to = 100, by = 20), 
                       limits = c(20, 100), 
                       expand = c(0, 0), 
                       name = "RPS6 protein (AU)") +
    publr_theme() +
    presentation_theme() +
    theme(legend.position = c(0.6, 0.85), 
          legend.background = element_rect(color = "white"), 
          legend.title = element_blank(), 
          axis.title.x = element_blank()) +
    labs(title = "RPS6 protein", 
         caption = " \n ")
    
  
  detrainUBF <- emmeans(ref_grid(ubfd), specs = ~"cond|detrain") %>%
    data.frame() %>%
    mutate(cond = factor(cond, levels = c("const", "var"), 
                         labels = c("Constant", "Variable")),
           expression = exp(emmean) / 10^6, 
           uci = exp(upper.CL) / 10^6 ,
           lci = exp(lower.CL) / 10^6, 
           detrain = factor(detrain, levels = c("train", "detrain"), 
                            labels = c("48-h\nPost-training",
                                       "8 days\n Post-training"))) %>%
    ggplot(aes(detrain, expression, fill = cond)) +
    geom_errorbar(aes(ymin = lci, ymax = uci), width = 0, 
                  position = position_dodge(width = 0.2)) +
    geom_point(shape = 21,
               size = 3,
               position = position_dodge(width = 0.2)) +
    scale_fill_manual(values = c(const_color, var_color)) +
    scale_y_continuous(breaks = seq(from = 0, to = 100, by = 20), 
                       limits = c(0, 100), 
                       expand = c(0, 0), 
                       name = "UBF protein (AU)") +
    publr_theme() +
    presentation_theme() +
    theme(legend.position = "none", 
          legend.background = element_rect(color = "white"), 
          legend.title = element_blank(), 
          axis.title.x = element_blank()) +
    annotate("text", x = 2.05, y = 65, label = "\U2020") +
    annotate("text", x = 1.95, y = 55, label = "*") +
    labs(title = "UBF protein", 
         caption = "\U2020, change from 48-h p < 0.05 compared to Constant\n *, p < 0.05 compared to 48-h")
  

  
  detrainin_protein <- plot_grid(detrainS6, detrainUBF, ncol = 2)
  
  ggsave("figures/detraining_protein.tiff", device = "tiff", 
         plot = detrainin_protein, 
         width = 22.65/1.8, 
         height = 15.1/1.8,
         dpi = dpi,
         units = "cm")
  
  
  
 
  
  #### Correlation between individual slopes and hypertrophy ###
  
  rna %>%
    filter(detrain == "train") %>%
    filter(time.c %in% c(1, 4, 5, 8, 9, 12)) %>%
    filter(participant %in% paste0(rep("P"), seq(1:7))) %>%
    group_by(participant, leg, time.c) %>%
    summarise(RNA.per.tissue = mean(RNA.per.tissue)) %>%
    ggplot(aes(time.c, RNA.per.tissue, group = paste(participant))) +
    geom_point() + geom_smooth(method = "lm", se = FALSE)
  
  
slopes <-   rna %>%
    filter(detrain == "train") %>%
    filter(time.c %in% c(0, 1, 4, 5, 8, 9, 12)) %>%
    filter(participant %in% paste0(rep("P"), seq(1:7))) %>%
    group_by(participant, leg) %>%
    summarise(slope = coef(lm(RNA.per.tissue ~ time.c))[2]) %>%
    print()
  




combined <- us_long %>%
  filter(participant %in% paste0(rep("P"), seq(1:7))) %>%
  spread(time, thickness) %>%
  mutate(increase_post = post/pre, 
         increase_post1w = post1w/pre) %>%
  gather(time, increase, increase_post:increase_post1w) %>%
  group_by(participant, leg, time) %>%
  summarise(increase = mean(increase)) %>%
  inner_join(slopes) %>%
  
  print()
  
mx <- combined %>%
  lme(increase ~ slope * time, 
          random = list(participant = ~1), 
          data = .)
 

summary(mx)  
plot(mx)


combined %>%
ggplot(aes(slope, increase)) + geom_point() +
  geom_smooth(method = "lm")
  print()
  
  
  
  
