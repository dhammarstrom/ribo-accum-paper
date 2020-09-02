############# Figure 1 ##################################



# Note: To select estimates set bayes to TRU for Bayesian inference and FALSE
# for frequentist inference

bayes <- TRUE




## Panels
# A: Study design
# B: Performed volume per condition
# C/D: Strength and body composition changes. 


source("./R/figure-source.R")
source("./R/libs.R")




# Colors for plot



control_color <- "white"

const_color <- "#56B4E9"
var_color <- "#009E73"
experiment_color <- "#4abdcb"



line_size <- 0.2
textsize <- 2.5


# design figure

biopsy_glyph <- "\U2193"
strength_glyph <-  "\U2020"
us_glyph <- "\U2021"

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
           size = 0) +
  
  
  
  scale_x_continuous(breaks = seq(1:12), limits = c(-4,15.5), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = c(0, 3, 6, 9, 12), limits = c(0,12), expand = c(0,0), 
                     labels = c("", "3", "6", "9", "")) +
  scale_fill_manual(values = c(color.scale[2], color.scale[3])) +

  annotate("text", x = c(0, 0, -2), y = c(9.1, 8.2, 6), 
           label = c("Muscle biopsy", "Muscle thickness", "Strength tests"), 
           vjust = 0, hjust = 1, alpha = 0.8, size = textsize) + 
  
  # Biopsies
  annotate("text", x = c(0.5, 1.5, 4.5, 5.5, 8.5, 9.5, 12.5, 14.5), y = rep(9.2, 8), 
           label = rep(biopsy_glyph, 8), vjust = 0) +
  
  annotate("text", x = c(0.5, 12.5, 14.5), y = rep(8.2, 3), 
           label = rep(us_glyph, 3), vjust = 0) +
  # Strength tests
  annotate("text", x = c(-1.5, 13, 15), y = rep(6, 3), 
           label = rep(strength_glyph, 3), vjust = 0) +
  
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
  
  
  annotate("segment", x = 12,   xend = 14.5, y = 10.99, yend = 10.99, size = line_size) +
  annotate("segment", x = 14.5, xend = 14.5, y = 10.99, yend = 10.79, size = line_size) +
  annotate("segment", x = 12,   xend = 12,   y = 10.99, yend = 10.79, size = line_size) +
  annotate("segment", x = 13.5, xend = 13.5, y = 10.99, yend = 11.19, size = line_size) +
  annotate("text", x = 13.5, y = 11.75, label = c("8 days"), size = textsize) +
  
  # biopsies vs strength tests
  annotate("segment", x = 14.5, xend = 15, y = 9.95, yend = 9.95, size = line_size) +
  annotate("segment", x = 14.5, xend = 14.5, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 15, xend = 15, y = 9.95, yend = 9.75, size = line_size) +
  annotate("segment", x = 14.75, xend = 14.75, y = 9.95, yend = 10.15, size = line_size) +
  annotate("text", x = 14.75, y = 10.5, label = c("24 h"), size = textsize) +  
  
  
  xlab("Session") + 
  ylab("Number of sets per session") + 

  
  plot_theme() +
  
  
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
  
  annotate("text", x = c(0, 0, -2), y = c(9.1, 8.2, 6), 
           label = c("Muscle biopsy", "Muscle thickness", "Strength tests"), 
           vjust = 0, hjust = 1, alpha = 0.8, size = textsize) + 
  
  # Biopsies
  annotate("text", x = c(0.5, 1.5, 12.5), y = rep(9, 3), 
           label = rep(biopsy_glyph, 3), vjust = 0) +
  
  annotate("text", x = c(0.5, 12.5), y = rep(8.2, 2), 
           label = rep(us_glyph, 2), vjust = 0) +
  # Strength tests
  annotate("text", x = c(-1.5, 13), y = rep(6, 2), 
           label = rep(strength_glyph, 2), vjust = 0) +
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
  

  plot_theme() +
  
  
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_blank()) 





design_full_biopsies <- plot_grid(design_exp, 
                                  design_ctrl, 
                                  ncol = 1, 
                                  rel_heights = c(1, 0.6)) 


# Panel X. Training data plot ##########################



### Exploratory graphs #
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
  scale_color_manual(values = c(color.scale[2], color.scale[3], color.scale[4])) +
  geom_text(data = data.frame(session = c(1, 1),
                              total.load = c(5000, 5000), 
                              
                              lab = c("Constant\nvolume","Variable\nvolume"),
                              cond = factor(c("const", "var"), 
                                            levels = c("const", "var"))), 
            aes(label = lab, group = NULL, color = NULL), 
            size = textsize + 0.5, 
            hjust = 0) +
  
  scale_y_continuous(limits = c(0, 6000), breaks = c(0, 2000, 4000, 6000), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 12), 
                     breaks = c(0,seq(1:12)),
                     labels = c(0, "", "", 3, "", "", 6, "", "", 9, "", "",12),
                     expand = c(0.1, 0.1)) +
  ylab("Training load (kg \U00D7 repetitions)") +
  xlab("Session") +
  facet_grid(. ~ cond) +
  plot_theme() +
  theme(panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        strip.text.x = element_blank(), 
        legend.position = "none")





# Panel X: Strength tests -- ######################################

# Retrieve estimates from modelling 
if(bayes == TRUE) {
  
comp_strength <-  readRDS("./data/derivedData/strength-bayes/comp_strength.RDS")
  
}

if(bayes == FALSE){
  
  comp_strength <- readRDS("./data/derivedData/strength-freq/comp_strength.RDS")
}




  

strength <- read_excel("./data/tr010_humac.xlsx") %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  mutate(time = if_else(timepoint %in% c("B1", "B2", "fam"), "baseline", timepoint), 
         time = if_else(time == "post_ctrl", "post", time),
         cond = if_else(cond == "ctrl_leg", "ctrl", cond)) %>%
  dplyr::select(participant, sex, time, leg, group,cond, isokinetic_torque, isometric_torque) %>%
  group_by(participant,sex,leg, time, group, cond) %>%
  
  # For plotting purposes, the baseline averaged over all attempts
  
  summarise(isok = mean(isokinetic_torque, na.rm = TRUE), 
            isom = mean(isometric_torque, na.rm = TRUE)) %>%
  pivot_longer(names_to = "type", values_to= "torque", cols = isok:isom) %>%
  pivot_wider(names_from = time, values_from = torque) %>%
  mutate(post = post - baseline, 
         post1w = post1w - baseline) %>%
  pivot_longer(names_to = "time", values_to = "change", cols = post:post1w) %>%
  filter(!(is.na(change))) %>%
  group_by(type) %>%
  mutate(baseline = baseline - mean(baseline, na.rm = TRUE)) %>%
  mutate(time_cond = paste0(time, "_", cond), 
         time_group = paste0(time, "_", group), 
         time_group = factor(time_group, levels = c("post_con", 
                                                    "post_int", 
                                                    "post1w_int"))) %>%
  print()



raw_scores <- strength %>%
  ggplot(aes(time_group, change, fill = type)) + 
  
  geom_hline(yintercept = 0, color = "gray90") +
  
  geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.25), 
             shape = 21, size = 1.2, alpha = 0.4) + 
  
  # Isokinetic strength
  geom_errorbar(data = comp_strength[1:3,], 
                aes(time_group, estimate, ymin = lower.CL, ymax = upper.CL, fill = NULL), 
                position = position_nudge(x = - 0.25), width = 0) +
 
   geom_point(data = comp_strength[1:3,], 
              aes(time_group, estimate, fill = NULL), 
             position = position_nudge(x = -0.25), color = "white", size = 0.5, shape = 15) + 
  # Isometric strength
  geom_errorbar(data = comp_strength[6:8,], 
                aes(time_group, estimate, ymin = lower.CL, ymax = upper.CL,  fill = NULL), 
                position = position_nudge(x =  0.25), width = 0) +
  
  geom_point(data = comp_strength[6:8,], 
             aes(time_group, estimate, fill = NULL), 
             position = position_nudge(x = 0.25), color = "white", size = 0.5, shape = 15) + 
  
  
  
  labs(y = "\U0394 Torque (Nm)") +
  
  plot_theme() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "none")



# Strength compare plot




comp_panel <- comp_strength %>%

  filter(time_group %in% c("inter:post_int", "inter:post1w_int")) %>%
  
  mutate(time_group = gsub("inter:", "", time_group)) %>%
  dplyr::select(estimate, lower.CL, upper.CL, type, time_group) %>%


  rbind(data.frame(estimate = NA, lower.CL = NA, upper.CL = NA, type = c("isom", "isok"),
                   time_group = "post_con")) %>%

  mutate(time_group = factor(time_group,  
                             levels = c("post_con", 
                                        "post_int", 
                                        "post1w_int"), 
                             labels =c("CON", "EXP\nS12", "EXP\nDe-train")), 
         type = factor(type, levels = c("isok", "isom"), 
                       labels = c("Isokinetic",
                                  "Isometric"))) %>%
  
  ggplot(aes(time_group, estimate, fill = type)) + 
  
  geom_hline(yintercept = 0, color = "gray90") +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL, color = type), width = 0.1, 
                position = position_dodge(width = 0.25)) + 
  geom_point(size = 1.5, 
             shape = 21,
             position = position_dodge(width = 0.25)) + 
  
  scale_y_continuous(limits = c(-20, 40), 
                     breaks = c(-20, -10, 0, 10, 20, 30, 40), 
                     labels = c(-20, "",   0,  "", 20,  "", 40), 
                     expand = c(0, 0)) +
  
  labs(y = "\U0394 EXP \U2212 \U0394 CON") +
  
  plot_theme() +
  
  theme(axis.title.x = element_blank(), 
        legend.position = c(0.25, 0.85), 
        legend.spacing.y = unit(-0.5, 'cm'),
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0.3), size = 7))



### Save UL results as graph
strength_results <- plot_grid(raw_scores, comp_panel, nrow = 2, align = "v")




# Panel X: Ultra sound -- muscle thickness ########################


#  Load model predictions 
if(bayes == TRUE) {
  
 comp_us <-  readRDS("./data/derivedData/us-freq-bayes/comp_us_bayes.RDS") %>%
   dplyr::select(time_group, estimate, lower.CL, upper.CL)
  
}

if(bayes == FALSE) {
  
  comp_us <-  readRDS("./data/derivedData/us-freq-bayes/comp_us_freq.RDS") %>%
    mutate(time_group = contrast) %>%
    dplyr::select(time_group, estimate, lower.CL, upper.CL)
  
}


us_data <- rbind(read_excel("./data/ultrasound/ultrasound_data.xlsx") %>%
                   inner_join(read_csv("./data/ultrasound/ultrasound_codekey.csv")) %>%
                   mutate(leg = gsub("VL", "", leg)) %>%
                   inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
                   dplyr::select(participant, time, leg, sex, cond, code, length) %>%
                   mutate(group = if_else(participant %in% paste("P", 1:7, sep = ""), "experiment", "control")) %>%
                   group_by(participant, time, leg, sex, cond, group) %>%
                   summarise(thickness = mean(length, na.rm = TRUE)) %>%
                   ungroup(), 
                 read_excel("./data/ultrasound/ultrasound_data_2019.xlsx") %>%
                   inner_join(read_csv("./data/ultrasound/ultrasound_codekey_2019.csv")) %>%
                   mutate(leg = gsub("VL", "", leg)) %>%
                   inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
                   dplyr::select(participant, time, leg, sex, cond, code, length) %>%
                   mutate(group = if_else(participant %in% paste("P", c(1:7, 19:23), sep = ""), "experiment", "control")) %>%
                   group_by(participant, time, leg, sex, cond, group) %>%
                   summarise(thickness = mean(length, na.rm = TRUE)) %>%
                   ungroup()) %>%
  
  mutate(time_pp = if_else(time == "post1w", "post", time),
         # The de-training period get its own coefficient
         detrain = if_else(time == "post1w" & group == "experiment", "detrain", "train"),
         # The effect of de training will be added to the model --within-- the intervention group 
         detrain = factor(detrain, levels = c("train", "detrain")), 
         time = factor(time, levels = c("pre", "post", "post1w")), 
         time_pp = factor(time_pp, levels = c("pre", "post"))) %>%
  
  print()

# 





 us_temp <- us_data %>%
   dplyr::select(participant:sex, group, thickness) %>%
   pivot_wider(names_from = time, values_from = thickness) %>%
   mutate(post = post - pre, 
          post1w = post1w - pre) %>%
   pivot_longer(names_to = "time", values_to = "change", cols = post:post1w) %>%
   filter(!is.na(change)) %>%
   mutate(group = if_else(group == "experiment", "int", "con"), 
          time_group = paste0(time, "_", group)) %>%
   
   print()
 
 



# 
# 
raw_scores <- us_temp %>%
  ggplot(aes(time_group, change)) + 
  
  geom_hline(yintercept = 0, color = "gray90") +
  
  geom_jitter(width = 0.05, size = 1.2, alpha = 0.4) + 
  geom_errorbar(data = comp_us[1:3,], 
                aes(time_group, estimate, ymin = lower.CL, ymax = upper.CL), 
                position = position_nudge(x = 0.25), width = 0) +
  geom_point(data = comp_us[1:3,], 
             aes(time_group, estimate), 
             position = position_nudge(x = 0.25), color = "white", size = 0.5, shape = 15) + 
  
  labs(y = "\U0394 Muscle thickness (mm)") +
  
  plot_theme() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())

# scale_y_continuous(limits = c(-8, 8), 
#                    breaks = c(-2.5, 0, 2.5, 5), 
#                    labels = c(-2.5,   0,  2.5, 5), 
#                    expand = c(0, 0))



comp_panel <- comp_us %>%
  
  filter(time_group %in% c("inter:post_int", "inter:post1w_int")) %>%
  
  mutate(time_group = gsub("inter:", "", time_group)) %>%
  
  rbind(data.frame(time_group = "post_con", estimate = NA, lower.CL = NA, upper.CL = NA )) %>%
  
  mutate(time_group = factor(time_group,  
                             levels = c("post_con", 
                                        "post_int", 
                                        "post1w_int"), 
                             labels = c("CON", "EXP\nS12", "EXP\nDe-train"))) %>%
  
  ggplot(aes(time_group, estimate)) + 

  
  geom_hline(yintercept = 0, color = "gray90") +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) + 
  geom_point(size = 2, shape = 21, fill = "white") + 
  
  scale_y_continuous(limits = c(-0.5, 2.5), 
                     breaks = c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5), 
                     labels = c("",   0,  "", 1,  "", 2,  ""), 
                     expand = c(0, 0)) +
  
  labs(y = "\U0394 EXP \U2212 \U0394 CON") +
  
  plot_theme() +
  
  theme(axis.title.x = element_blank())


### Save UL results as graph
ul_results <- plot_grid(raw_scores, comp_panel, nrow = 2, align = "v")






############## Combine panels  ###################



figure1 <- plot_grid(design_full_biopsies,
                     plot_grid(training_load, strength_results,  ul_results, ncol = 3, rel_widths = c(0.5, 0.25, 0.25)), 
                     rel_heights = c(2, 1),
                     ncol = 1) +
  
  
  draw_plot_label(label=c("A",  "B",  "C", "D"),
                  x =   c(0.02, 0.02, 0.5, 0.75), 
                  y =   c(0.98, 0.33, 0.33, 0.33),
                  hjust=.5, vjust=.5, size = label.size)






# Width of figure = 1x columns 8.9 cm
# height of figure = full page = 23 cm


ggsave("figures/figure1.pdf", plot = figure1, width = 8.9*2, height = 23 * 0.75, 
       dpi = 600,
       units = "cm", device=cairo_pdf)







