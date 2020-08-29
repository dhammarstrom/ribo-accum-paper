############# Figure 1 ##################################


## Panels
# A: Study design
# B: Performed volume per condition
# C/D: Strength and body composition changes. 


source("./R/figure-source.R")
source("./R/libs.R")


# 45S primers graphics #############################

### rRNA and primers #########


# Sequences can be found : https://www.ncbi.nlm.nih.gov/nuccore/NR_146144.1
# Add primer sequences to plot and annotate mature segments

# Primer locations
primers <- read_excel("./data/tr010_primers.xlsx") %>%
  filter(primer_id %in% c("rRNA47S F1R1", 
                          "rRNA45S F5R5", 
                          "rRNA45SITS F12R12", 
                          "rRNA5.8S F2R2", 
                          "rRNA28S F2R2", 
                          "rRNA18S F2R2")) %>%
  print()



primers_fig <- data.frame(x = c(1, 13500), 
           y = c(1, 2)) %>%
  ggplot(aes(x, y)) +
  
  # Full 45S segment 
  annotate("segment", x = 1, xend = 13351, y = 1.5, yend = 1.5, size = 0.5) +
  # 18S segment
  annotate("segment", x = 3655, xend = 5523, y = 1.5, yend = 1.5, size = 2.5) +
  # 5.8S segment
  annotate("segment", x = 6601, xend = 6757, y = 1.5, yend = 1.5, size = 2.5) +
  # 28S segment
  annotate("segment", x = 7925, xend = 12990, y = 1.5, yend = 1.5, size = 2.5) +
  
  geom_segment(data = primers, aes(y = c(1.51, 1.52, 1.53, 1.54, 1.55, 1.56), 
                                   yend = c(1.51, 1.52, 1.53, 1.54, 1.55, 1.56), 
                                   x= on45S_start, xend = on45S_end), 
               size = 4) +
  
  geom_text_repel(data = primers, aes(y = c(1.51, 1.52, 1.53, 1.54, 1.55, 1.56),
                                      x= on45S_start, label = symbol_ext), 
                  size = 2.2, 
                  nudge_y = 0.01) +
  
  scale_y_continuous(limits = c(1.48, 1.58)) +
  plot_theme() + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank())




# Total RNA analysis and plot #################################################

# load full data set
rna <- readRDS("./data/derivedData/tot-rna/tot-rna.RDS")


rna_complete <- rna %>%
 # filter(outlier == "in") %>% # removes outlier based on RNA to tissue weight relationship
  dplyr::select(participant, series, sample, leg, time, cond, tissue_weight, rna) %>%
  mutate(time.c = gsub("S", "", time), 
         detrain = if_else(time.c == "post1w", "detrain", "train"),
         time.c = if_else(time.c == "post1w", "12", time.c),
         time.c = if_else(time.c == "postctrl", "12", time.c),
         time.c = gsub("c", "", time.c), 
         time.c = as.numeric(time.c), 
         time = factor(time, levels = c("S0", "S1","S1c", 
                                        "S4", "S5", "S8", 
                                        "S9", "S12", "post1w", "postctrl"))) %>%
  print()






# Time course data in experimental group #
time_course <- rna_complete %>%
  filter(cond != "ctrl_leg") %>%
  filter(time != "post1w") %>%
  
  group_by(participant,leg, time,time.c, cond, ) %>%
  summarise(rna = mean(rna, na.rm = TRUE), 
            tissue_weight = mean(tissue_weight, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(participant = factor(participant), 
         cond = factor(cond), 
         rna.tissue = log(rna/tissue_weight), 
         tw = tissue_weight - mean(tissue_weight), 
         time.c.cent = time.c - 4) %>%
  print()


# Splines per group
gam2 <- gamm(rna ~ tw +  s(time.c,  k = 7, by = cond), 
           random = list(participant = ~ 1, leg = ~ 1), 
          data = time_course, method = "ML")




# Control vs. experimental group 



rna_control <- rna_complete %>%
  filter(time %in% c("S0","S1", "S1c", "S12", "postctrl", "post1w")) %>%
  mutate(time_pp = if_else(time == "S0", "pre", 
                        if_else(time %in% c("S1", "S1c"), "S1",  "post")), 
         group = if_else(cond == "ctrl_leg", "con", "int")) %>%
  
 # mutate(leg = if_else(group == "con",  "R", leg)) %>%
  
  filter(!(participant == "P10" & time_pp == "S1")) %>%

  group_by(participant, leg, time, time_pp, group, detrain) %>%
  summarise(rna = mean(rna, na.rm = TRUE), 
            tissue_weight = mean(tissue_weight, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(tw = tissue_weight - mean(tissue_weight), 
         detrain = factor(detrain, levels = c("train", "detrain")), 
         time_pp = factor(time_pp, levels = c("pre", "S1", "post"))) %>%
  print()



rna_control %>%
  filter(group == "con") %>%
  group_by(time, participant) %>%
  summarise(tw = mean(tw), 
            tissue_weight = mean(tissue_weight), 
            rna = mean(rna)) %>%
  
ggplot(aes(tissue_weight, rna, label = paste(participant, time))) + geom_point() + geom_text_repel()
  
  
  ggplot(aes(time, rna/tissue_weight, 
             group = participant, 
             label = paste(participant, time))) + geom_line() + geom_text_repel()




# Testing random effects in the control vs. experimental group model.

m1 <- lmer(rna ~ tw + time_pp * group + time_pp:group:detrain + (time_pp|participant), 
           data = rna_control)

summary(rePCA(m1))
# Basically no variation in slope terms, removing correlation
m2 <- lmer(rna ~ tw + time_pp * group + time_pp:group:detrain + (dummy(time_pp, "post")||participant), 
           data = rna_control)

summary(rePCA(m2))

# Results in singular fit, no information in the slope... A simpler model 
m3 <- lmer(rna ~ tw + time_pp * group + time_pp:group:detrain + (1|participant), 
           data = rna_control)


anova(m1, m2, m3)
# Adding the slope term givs additional fit to the model, model 2 is preferred. 

summary(m2)


control_em <- emmeans(m2, specs = ~ time_pp |  group + detrain) %>%
  data.frame() %>%
  mutate(emmean = emmean / mean(rna_control$tissue_weight), 
         lower.CL = lower.CL / mean(rna_control$tissue_weight), 
         upper.CL = upper.CL / mean(rna_control$tissue_weight), 
         time.c = if_else(time_pp == "pre", 0, 
                          if_else(time_pp == "S1", 1, 
                                  if_else(time_pp == "post", 12, 0))), 
         time.c = if_else(detrain == "detrain", 14, time.c)) %>%
  print()




## Predicted values from 
total_rna_fig <- emmeans(gam2, specs =  ~ cond|time.c, at = list(time.c = c(0, 1,2, 3,4, 5,6,7, 8, 9,10, 11, 12)), 
        data = data.frame(time_course)) %>%
  data.frame() %>%
  mutate(emmean   = emmean   / mean(time_course$tissue_weight), 
         lower.CL = lower.CL / mean(time_course$tissue_weight),
         upper.CL = upper.CL / mean(time_course$tissue_weight)) %>%
  ggplot(aes(time.c, emmean, group = cond, fill = cond, color = cond)) + 
  
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, color = NULL), alpha = 0.2) +
  

  geom_errorbar(data = control_em, 
                aes(time.c, emmean, ymin = lower.CL, ymax = upper.CL, 
                    fill = NULL,
                    color = NULL, 
                    group = group), 
                position = position_dodge(width =  0.2), 
                width = 0.2, 
                color = "gray50") +
  
  geom_point(data = control_em, 
             aes(time.c, emmean, fill = group, color = NULL, group = NULL), 
             shape = 21, 
          
             position = position_dodge(width = 0.2)) +
                
  
  
  geom_line()  +
  
  plot_theme() + 
  theme(legend.position = c(0.2, 0.8))





############## Combine panels  ###################



figure2 <- plot_grid(plot_grid(total_rna_fig,NULL, ncol = 2),
                     plot_grid(primers_fig, NULL,  NULL, ncol = 3, rel_widths = c(0.5, 0.25, 0.25)), 
                     rel_heights = c(1, 1),
                     ncol = 1) +
  
  
  draw_plot_label(label=c("A",  "B",  "C", "D"),
                  x =   c(0.02, 0.02, 0.5, 0.75), 
                  y =   c(0.98, 0.33, 0.33, 0.33),
                  hjust=.5, vjust=.5, size = label.size)






# Width of figure = 1x columns 8.9 cm
# height of figure = full page = 23 cm


ggsave("figures/figure2.pdf", plot = figure2, width = 8.9*2, height = 23 * 0.75, 
       dpi = 600,
       units = "cm", device=cairo_pdf)

