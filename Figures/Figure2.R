############# Figure 2 ##################################


## Panels



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

 # group_by(participant, leg, time, time_pp, group, detrain) %>%
 # summarise(rna = mean(rna, na.rm = TRUE), 
 #           tissue_weight = mean(tissue_weight, na.rm = TRUE)) %>%
  mutate(sample = paste(participant, leg, time, sep = "_")) %>%
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

m1 <- lmer(log(rna) ~ tw + time_pp * group + time_pp:group:detrain + 
             (time_pp|participant) + 
             (1|sample), 
           data = rna_control)
plot(m1)
summary(rePCA(m1))
# Basically no variation in slope terms, removing correlation
m2 <- lmer(log(rna) ~ tw + time_pp * group + time_pp:group:detrain + 
             (1 + dummy(time_pp, "post")||participant) +
             (1|sample), 
           data = rna_control)

summary(rePCA(m2))

# Results in singular fit, no information in the slope... A simpler model 
m3 <- lmer(log(rna) ~ tw + time_pp * group + time_pp:group:detrain + 
             (1|participant) + (1|sample), 
           data = rna_control)


anova(m1, m2, m3)
# Adding the slope term givs additional fit to the model, model 2 is preferred, however model
# 2 is probably to complex borderline singular fit. 

summary(m3)


# Get percentage change in each time-point compared to control using a custom 

# Custom contrast matrix k

head(model.matrix(m3))


con.s1 <- matrix(c(1, 0, 1, 0, 0, 0, 0, 0),1) - matrix(c(1, 0, 0, 0, 0, 0, 0, 0),1)
con.post <- matrix(c(1, 0, 0, 1, 0, 0, 0, 0),1) - matrix(c(1, 0, 0, 0, 0, 0, 0, 0),1)

int.s1       <- matrix(c(1, 0, 1, 0, 1, 1, 0, 0),1) - matrix(c(1, 0, 0, 0, 1, 0, 0, 0),1)
int.post     <- matrix(c(1, 0, 0, 1, 1, 0, 1, 0),1) - matrix(c(1, 0, 0, 0, 1, 0, 0, 0),1) 
int.detrain  <- matrix(c(1, 0, 0, 1, 1, 0, 1, 1),1) - matrix(c(1, 0, 0, 0, 1, 0, 0, 0),1) 

summary(m3)

k <-  rbind(con.s1, 
            con.post,
            int.s1,     
            int.post,   
            int.detrain,
            int.s1 - con.s1, 
            int.post - con.post, 
            int.detrain - con.post) 
# Set rownames
rownames(k) <- c("con_s1", 
                 "con_post", 
                 "int_s1", 
                 "int_post", 
                 "int_detrain", 
                 "inter:int_s1",
                 "inter:int_post",
                 "inter:int_detrain") 

## These confidence intervals are not adjusted. 
ci <- confint(glht(m3, linfct = k), calpha = univariate_calpha())


# Calculate average change scores per participant


ci$confint %>%
  data.frame() %>%
  print()


rna_control %>%

  group_by(participant, leg, time_pp, group, detrain) %>%
  summarise(rna.weight = mean(rna, na.rm = TRUE) / mean(tissue_weight, na.rm = TRUE)) %>%
  
  mutate(time = paste0(time_pp, "_", detrain)) %>%
  ungroup() %>%
  dplyr::select(participant, leg, time,group, rna.weight) %>%

  
  pivot_wider(names_from = time, values_from = rna.weight) %>%

  mutate(S1 = 100 * ((S1_train/pre_train) -1), 
         post = 100 * ((post_train/pre_train) -1), 
         detrain = 100 * ((post_detrain/pre_train) -1)) %>%
  
  dplyr::select(participant, leg, group, S1:detrain) %>%
  pivot_longer(names_to = "time", 
               values_to = "rel_change", 
               cols = S1:detrain) %>%
  filter(!(is.na(rel_change))) %>%
  mutate(time_group = paste0(group, "_", time)) %>%
  

  
  ggplot(aes(time_group, rel_change)) + geom_point()



###### Panel X Intervention vs. control qPCR ##################


qpcr_res_int_con <- readRDS("./data/derivedData/qpcr-analysis-bayes/qpcr_res_int_con.RDS")




interaction_effects <- qpcr_res_int_con %>%
  
  filter(comparison %in% c("inter:S1", "inter:post", "inter:post1w","")) %>%
  
  mutate(comparison = gsub("inter:", "", comparison), 
         comparison = factor(comparison, levels = c("S1", "post", "post1w"))) %>%
  ggplot(aes(target, Estimate)) + 
  geom_point() +
  geom_errorbar(aes(ymin = CI.Lower, ymax = CI.Upper)) + 
  facet_wrap(  comparison ~ .) + coord_flip()


fold_changes <- qpcr_res_int_con %>%
  
  filter(!(comparison %in% c("inter:S1", "inter:post", "inter:post1w",""))) %>%
  separate(comparison, into = c("group", "time"), sep = "_") %>%
  
  mutate(Estimate = exp(Estimate), 
         CI.Lower = exp(CI.Lower), 
         CI.Upper = exp(CI.Upper), 
         time = factor(time, levels = c("S1", "post", "post1w"))) %>%
  ggplot(aes(time, Estimate, fill = group)) + 
  
  geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.15) + 
  
  geom_errorbar(aes(ymin = CI.Lower, 
                    ymax = CI.Upper), 
                position = position_dodge(width = 0.3), 
                width = 0.1) + 
  
  facet_grid(target ~ .)


cowplot::plot_grid(fold_changes, interaction_effects)














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

