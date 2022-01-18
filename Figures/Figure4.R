##-------------------------------------
## figure4.R
##
## Title: Figure 4
## Purpose: Show relationship between Total RNA increase and muscle hypertrophy
## Author:
##
##
##
##-------------------------------------
## Notes:

# For analysis see ubf-tot-rna-model.R
#
#
#
#
#
#
#
## ------------------------------------


source("./R/figure-source.R")
library(brms)
library(readxl)


#### Figure 4B Leave one out analysis #########################



loo.results <- readRDS("./data/derivedData/ubf-tot-rna-model/leave-one-out.RDS")

leg <- read_excel("./data/leg_randomization.xlsx")


res <- loo.results$loo.sample %>%
  inner_join(leg) %>%
  
  mutate(Participant = factor(participant, 
                              levels = c("Full\nmodel", paste0("P", c(seq(1:7), 19, 21, 22, 23))), 
                              labels = c("Full\nmodel", paste0("P", seq(1:11)))), 
         Participant = fct_rev(Participant)) %>%

  
  pivot_longer(names_to = "stat", values_to = "estimates", estimate:upr) %>%
  print()
  
res2 <- loo.results$loo.participant %>%
  inner_join(leg) %>%
  mutate(Participant = factor(participant, 
                              levels = c("Full\nmodel", paste0("P", c(seq(1:7), 19, 21, 22, 23))), 
                              labels = c("Full\nmodel", paste0("P", seq(1:11)))), 
         Participant = fct_rev(Participant)) %>%
  
  
  filter(model == "m3", !(coef %in% c("sex", "pre"))) %>%
  mutate(coef = factor(coef, levels = c("slope", "intercept"), 
                       labels = c("Total RNA increase\nper session", 
                                  "Average total RNA\n(Session 6)"))) %>%
    print()
  

## Predictions from the model

# get predictions from the model
hyp.pre.m1 <- readRDS("./data/derivedData/ubf-tot-rna-model/hyp_pre_m3.RDS")

m.coefs <- summary(hyp.pre.m1)

full_model_estimates <- data.frame(Participant = "Full\nmodel",
                                   stat = c("estimate", 
                                            "lwr", 
                                            "upr", 
                                            "estimate", 
                                            "lwr", 
                                            "upr"),
                                    coef = c(rep("slope",3), 
                                             rep("intercept", 3)), 
                                    estimates = c(m.coefs$fixed[5, 1], 
                                                  m.coefs$fixed[5, 3], 
                                                  m.coefs$fixed[5, 4], 
                                                  m.coefs$fixed[4, 1], 
                                                  m.coefs$fixed[4, 3], 
                                                  m.coefs$fixed[4, 4])) %>%
  mutate(coef = factor(coef, levels = c("slope", "intercept"),
                       labels = c("Total RNA increase\nper session", 
                                  "Average total RNA\n(Session 6)"))) %>%
  pivot_wider(names_from = stat, 
              values_from = estimates) %>%
  mutate(Participant = factor(Participant, 
                              levels = c("Full\nmodel", paste0("P", c(seq(1:7), 19, 21, 22, 23))), 
                              labels = c("Full\nmodel", paste0("P", seq(1:11)))), 
         Participant = fct_rev(Participant)) %>%
  print()



loo_panel <- res %>%
  filter(stat != "t_val", model == "m3", !(coef %in% c("sex", "pre"))) %>%
  
  mutate(coef = factor(coef, levels = c("slope", "intercept"), 
                       labels = c("Total RNA increase\nper session", 
                                  "Average total RNA\n(Session 6)"))) %>%
  
  ggplot(aes(estimates, Participant, fill = stat)) +

  
  geom_vline(xintercept = 0, lty = 2, color = "gray50") +
    
#  geom_density_ridges(color = "white") +
  
  geom_point(position = position_nudge(y = 0.2),
             size = 1.2,
             shape = 21, 
             alpha = 0.4) +
  
  facet_wrap( ~ coef, scales = "free") +
  
  scale_fill_manual(values = c(color.scale[5], color.scale[4], color.scale[4])) +


  geom_errorbarh(data = res2, aes(fill = NULL, x = estimate, y = Participant, 
                                  xmax =  upr, xmin = lwr), 
                 height = 0) + 
  
  geom_point(data = res2, aes(estimate, Participant), 
             fill = color.scale[3], shape = 22) +
  
  
  # Adds full model estimates
  
  geom_errorbarh(data = full_model_estimates, 
                 aes(fill = NULL, x = estimate, y = Participant, 
                     xmax =  upr, xmin = lwr), 
                 height = 0) + 
  
  
  geom_point(data = full_model_estimates, 
             aes(estimate, Participant, fill = NULL), 
             fill = color.scale[3], 
             size = 3,
             shape = 21) +

  
  
 
   labs(x = "Estimated change in muscle thickness (mm)\nper one unit change in predictor", 
        y = "Participant") +
  scale_y_discrete(limits = c("Full\nmodel", paste0("P", seq(1:11)))) +
  
  plot_theme() + 
  theme(strip.background = element_rect(color = "white", fill = "white"),
        
        strip.text = element_text(size = 7, hjust = 0),
        axis.title.y = element_text(size = 7), 
        legend.position = "none")


### Panel of predicted values against observed 



# Read data frame of predicted values
predict_df <- readRDS("./data/derivedData/ubf-tot-rna-model/predict_df.RDS")

# get predictions from the model
hyp.pre.m1 <- readRDS("./data/derivedData/ubf-tot-rna-model/hyp_pre_m3.RDS")
summary(hyp.pre.m1)
effects_m1 <- conditional_effects(hyp.pre.m1)

mean_diff <- (effects_m1$sex[2, 9] - effects_m1$sex[1, 9]) / 2


prediction_slope <- predict_df %>%
  
  inner_join(leg) %>%
  mutate(cond = factor(cond, levels = c("const", "var"), 
                       labels = c("CONST", 
                                  "VAR")))  %>%

  group_by(participant, leg, cond) %>%
  summarise(slope = mean(slope), 
            mm_incr = mean(mm_incr)) %>%
  ungroup() %>%
   mutate(Participant = factor(participant, 
                              levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                              labels = paste0("P", seq(1:11))), 
         Participant = fct_rev(Participant)) %>%
  ggplot(aes(slope, mm_incr, fill = cond))  + 
  geom_ribbon(data =   effects_m1$slope, 
              aes(ymin = lower__ + mean_diff, 
                  ymax = upper__ + mean_diff), 
              fill = "gray85") +
  geom_line(data =   effects_m1$slope, 
            aes(slope, estimate__ + mean_diff, fill = NULL)) +
  geom_point(shape = 21, size = 2) + 
  geom_text_repel(aes(label = Participant), 
                  size = 2.5) +
  
  scale_y_continuous(limits = c(-3, 5), 
                     breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5), 
                     labels = c("", -2, "", 0, "", 2, "", 4, "")) +
  
  scale_fill_manual(values = c(color.scale[2], color.scale[4])) +
  labs(x = "Total RNA increase per session (%)", 
       y = "Muscle thickness change (mm)") +
  plot_theme() + 
  theme(axis.title.y = element_text(size = 7), 
        legend.margin = margin(t = 1, r = 1, b = 2, l = 1, unit = "pt"),
        legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 7),
        legend.key.size = unit(0.3, "cm"), 
        legend.position = c(0.15, 0.90))



prediction_intercept <- predict_df %>%
  
  inner_join(leg) %>%
  mutate(cond = factor(cond, levels = c("const", "var"), 
                       labels = c("Constant volume", 
                                  "Variable volume")))  %>%
  
  
  group_by(participant, leg, cond) %>%
  summarise(intercept = mean(intercept), 
            mm_incr = mean(mm_incr)) %>%
  ungroup() %>%
  mutate(Participant = factor(participant, 
                              levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                              labels = paste0("P", seq(1:11))), 
         Participant = fct_rev(Participant)) %>%
  ggplot(aes(intercept, mm_incr, fill = cond))  + 
  geom_ribbon(data =   effects_m1$intercept, 
              aes(ymin = lower__ + mean_diff, 
                  ymax = upper__ + mean_diff), 
              fill = "gray85") +
  geom_line(data =   effects_m1$intercept, 
            aes(intercept, estimate__ + mean_diff, fill = NULL)) +
  geom_point(shape = 21, size = 2) + 
  geom_text_repel(aes(label = Participant), 
                  size = 2.5) +
  scale_fill_manual(values = c(color.scale[2], color.scale[4])) +
  labs(x = "Average Total RNA at Session 6\n(Standard deviations from the mean)", 
       y = "Muscle thickness change (mm)") +
  
  scale_y_continuous(limits = c(-3, 5), 
                     breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5), 
                     labels = c("", -2, "", 0, "", 2, "", 4, "")) +
  plot_theme() + 
  theme(axis.title.y = element_text(size = 7), 
        legend.position = "none")


# Combining data to show estimates of each predictor when leaving one out
combined_df <- readRDS("./data/derivedData/ubf-tot-rna-model/predict_combined_df.RDS")


# A explorative plot of linear fits for each participant (RNA to number of sessions).

rna_to_time_estimates <- combined_df %>%
  filter(!(time %in% c("post1w"))) %>% # Removing the de-training estimate as doeas not represent training induced increase.
  
  mutate(Participant = factor(participant, 
                              levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                              labels = paste0("P", seq(1:11)))) %>%

  inner_join(leg) %>%
  mutate(cond = factor(cond, levels = c("const", "var"), 
                       labels = c("Constant volume", 
                                  "Variable volume")))  %>%
  
  # filter(!(participant == "P21" & time == "S12" & leg == "R")) %>%
  
  mutate(time.c = time.c) %>% # this centers number of sessions and the intercept becomes
  # the estimate total RNA (on log scale) in the mid of the training intervention.
  
  ggplot(aes(time.c, rna.mg, group = paste(participant, cond), 
             color = cond, fill = cond)) + 
  geom_point(shape = 21, size = 1.5) + 
  geom_smooth(method = "lm", se = FALSE, size = 0.75) +
  facet_wrap(~ Participant) + 
  
  scale_x_continuous(limits = c(0, 12), 
                     breaks = c(0, 3, 6, 9, 12), 
                     labels = c(0, "", 6, "", 12)) +
  
  scale_fill_manual(values = c(color.scale[2], color.scale[4])) +
  scale_color_manual(values = c(color.scale[2], color.scale[4])) +
  
  
  labs(y = bquote('Total RNA (mg'~'ng'^-1*')'), 
       x = "Sessions") +
  
  plot_theme() + 
  theme(strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(size = 8, hjust = 0),
        axis.title.y = element_text(size = 7), 
        legend.position = "none")




# Combine all elements of plot 



### Whole page width 
figure4 <- plot_grid( 
    plot_grid(prediction_slope, prediction_intercept, align = "h"), 
    plot_grid(rna_to_time_estimates, loo_panel, ncol = 2),
    ncol = 1, 
    rel_heights = c(1, 1)) +
  

  draw_plot_label(label=c("A",  "B",  "C", "D"),
                  x =   c(0.02, 0.53, 0.02, 0.53), 
                  y =   c(0.98, 0.98, 0.48, 0.48),
                  hjust=.5, vjust=.5, size = label.size)

# Half page width 
# figure4 <- 
#   plot_grid(prediction_slope, 
#             prediction_intercept, 
#             rna_to_time_estimates, 
#             loo_panel, 
#   ncol = 1, 
#   rel_heights = c(1, 1, 1, 1)) +
#   
#   
#   draw_plot_label(label=c("A",  "B",  "C", "D"),
#                   x =   c(0.02, 0.02, 0.02, 0.02), 
#                   y =   c(0.98, 0.73, 0.48, 0.23),
#                   hjust=.5, vjust=.5, size = label.size)
#



saveRDS(figure4, "./Figures/rds/figure4_x.RDS")


# Width of figure = 1x columns 8.9 cm
# height of figure = full page = 23 cm


ggsave("figures/figure4.pdf", plot = figure4, width = 7.7 * 2, height = 23 * 0.75, 
       dpi = 600,
       units = "cm", device=cairo_pdf)


  





