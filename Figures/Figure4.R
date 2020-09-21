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


#### Figure 4B Leave one out analysis #########################



loo.results <- readRDS("./data/derivedData/ubf-tot-rna-model/leave-one-out.RDS")

res <- loo.results$loo.sample %>%
  
  mutate(Participant = factor(participant, 
                              levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                              labels = paste0("P", seq(1:11))), 
         Participant = fct_rev(Participant)) %>%

  
  pivot_longer(names_to = "stat", values_to = "estimates", estimate:upr) %>%
  print()
  
res2 <- loo.results$loo.participant %>%
  
  mutate(Participant = factor(participant, 
                              levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                              labels = paste0("P", seq(1:11))), 
         Participant = fct_rev(Participant)) %>%
  
  
  filter(model == "m3", coef != "sex") %>%
  mutate(coef = factor(coef, levels = c("slope", "intercept"), 
                       labels = c("Total RNA increase\nper session", 
                                  "Average total RNA\n(Session 6)"))) %>%
    print()
  


loo_panel <- res %>%
  filter(stat != "t_val", model == "m3", coef != "sex") %>%
  
  mutate(coef = factor(coef, levels = c("slope", "intercept"), 
                       labels = c("Total RNA increase\nper session", 
                                  "Average total RNA\n(Session 6)"))) %>%
  
  ggplot(aes(estimates, Participant, fill = stat)) +

  
  geom_vline(xintercept = 0, lty = 2, color = "gray50") +
    
  geom_density_ridges(color = "white") +
  facet_wrap( ~ coef, scales = "free") +
  

  
  
  

  geom_point(data = res2, aes(estimate, Participant, fill = NULL)) +
  geom_errorbarh(data = res2, aes(fill = NULL, x = estimate, y = Participant, 
                                  xmax =  upr, xmin = lwr), 
                 height = 0.5) + 
 
   labs(x = "Estimated change in muscle thickness (mm)\nper one unit change in predictor", 
        y = "Participant") +
  
  plot_theme() + 
  theme(strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(size = 8, hjust = 0),
        axis.title.y = element_markdown(size = 7), 
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

  group_by(participant, leg) %>%
  summarise(slope = mean(slope), 
            mm_incr = mean(mm_incr)) %>%
  ungroup() %>%
   mutate(Participant = factor(participant, 
                              levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                              labels = paste0("P", seq(1:11))), 
         Participant = fct_rev(Participant)) %>%
  ggplot(aes(slope, mm_incr))  + 
  geom_ribbon(data =   effects_m1$slope, 
              aes(ymin = lower__ + mean_diff, 
                  ymax = upper__ + mean_diff), 
              fill = "gray75") +
  geom_line(data =   effects_m1$slope, 
            aes(slope, estimate__ + mean_diff)) +
  geom_point(shape = 21, size = 2, fill= "blue") + 
  geom_text_repel(aes(label = Participant), 
                  size = 2.5) +
  
  scale_y_continuous(limits = c(-3, 5), 
                     breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5), 
                     labels = c("", -2, "", 0, "", 2, "", 4, "")) +
  
  
  labs(x = "Total RNA increase per session (%)", 
       y = "Muscle thickness change (mm)") +
  plot_theme() + 
  theme(axis.title.y = element_markdown(size = 7), 
        legend.position = "none")



prediction_intercept <- predict_df %>%
  
  group_by(participant, leg) %>%
  summarise(intercept = mean(intercept), 
            mm_incr = mean(mm_incr)) %>%
  ungroup() %>%
  mutate(Participant = factor(participant, 
                              levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                              labels = paste0("P", seq(1:11))), 
         Participant = fct_rev(Participant)) %>%
  ggplot(aes(intercept, mm_incr))  + 
  geom_ribbon(data =   effects_m1$intercept, 
              aes(ymin = lower__ + mean_diff, 
                  ymax = upper__ + mean_diff), 
              fill = "gray75") +
  geom_line(data =   effects_m1$intercept, 
            aes(intercept, estimate__ + mean_diff)) +
  geom_point(shape = 21, size = 2, fill= "blue") + 
  geom_text_repel(aes(label = Participant), 
                  size = 2.5) +
  labs(x = "Average Total RNA at Session 6\n(Standard deviations from the mean)", 
       y = "Muscle thickness change (mm)") +
  
  scale_y_continuous(limits = c(-3, 5), 
                     breaks = c(-3, -2, -1, 0, 1, 2, 3, 4, 5), 
                     labels = c("", -2, "", 0, "", 2, "", 4, "")) +
  plot_theme() + 
  theme(axis.title.y = element_markdown(size = 7), 
        legend.position = "none")


  




# Combining data to show estimates of each predictor when leavining one out


combined_df <- readRDS("./data/derivedData/ubf-tot-rna-model/combined_df.RDS")


# A exploitative plot of linear fits for each participant (RNA to number of sessions).

rna_to_time_estimates <- combined_df %>%
  filter(!(time %in% c("post1w"))) %>% # Removing the de-training estimate as doeas not represent training induced increase.
  
  mutate(Participant = factor(participant, 
                              levels = paste0("P", c(seq(1:7), 19, 21, 22, 23)), 
                              labels = paste0("P", seq(1:11)))) %>%

  
  
  # filter(!(participant == "P21" & time == "S12" & leg == "R")) %>%
  
  mutate(time.c = time.c) %>% # this centers number of sessions and the intercept becomes
  # the estimate total RNA (on log scale) in the mid of the training intervention.
  
  ggplot(aes(time.c, rna.mg, group = paste(participant, leg), 
             color = leg)) + geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Participant) + 
  
  scale_x_continuous(limits = c(0, 12), 
                     breaks = c(0, 3, 6, 9, 12), 
                     labels = c(0, "", 6, "", 12)) +
  
  
  labs(y = "Total RNA (ng \U00D7 mg <sup>-1</sup>)", 
       x = "Sessions") +
  
  plot_theme() + 
  theme(strip.background = element_rect(color = "white", fill = "white"),
        strip.text = element_text(size = 8, hjust = 0),
        axis.title.y = element_markdown(size = 7), 
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






# Width of figure = 1x columns 8.9 cm
# height of figure = full page = 23 cm


ggsave("figures/figure4.pdf", plot = figure4, width = 8.9 * 2, height = 23 * 0.75, 
       dpi = 600,
       units = "cm", device=cairo_pdf)


  





