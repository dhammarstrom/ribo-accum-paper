##-------------------------------------
## ubf_tot_rna_model.R
##
## Title: 
## Purpose: Model total RNA over time and influence of UBF levels
## Author:
##
##
##
##-------------------------------------
## Notes:
#
#
#
#
#
#
#
#
## ------------------------------------

# Load packages
library(tidyverse)
library(brms)
library(bayesplot)
library(lme4)





# Load RNA data 

rna <- readRDS("./data/derivedData/tot-rna/tot-rna.RDS")



rna_complete <- rna %>%
  filter(outlier == "in") %>% # removes outlier based on RNA to tissue weight relationship
  dplyr::select(participant, series, sample, leg, time, cond, tissue_weight, rna) %>%
  mutate(time.c = gsub("S", "", time), 
         detrain = if_else(time.c == "post1w", "detrain", "train"),
         time.c = if_else(time.c == "post1w", "12", time.c),
         time.c = if_else(time.c == "postctrl", "12", time.c),
         time.c = gsub("c", "", time.c), 
         time.c = as.numeric(time.c), 
         time = factor(time, levels = c("S0", "S1","S1c", 
                                        "S4", "S5", "S8", 
                                        "S9", "S12", "post1w", "postctrl")), 
         
         # Create segments for segmented model
         
         time1 = time.c,
         time2 = if_else(time %in% c("S0", "S1", "S4"), 0, time.c - 4),
         time3 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8) ) %>%
  print()


# Load western blod data 

western_data <- readRDS("./data/derivedData/western-compile/western_data.RDS")


west_cont <- western_data %>%
  
  mutate(time.c = gsub("S", "", time), 
         detrain = if_else(time.c == "post1w", "detrain", "train"),
         time.c = if_else(time.c == "post1w", "12", time.c),
         time.c = if_else(time.c == "postctrl", "12", time.c),
         time.c = gsub("c", "", time.c), 
         time.c = as.numeric(time.c), 
         time = factor(time, levels = c("S0", "S1","S1c", 
                                        "S4", "S5", "S8", 
                                        "S9", "S12", "post1w", "postctrl")), 
         
         # Create segments for segmented model
         
         time1 = time.c,
         time2 = if_else(time %in% c("S0", "S1", "S4"), 0, time.c - 4),
         time3 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8) ) %>%
  
  group_by(target, participant, time.c, cond, time, time1, time2, time3) %>%
 # Summarise duplicate values (if any) keep scaled expression for modeling
   summarise(expression = mean(expression), 
            expression.scaled = mean(expression.scaled)) %>%

  filter(!is.na(expression)) %>%
  
  print()

west_cont %>%
  filter(participant == "P1", cond == "const") %>%
  dplyr::select(time, time1, time2, time3)

# Combined data from both Total-RNA and western blot analysis


combined_df  <- west_cont %>%
  group_by(target, participant, cond, time) %>%
  summarise(ubf.scaled = mean(expression.scaled, na.rm = TRUE)) %>%
  filter(target != "t-s6") %>%
  ungroup() %>%
  dplyr::select(participant, cond, time, ubf.scaled) %>%
  inner_join(rna_complete %>%
               mutate(rna.mg = rna/tissue_weight) %>%
               group_by(participant, leg, time,time.c, time1, time2, time3, cond) %>%
               summarise(rna.mg = mean(rna.mg, na.rm = TRUE))) %>%
  mutate(detrain = factor(if_else(time == "post1w", "detrain", "train"), 
                          levels = c("train", "detrain")), 
         log.rna = log(rna.mg)) %>%
  filter(cond != "ctrl_leg") %>%
  print()





### Combined df in brms models ##############################


comb.m1 <- brm(bf(log.rna ~ ubf.scaled + (time1 +time2 + time3) + detrain + (1 + time1|participant)),
              
              warmup = 1000, # number of samples before sampling
              iter = 4000,  # number of mcmc iterations 
              cores = 4, # number of cores used in sampling
              chains = 4, # number of chains
              seed = 5, # seed for reproducible results
              #    control = list(adapt_delta = 0.95), 
              data = combined_df)




pp_check(comb.m1, type = "ecdf_overlay")

pp_check(comb.m1, type = "stat", stat = "median")

loo(comb.m1)

summary(comb.m1)


hypothesis(comb.m1, c("time1 = 0", 
                      "time1 + time2 = 0", 
                      "time1 + time2 + time3 = 0"))



conditional_effects(comb.m1)



bayes_R2(comb.m1)



pred1 <- expand_grid(ubf.scaled = 0, 
                      time1 = c(0, 1, 4), 
                      time2 = c(0, 0, 0), 
                      time3 = c(0, 0, 0), 
                      detrain = "train",
                      participant = paste0("P", c(seq(1:7),19, 21, 22, 23)))



pred2 <- data.frame(ubf.scaled = rep(0.5, 8), 
                    time1 = c(0, 1, 4, 5, 8, 9, 12, 12), 
                    time2 = c(0, 0, 0, 1, 4, 5, 8, 8), 
                    time3 = c(0, 0, 0, 0, 0, 1, 4, 4), 
                    detrain = c(rep("train", 7), "detrain" ), 
                    participant = "P23")

temp <- cbind(pred1, predict(comb.m1, newdata = pred1)) %>%
  print()

 
temp %>%
  ggplot(aes(time1, Estimate, group = participant)) + geom_line()




 
cbind(combined_df, predict(comb.m1, re_formula = NA))  %>%
  filter(detrain == "train") %>%
  ggplot(aes(time.c, exp(Estimate), group = participant)) + geom_line() + 
  geom_ribbon(aes(ymin = exp(Q2.5), ymax = exp(Q97.5), color = NULL), 
              alpha = 0.1)






# Save model 
saveRDS(comb.m1, "./data/derivedData/ubf-tot-rna-model/comb_m1.RDS")




 ### 






















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


predict_df <- combined_df %>%
  filter(time != "post1w") %>%
  filter(time.c < 13) %>%

  group_by(participant, leg) %>%
  do(slope = coef(lm(log(rna.mg) ~ time.c, data = .))[2]) %>%
  unnest(slope) %>%
  inner_join(us_data %>%
               filter(cond != "ctrl_leg") %>%
               dplyr::select(-time_pp, -detrain) %>%
               pivot_wider(names_from = time, values_from = thickness) %>%
               mutate(post = post - pre, 
                      post1w = post1w - pre) %>%  
               group_by(participant, leg) %>%
              pivot_longer(names_to = "time", values_to = "mm_incr", cols = post:post1w)) %>%

  mutate(slope = 100 * (exp(slope) - 1)) %>%
  print()
  

predict_df %>%
  lmer(mm_incr  ~ time + slope  + (1|participant/leg), data = .) %>%
  summary()
  
  

combined_df %>%
  inner_join(us_data %>%
               filter(cond != "ctrl_leg") %>%
               dplyr::select(-time_pp, -detrain) %>%
               pivot_wider(names_from = time, values_from = thickness) %>%
               mutate(post = post/pre, 
                      post1w = post1w/pre)) %>%
  ggplot(aes(rna.mg, post1w)) + geom_point() + facet_wrap(~ time) + geom_smooth(method = "lm")
  print()








combined_df %>%

  lmer(log.rna ~  (time1 +time2 + time3) + detrain + (1 |participant), data = .)  %>%
  summary()

anova(mx)

ggpredict(mx)





exp(0.05591)
exp(0.1373)

exp(0.09669)





########################################################################

##### Drop in UBF in detraining period #################################



s6.m1 <- readRDS( "./data/derivedData/western-analysis/s6_tx_m1.RDS")
ubf.m1 <- readRDS( "./data/derivedData/western-analysis/ubf_tx_m1.RDS")





