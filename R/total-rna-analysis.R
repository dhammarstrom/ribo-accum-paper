######## Total RNA analysis ########################

# Load libraries and data 

library(tidyverse)
library(tidybayes)
library(lme4)
library(brms)
library(readxl)
library(emmeans)
library(gamm4)

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
                                        "S9", "S12", "post1w", "postctrl"))) %>%
  print()




# Control vs. intervention in total RNA per mg tissue

# The same setup for the data set is used as with other variables (qPCR, US, strength).
# The tx variable is used to model combinations of time, group and detrain factors. 

rna_data2 <- rna_complete %>%
  filter(time %in% c("S0", "S1", "S1c",  "S12", "postctrl", "post1w")) %>%
  mutate(Time = if_else(time == "S1c", "S1", 
                        if_else(time %in% c("post1w", "postctrl", "S12"), 
                                "post", as.character(time))), 
         tw = tissue_weight - mean(tissue_weight), 
         time = factor(Time, levels = c("S0", "S1", "post")), 
         group = if_else(cond == "ctrl_leg", "con", "int"), 
         tx = paste0(time, "_", group, "_", detrain), 
         tx = factor(tx, levels = c("S0_con_train", 
                                    "S1_con_train", 
                                    "post_con_train", 
                                    "S0_int_train", 
                                    "S1_int_train", 
                                    "post_int_train", 
                                    "post_int_detrain"))) %>%

  # Values are averaged over leg/time-point
  
  
    group_by(participant, time, leg, group, tx) %>%
  summarise(rna = mean(rna, na.rm = TRUE), 
            tissue_weight = mean(tissue_weight)) %>%

  print()
  




m1 <- lmer(log(rna/tissue_weight) ~ tx + (1|participant), 
           data = rna_data2)


# Combine data for plotting 

# Emmeans can be used as the contrats specified are combinations of factors. 
# emmeans::emmeans is used to create a emmGrid object:

em <- emmeans(m1, specs =  ~ tx)

# To create contrasts, each levels can be specified by a indicator vector. 

S0_con_train    =  c(1, 0, 0, 0, 0, 0, 0)
S1_con_train    =  c(0, 1, 0, 0, 0, 0, 0) 
post_con_train  =  c(0, 0, 1, 0, 0, 0, 0) 

S0_int_train   =  c(0, 0, 0, 1, 0, 0, 0)
S1_int_train   =  c(0, 0, 0, 0, 1, 0, 0) 
post_int_train =  c(0, 0, 0, 0, 0, 1, 0) 
post_int_detrain =  c(0, 0, 0, 0, 0, 0, 1) 


# these can be combined in a named list to produce the contrasts of interest. 


comp_rna <- contrast(em, 
                    method = list("con_S1"     = S1_con_train -  S0_con_train, 
                                  "con_post"    = post_con_train - S0_con_train, 
                                  "int_S1"      =  S1_int_train - S0_int_train,
                                  "int_post"    = post_int_train - S0_int_train,
                                  "int_post1w"  = post_int_detrain - S0_int_train,
                                  "inter:S1"    = (S1_int_train - S0_int_train) -  (S1_con_train - S0_con_train), 
                                  "inter:post"  = (post_int_train - S0_int_train) -  (post_con_train - S0_con_train),
                                  "inter:post1w" = (post_int_detrain - S0_int_train) -  (post_con_train - S0_con_train))) %>%
  confint() %>%
  data.frame() %>%
  print()


#### brms solution ######################################


m1.brm <- brm(log(rna/tissue_weight) ~ tx + (1|participant), 
             
             warmup = 1000, # number of samples before sampling
             iter = 4000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 123, # seed for reproducible results
             save_all_pars = TRUE,
             data = rna_data2)


plot(m1.brm)
pp_check(m1.brm)


summary(m1.brm)

h.rna <- hypothesis(m1.brm, c("txS1_con_train = 0", 
                             "txpost_con_train = 0",
                             "txS1_int_train = txS0_int_train", 
                             "txpost_int_train = txS0_int_train", 
                             "txpost_int_detrain = txS0_int_train", 
                             "txS1_int_train - txS0_int_train  = txS1_con_train", 
                             "txpost_int_train - txS0_int_train  = txpost_con_train",
                             "txpost_int_detrain - txS0_int_train  = txpost_con_train"))

  

comp_rna_bayes <- data.frame(time_group = c("con_S1",                
                                            "con_post",              
                                            "int_S1",                
                                            "int_post",              
                                            "int_post1w",      
                                            "inter:S1",          
                                            "inter:post",        
                                            "inter:post1w")) %>%
  cbind(
    data.frame(h.rna$hypothesis)
  ) %>%
  mutate(estimate = Estimate, lower.CL = CI.Lower, upper.CL = CI.Upper) %>%
  dplyr::select(time_group, estimate, lower.CL, upper.CL) %>%
  print()



# Save data for plotting 
saveRDS(comp_rna_bayes, "./data/derivedData/total-rna-analysis/comp_rna_bayes.RDS")
saveRDS(comp_rna, "./data/derivedData/total-rna-analysis/comp_rna_freq.RDS")





# Time course data in experimental group #

  
time_course <- rna_complete %>%
  filter(cond != "ctrl_leg") %>%
  
  
  
  group_by(participant,leg, time,time.c, cond, ) %>%
  summarise(rna = mean(rna, na.rm = TRUE), 
            tissue_weight = mean(tissue_weight, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(participant = factor(participant), 
         cond = factor(cond), 
         rna.tissue = rna/tissue_weight, 
         tw = tissue_weight - mean(tissue_weight), 
         time.c.cent = time.c - 4, 
         detrain = factor(if_else(time == "post1w", "detrain", "train"), 
                          levels = c("train", "detrain"))) %>%
  # Create dummy variables for segmented model 
  
  mutate(time1 = time.c,
         time2 = if_else(time %in% c("S0", "S1", "S4"), 0, time.c - 4),
         time3 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8))  %>%
  print()
  
  
 


# gam models fitted in brms to visualize general patterns in the data. 
# a piecewise model is fitted to estimate slopes in different parts of the time course.   

# Piecewise model "full" model
tc.m1 <- brm(bf(log(rna.tissue) ~ cond  * (time1 + time2 + time3) + detrain + (1|participant)), 
     data = time_course, 
     warmup = 1000, # number of samples before sampling
     iter = 4000,  # number of mcmc iterations 
     cores = 4, # number of cores used in sampling
     chains = 4, # number of chains
     seed = 123, # seed for reproducible results
     save_all_pars = TRUE, 
     control = list(adapt_delta = 0.99))
     

pp_check(tc.m1, type = "ecdf_overlay")
summary(tc.m1)

tc.m2 <- brm(bf(log(rna.tissue) ~ (time1 + time2 + time3) + detrain + (1|participant)), 
             data = time_course, 
             warmup = 1000, # number of samples before sampling
             iter = 4000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 123, # seed for reproducible results
             save_all_pars = TRUE, 
             control = list(adapt_delta = 0.99))


# No interaction effect between any segment and cond. Reduce model to 
# only containing time effects. 
# Slopes in the different segments can be calculated through the hypothesis 
# function. 

hypothesis(tc.m2, c("time2 + time1 = 0")) # The slope (on log scale) between session 4 and 8
hypothesis(tc.m2, c("time2 + time1 +time3 = 0")) # The slope on log scale between session 8 and 12

# Save the model for stats in text.
saveRDS(tc.m2, "./data/derivedData/total-rna-analysis/tc.m2.RDS")

# Save the gam model 
saveRDS(tc.m1, "./data/derivedData/total-rna-analysis/tc.m1.RDS")


pp_check(tc.m2, type = "ecdf_overlay")
summary(tc.m2)


exp(0.05)

### A model for "descriptive" analysis of available data points including other studies using priors
prior_summary(tc.m3)
prior_m3 <- set_prior("normal(0.0168, 0.0426)", class = "b", coef = "time1")



tc.m3 <- brm(bf(log(rna.tissue) ~ (time1 + time2 + time3) + detrain + (1|participant)), 
             data = time_course, 
             prior = prior_m3,
             warmup = 1000, # number of samples before sampling
             iter = 4000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 123, # seed for reproducible results
             save_all_pars = TRUE, 
             control = list(adapt_delta = 0.99))



### Plotting the effect of prior information #####################


preds  <- tc.m2 %>%
  spread_draws(b_Intercept,
               b_time1, 
               b_time2, 
               b_time3,
               b_detraindetrain) %>%
  mutate(
    s0  = b_Intercept, 
    s1  = b_Intercept + (b_time1 * 1), 
    s2  = b_Intercept + (b_time1 * 2), 
    s3  = b_Intercept + (b_time1 * 3),
    s4  = b_Intercept + (b_time1 * 4), 
    s5  = b_Intercept + (b_time1 * 5) + (b_time2 * 1), 
    s6  = b_Intercept + (b_time1 * 6) + (b_time2 * 2), 
    s7  = b_Intercept + (b_time1 * 7) + (b_time2 * 3), 
    s8  = b_Intercept + (b_time1 * 8) + (b_time2 * 4), 
    s9  = b_Intercept + (b_time1 * 9) + (b_time2 * 5) + (b_time3 * 1),
    s10 = b_Intercept + (b_time1 * 10) + (b_time2 * 6) + (b_time3 * 2),
    s11 = b_Intercept + (b_time1 * 11) + (b_time2 * 7) + (b_time3 * 3),
    s12 = b_Intercept + (b_time1 * 12) + (b_time2 * 8) + (b_time3 * 4), 
    s13 = b_Intercept + (b_time1 * 12) + (b_time2 * 8) + (b_time3 * 4) + b_detraindetrain) %>%
  dplyr::select(.draw, s0:s13) %>%
  pivot_longer(names_to =  "time", values_to = "prediction", 
               cols = s0:s13) %>%
  mutate(time.c = as.numeric(gsub("s", "", time)), 
         group = if_else(time.c == 13, paste0("detrain"), paste0("train"))) %>%
  group_by(time.c, group) %>%
  summarise(m = median(prediction), 
            lwr =  quantile(prediction, 0.025), 
            upr = quantile(prediction, 0.975)) %>%
  print()


preds %>%
  filter(time.c < 12.5) %>% 
  ggplot(aes(time.c, exp(m), group = paste(group))) + 
  
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)), alpha = 0.15) +
  geom_line(aes(), size = 1.5) + 
  
  geom_errorbar(data = filter(preds, time.c > 12.5),
                aes(time.c, exp(m), ymin = exp(lwr), ymax = exp(upr)), width = 0.2) + 
  
  geom_point(data = filter(preds, time.c > 12.5),
             aes(time.c, exp(m)), size = 2.5, shape = 21) 








summary(tc.m3)
summary(tc.m2)

pp_check(tc.m3)

saveRDS(cond_eff_rna_tc, "./data/derivedData/total-rna-analysis/cond_eff_rna_tc.RDS")




