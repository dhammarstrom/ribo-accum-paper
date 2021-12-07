### Western blot analysis  ######

# Load packages ###########

library(tidyverse)
library(brms)
library(bayesplot)



# Western blot data are compiled in western_compile.R and stored
# gel information is included in the raw data but averaged over all gels
# below. 


western_data <- readRDS("./data/derivedData/western-compile/western_data.RDS")

# Re-format time variables to match comparison between control and 
# experimental groups. The "tx" factor enables these comparisons. The 
# tx factor is needed to avoid "rank deficiency" in the comparison.

west.dat <- western_data %>%
 filter(time %in% c("S0", "S1c", "S1", "S12", "postctrl", "post1w")) %>%
  mutate(group = if_else(cond == "ctrl_leg", "con", "int"), 
         detrain = if_else(time == "post1w", "detrain", "train"), 
         time = if_else(time == "S0", "pre", 
                        if_else(time == "S1c", "S1", 
                                if_else(time %in% c("postctrl", "post1w", "S12"), "post", time))), 
         time = factor(time, levels = c("pre", "S1", "post")), 
         tx = paste(time, group, detrain, sep = "_"), 
         tx = factor(tx, levels = c("pre_con_train",
                                    "S1_con_train",
                                    "post_con_train",
                                    "pre_int_train",
                                    "S1_int_train",
                                    "post_int_train", 
                                    "post_int_detrain")), 
         sample2 = paste0(participant, leg, time)) %>%
  group_by(participant, tx, target) %>%
  summarise(expression = mean(expression, na.rm = TRUE)) %>%
  mutate(log.expr = log(expression)) %>%
  filter(!is.na(expression)) %>%
  print()
  

s6.m1 <- brm(bf(log.expr ~ tx + (1|participant)),

          warmup = 1000, # number of samples before sampling
          iter = 4000,  # number of mcmc iterations 
          cores = 4, # number of cores used in sampling
          chains = 4, # number of chains
          seed = 5, # seed for reproducible results
          save_pars = save_pars(all = TRUE),
      #    control = list(adapt_delta = 0.95), 
          data = west.dat[west.dat$target == "t-s6",])


ubf.m1 <- brm(bf(log.expr ~ tx + (1|participant)),
             
             warmup = 1000, # number of samples before sampling
             iter = 4000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 5, # seed for reproducible results
             save_pars = save_pars(all = TRUE),
             #    control = list(adapt_delta = 0.95), 
             data = west.dat[west.dat$target == "t-UBF",])





# Save model 
saveRDS(s6.m1, "./data/derivedData/western-analysis/s6_tx_m1.RDS")
saveRDS(ubf.m1, "./data/derivedData/western-analysis/ubf_tx_m1.RDS")


# s6.m1x <- readRDS( "./data/derivedData/western-analysis/s6_tx_m1.RDS")
# ubf.m1x <- readRDS( "./data/derivedData/western-analysis/ubf_tx_m1.RDS")



# Model checks

# predictive posterior checks
pp_check(s6.m1, type = "stat", stat = "median")
pp_check(s6.m1, type = "ecdf_overlay")

pp_check(ubf.m1, type = "stat", stat = "median")
pp_check(ubf.m1, type = "ecdf_overlay")

# Leave one out diagnostics
loo(s6.m1, moment_match = TRUE)
loo(ubf.m1, moment_match = TRUE)

# Note: pp-checks looks ok. 
# Some observations had an "ok" pareto diagnostic value. Otherwise "good"

# Continuous data between conditions ##########

# Transform the time variable to a continuous variable covering the number of sets
# and further create variables for the piece-wise model.
# Expression values are averaged over gels.

west_cont <- western_data %>%
  filter(cond != "ctrl_leg", 
         time %in% c("S0", 
                     "S1", 
                     "S4", 
                     "S5", 
                     "S8", 
                     "S9", 
                     "S12")) %>%
  mutate(participant = factor(participant), 
       cond = factor(cond), 
       target = factor(target),
       gel = factor(gel),
       time = factor(time, levels = c("S0", 
                                      "S1", 
                                      "S4", 
                                      "S5", 
                                      "S8", 
                                      "S9", 
                                      "S12")), 
       time.c = as.numeric(gsub("S", "", time)), 
       time1 = time.c,
       time2 = if_else(time %in% c("S0", "S1", "S4"), 0, time.c - 4),
       time3 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8) ) %>%
  
  group_by(target, participant, time.c, cond, time, time1, time2, time3) %>%
  summarise(expression = mean(expression)) %>%
  mutate(log.expr = log(expression)) %>%
  filter(!is.na(expression)) %>%
  
  print()

# Model fits

s6.m1 <- brm(bf(log.expr ~ cond * (time1 + time2 + time3) + (1|participant)), 
             data = west_cont[west_cont$target == "t-s6",], 
             warmup = 1000, # number of samples before sampling
             iter = 4000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 123, # seed for reproducible results
             save_pars = save_pars(all = TRUE),
             control = list(adapt_delta = 0.99))

s6.m2 <- brm(bf(log.expr ~  (time1 + time2 + time3) + (1|participant)), 
             data = west_cont[west_cont$target == "t-s6",], 
             warmup = 1000, # number of samples before sampling
             iter = 4000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 123, # seed for reproducible results
             save_pars = save_pars(all = TRUE),
             control = list(adapt_delta = 0.99))

ubf.m1 <- brm(bf(log.expr ~ cond * (time1 + time2 + time3) + (1|participant)), 
             data = west_cont[west_cont$target == "t-UBF",], 
             warmup = 1000, # number of samples before sampling
             iter = 4000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 123, # seed for reproducible results
             save_pars = save_pars(all = TRUE),
             control = list(adapt_delta = 0.99))



# Model checks

# predictive posterior checks
pp_check(s6.m1, type = "stat", stat = "median")
pp_check(s6.m1, type = "ecdf_overlay")

pp_check(ubf.m1, type = "stat", stat = "median")
pp_check(ubf.m1, type = "ecdf_overlay")

# Leave one out diagnostics
s6_diagnostics <- loo(s6.m1, moment_match = TRUE)
ubf_diagnostics <- loo(ubf.m1, moment_match = TRUE)




### Save models
saveRDS(s6.m1, "./data/derivedData/western-analysis/s6_segmented.RDS")
saveRDS(ubf.m1, "./data/derivedData/western-analysis/ubf_segmented.RDS")

s6.m1  <- readRDS("./data/derivedData/western-analysis/s6_segmented.RDS")
ubf.m1 <- readRDS( "./data/derivedData/western-analysis/ubf_segmented.RDS")





## Create a hypothesis vector for multiple comparisons
# All hypothesis may be "tested" from the same model as well as get 
# average estimates for slopes. 


comps_protein_segmented <- c(
  # Differences between slopes in the two volume conditions
  "condvar:time1 = 0",
  "time2 + time1  = time2 + time1  + condvar:time2 + condvar:time1", #  between session 4 and 8
  "time2 + time1 + time3 = time2 + time1 + time3 + condvar:time2 + condvar:time1 + condvar:time3",  # The slope on log scale between session 8 and 12
  # Differences in slopes average slopes in different segments
  "time1 + (condvar:time1/2) + time2 + (condvar:time2/2) = time1 + (condvar:time1/2)", 
  "time1 + (condvar:time1/2) + time2 + (condvar:time2/2) + time3 + (condvar:time3/2) = time1 + (condvar:time1/2)" ,
  # Estimates of each average slope 
  "time1 + (condvar:time1/2) = 0",
  "time1 + (condvar:time1/2) + time2 + (condvar:time2/2) = 0", 
  "time1 + (condvar:time1/2) + time2 + (condvar:time2/2) + time3 + (condvar:time3/2) = 0"
  ) 


# Create hypothesis output
s6.hyp <- hypothesis(s6.m1, comps_protein_segmented)
ubf.hyp <- hypothesis(ubf.m1, comps_protein_segmented)


# Compare average slopes to the non-interaction model to check this!
hypothesis(s6.m2, c("time1 = 0", 
                    "time1 + time2 = 0", 
                    "time1 + time2 + time3 = 0"))
# Average slopes are similar to the non interaction model (probably differs only by mcmc error)



# Compile data sets for data presentation 

segmented_protein_results <- cbind(data.frame(target = rep(c("rpS6", "UBF"), each = 8), 
           comparison = rep(c("seg1:condvar", 
                              "seg2:condvar", 
                              "seg3:condvar", 
                              "seg2 vs. seg1", 
                              "seg3 vs. seg1", 
                              "avg.est seg1", 
                              "avg.est seg2", 
                              "avg.est seg3"), 2), 
           description = rep(c("Diff btw cond in seg 1", 
                               "Diff btw cond in seg 2", 
                               "Diff btw cond in seg 3", 
                               "Diff in slopes btw seg 1 and 2", 
                               "Diff in slopes btw seg 1 and 3", 
                               "Avg slope in seg 1", 
                               "Avg slope in seg 2", 
                               "Avg slope in seg 3"))), 
      rbind(s6.hyp$hypothesis[,-1], 
            ubf.hyp$hypothesis[,-1]))


### Save compiled results 
saveRDS(segmented_protein_results, "./data/derivedData/western-analysis/segmented_results.RDS")

################ De-training model comparing cond and var ################


# A de-training model is used to estimate differences between 
# volume conditions 

west_detrain <- western_data %>%
  filter(cond != "ctrl_leg", 
         time %in% c("S12", "post1w")) %>%
  mutate(participant = factor(participant), 
         cond = factor(cond), 
         target = factor(target),
         time = factor(time, levels = c("S12", "post1w"))) %>%
  
  group_by(target, participant, cond, time) %>%
  summarise(expression = mean(expression)) %>%
  mutate(log.expr = log(expression)) %>%
  filter(!is.na(expression)) %>%
  
  print()


s6.m3 <- brm(bf(log.expr ~ time * cond + (1|participant)),
             
             warmup = 2000, # number of samples before sampling
             iter = 10000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 5, # seed for reproducible results
             thin = 10,
             control = list(adapt_delta = 0.95), 
             save_pars = save_pars(all = TRUE),
             data = west_detrain[west_detrain$target == "t-s6",])


ubf.m3 <- brm(bf(log.expr ~ time * cond + (1|participant)),
              
              warmup = 2000, # number of samples before sampling
              iter = 10000,  # number of mcmc iterations 
              cores = 4, # number of cores used in sampling
              chains = 4, # number of chains
              seed = 5, # seed for reproducible results
              thin = 10,
              control = list(adapt_delta = 0.95), 
              save_pars = save_pars(all = TRUE),
              data = west_detrain[west_detrain$target == "t-UBF",])



pp_check(ubf.m3, type = "ecdf_overlay")
pp_check(s6.m3, type = "ecdf_overlay")

# Leave one out diagnostics
s6_diagnostics <- loo(s6.m3, moment_match = TRUE)
ubf_diagnostics <- loo(ubf.m3, moment_match = TRUE)


# Save models

saveRDS(s6.m3, "./data/derivedData/western-analysis/s6_detrain.RDS")
saveRDS(ubf.m3, "./data/derivedData/western-analysis/ubf_detrain.RDS")


