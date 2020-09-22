##-------------------------------------
## ubf_tot_rna_model.R
##
## Title: 
## Purpose: Model total RNA over time and influence of UBF levels and effect of rRNA on muscle hypertrophy
## Author:
##
##
##
##-------------------------------------
## Notes:
# UBF is believed to be a contributor to rDNA transcription by recruiting 
# the SL1 transcription factor to the rDNA promoter. mTOR increases UBF 
# phosphorylation (rapamycin sensitive) increasing its activity and c-Myc induces
# UBTF transcription increasing protein levels (rapamycin insensitive). 
# 
# Below we want to test if rRNA levels are assiciated with UBF levels. Both 
# are increasing with time and therefore, the model controls for the effect of training (time).
#
# This script also contain a model of muscle growth were we estimate the effect of increase in rRNA
# on muscle growth (muscle thicknes change pre- to post training). Post training values are 
# from both S12 and post 1w (de-training). rRNA transcription is estimated as the average linear increase over 
# all sessions for each individual and leg. 
#
#
## ------------------------------------

# Load packages
library(tidyverse)
library(brms)
library(bayesplot)
library(lme4)
library(readxl)




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


# Load western blot data 

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


# Combined data from both Total-RNA and western blot analysis

combined_df  <- west_cont %>%
  group_by(target, participant, cond, time) %>%
  summarise(scaled = mean(expression.scaled, na.rm = TRUE), 
            raw = mean(expression)) %>%
  mutate(target = if_else(target == "t-s6", "rps6", "ubf")) %>%
  
  pivot_wider(names_from = target, values_from = c(scaled, raw)) %>%

  ungroup() %>%
  dplyr::select(participant, cond, time, scaled_rps6:raw_ubf) %>%
  inner_join(rna_complete %>%
               mutate(rna.mg = rna/tissue_weight) %>%
               group_by(participant, leg, time,time.c, time1, time2, time3, cond) %>%
               summarise(rna.mg = mean(rna.mg, na.rm = TRUE))) %>%
  mutate(detrain = factor(if_else(time == "post1w", "detrain", "train"), 
                          levels = c("train", "detrain")), 
         log.rna = log(rna.mg)) %>%
  filter(cond != "ctrl_leg") %>%
  print()


saveRDS(combined_df, "./data/derivedData/ubf-tot-rna-model/combined_df.RDS")



### Combined df in brms models ##############################

# This model combines several aspects of our data. We are interested in 
# the effect of UBF levels. UBF is scaled and should be interpreted in units of 
# standard deviation. The time variables are used to capture three segmenst of the data 
# session 1-4, 48, and 8-12. The detraining models the effect of detraining.

# Preliminary model did not reveal any interaction effect between UBF and time or detraining
# meaning that the effect of UBF levels are similar over the course of the study. 

### The ubf model
comb.m1 <- brm(bf(log.rna ~ scaled_ubf + (time1 +time2 + time3) + detrain + (1|participant/leg)),
              
              warmup = 1000, # number of samples before sampling
              iter = 5000,  # number of mcmc iterations 
              cores = 4, # number of cores used in sampling
              chains = 4, # number of chains
              seed = 5, # seed for reproducible results
              #    control = list(adapt_delta = 0.95), 
              data = combined_df)

#### The rpS6 model


comb.s6.m1 <- brm(bf(log.rna ~ scaled_rps6 + (time1 +time2 + time3) + detrain + (1|participant/leg)),
               
               warmup = 1000, # number of samples before sampling
               iter = 5000,  # number of mcmc iterations 
               cores = 4, # number of cores used in sampling
               chains = 4, # number of chains
               seed = 5, # seed for reproducible results
               #    control = list(adapt_delta = 0.95), 
               data = combined_df)



# A naive model not accounting for time 
# comb.m1 <- brm(bf(log.rna ~ ubf.scaled + time + (1 + time1|participant)),
#                
#                warmup = 1000, # number of samples before sampling
#                iter = 4000,  # number of mcmc iterations 
#                cores = 4, # number of cores used in sampling
#                chains = 4, # number of chains
#                seed = 5, # seed for reproducible results
#                #    control = list(adapt_delta = 0.95), 
#                data = combined_df)
# 

### Check how UBF affects total RNA levels 

# Models containg all participants
# Naive model (without time)

# m1 <- lmer(log.rna ~ ubf.scaled + (1|participant), 
#      data = combined_df)
# 
# # Continuous time
# m2 <- lmer(log.rna ~ ubf.scaled + (time1 + time2 + time3) + detrain + (1|participant), 
#            data = combined_df)
# # Time factor
# m3 <- lmer(log.rna ~ ubf.scaled + time  + (1|participant), 
#            data = combined_df)
# 
# 
# # Models with only TRAIN group
# m1x <- lmer(log.rna ~ ubf.scaled + (1|participant), 
#            data = combined_df[combined_df$cond != "ctrl_leg",])
# 
# # Continuous time
# m2x <- lmer(log.rna ~ ubf.scaled + (time1 + time2 + time3) + detrain + (1 + time1|participant), 
#             data = combined_df[combined_df$cond != "ctrl_leg",])
# # Time factor
# m3x <- lmer(log.rna ~ ubf.scaled + time  + (1|participant), 
#             data = combined_df[combined_df$cond != "ctrl_leg",])
# 
# summary(m2x)
# 
# 
# 
# # Combine estimates 
# 
# df <- data.frame(model = c("naive", "continuous", "factor", "naive_TRAIN", "continuous_TRAIN", "factor_TRAIN"), 
#            est = c(summary(m1)$coef[2,1], 
#                    summary(m2)$coef[2,1], 
#                    summary(m3)$coef[2,1], 
#                    summary(m1x)$coef[2,1], 
#                    summary(m2x)$coef[2,1], 
#                    summary(m3x)$coef[2,1]), 
#            
#                    lwr = c(confint(m1)[4, 1], 
#                            confint(m2)[4, 1], 
#                            confint(m3)[4, 1], 
#                            confint(m1x)[4, 1], 
#                            confint(m2x)[6, 1], 
#                            confint(m3x)[4, 1]), 
#                    upr = c(confint(m1)[4, 2], 
#                            confint(m2)[4, 2], 
#                            confint(m3)[4, 2], 
#                            confint(m1x)[4, 2], 
#                            confint(m2x)[6, 2], 
#                            confint(m3x)[4, 2]))
# 
# df %>%
#   ggplot(aes(x = est, y = model)) +
#            geom_errorbarh(aes(xmin = lwr, xmax = upr)) +
#            geom_point() 
# 
# plot(m2)
# 



# Model checks:

# Posterior predictive checks
pp_check(comb.m1, type = "ecdf_overlay")
pp_check(comb.m1, type = "stat", stat = "median")
# The data is in agreement with what the model predicts, good.

pp_check(comb.s6.m1, type = "ecdf_overlay")
pp_check(comb.s6.m1, type = "stat", stat = "median")

# Leave one out statistics
loo(comb.m1)
loo(comb.s6.m1)
# All Pareto estimates are OK


# Plot chains for convergance and RE estimates

plot(comb.m1)

# All effects show convergance and are interpretable (e.g. RE > 0).

summary(comb.m1)


conditional_effects(comb.m1)



# Save model 
saveRDS(comb.m1, "./data/derivedData/ubf-tot-rna-model/comb_m1.RDS")
saveRDS(comb.s6.m1, "./data/derivedData/ubf-tot-rna-model/comb_m1.RDS")




 ### Estimation of effects of total RNA on muscle growth ###########################

# A function to test for colinearity

# From https://jonlefcheck.net/2012/12/28/dealing-with-multicollinearity-using-variance-inflation-factors/

# Notes:

# Below we will fit models with the intercept and slope from the rna to number of sessions
# relationship for each leg. These could be related (greater intercept smaller slope). This in turn
# could lead to multicolinearity. 
# Currently there seems to be no method to estimate the variance inflation factor in 
# brms models. Paul Buerkner suggests fitting models in lm(/lmer) to estimate VIF to diagnose any problems
# https://discourse.mc-stan.org/t/alternative-to-car-vif-for-brmsfit-models/3970

# A function to calculate VIF from lmer model objects. 

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}



# Load muscle thickness data.

# The dependent variable will be mm change in muscle thickness. We will use
# estimates from both post and post1w (the model will average out/estimates potential effects).

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

# A exploitative plot of linear fits for each participant (RNA to number of sessions).

combined_df %>%
filter(!(time %in% c("post1w")), 
       cond != "ctrl_leg") %>% # Removing the de-training estimate as doeas not represent training induced increase.

  
 # filter(!(participant == "P21" & time == "S12" & leg == "R")) %>%
  
  mutate(time.c = time.c - 6) %>% # this centers number of sessions and the intercept becomes
  # the estimate total RNA (on log scale) in the mid of the training intervention.
  
  ggplot(aes(time.c, rna.mg, group = paste(participant, leg))) + geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ participant)


## UBF 
combined_df %>%
  filter(!(time %in% c("post1w")), 
         cond != "ctrl_leg") %>% # Removing the de-training estimate as doeas not represent training induced increase.
  
  
  # filter(!(participant == "P21" & time == "S12" & leg == "R")) %>%
  
  mutate(time.c = time.c - 6) %>% # this centers number of sessions and the intercept becomes
  # the estimate total RNA (on log scale) in the mid of the training intervention.
  
  ggplot(aes(time.c, ubf, group = paste(participant, leg))) + geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ participant)







predict_df <- combined_df %>%
  filter(!(time %in% c("post1w")), 
         cond != "ctrl_leg") %>% # Removing the de-training estimate as doeas not represent training induced increase.
  
  # Change the intercept to represent average total-RNA. This gives more meaning to the intercept 
  # and it could be use for prediction. It will also remove some of the potential correlation between 
  # the intercept and slope.
  
  mutate(time.c = time.c - 6 ) %>% # this centers number of sessions and the intercept becomes
  # the estimate total RNA (on log scale) in the mid of the training intervention.

  group_by(participant, leg) %>%
  
  # A simple linear model is used for each participant. 
  # A maximal and mean residual for each model is retained
  # for checking if this influences the model in any meaningfull 
  # way.
  
  do(intercept = coef(lm(log(rna.mg) ~ time.c, data = .))[1],
     slope =     coef(lm(log(rna.mg) ~ time.c, data = .))[2], 
     max.resid = max(resid(lm(log(rna.mg) ~ time.c, data = .))^2), 
     mean.resid =mean(resid(lm(log(rna.mg) ~ time.c, data = .))^2), 
     ubf.slope = coef(lm(log(ubf) ~ time.c, data = .))[2]) %>% 

  unnest(c(intercept, slope, max.resid, mean.resid, ubf.slope)) %>%
  
  
# Join with the ultra sound data 
  
  inner_join(us_data %>%
               filter(cond != "ctrl_leg") %>%
               dplyr::select(-time_pp, -detrain) %>%
               pivot_wider(names_from = time, values_from = thickness) %>%
               
               # Calculates differences in muscle thickness.
               
               mutate(post = post - pre, 
                      post1w = post1w - pre) %>%  
               group_by(participant) %>%
              
               
              pivot_longer(names_to = "time", values_to = "mm_incr", cols = post:post1w) ) %>%
  # Scales slope to represent percentage increase per session
  # Scales the intercept to represent sd from the mean on the natural scale.
mutate(slope =  100 * (exp(slope)-1), 
         intercept = scale(exp(intercept))) %>%
  print()
  

## Save the full df for plotting 

saveRDS(predict_df, "./data/derivedData/ubf-tot-rna-model/predict_df.RDS")



# Explorative plot of the relationship between slope and intercept 

predict_df %>%
  ggplot(aes(ubf.slope, slope)) + geom_point() +
  geom_smooth(method = "lm")

# No relationship between slope and intercept (RNA). This is a good sign to avoid 
# colinearity. But there is a tendency toward positive correlation 
# between UBF and RNA slopes. We next check if this is a problem in the biggest model and subsequent 
# simpler ones. 
# Fitting is done in lme4 first to check potential issues. 

# Find participantrs with large deviations from the linera model, do the affect 
# estimates


predict_df %>%
  ggplot(aes(ubf.slope, slope, label = paste(participant, leg))) + geom_point() +
  ggrepel::geom_text_repel() + geom_smooth(method = "lm")




predict_df %>%
  ggplot(aes(participant, max.resid)) + geom_point()  


m0 <- lmer(mm_incr  ~ pre +  ubf.slope + slope +  (1|participant), 
           data = predict_df, REML = FALSE)


m1 <- lmer(mm_incr  ~ pre  + (1|participant),
           data = predict_df, REML = FALSE) 


summary(m0)
anova(m0, m1)


m2 <- lmer(mm_incr  ~  pre +  intercept + slope + ubf.slope + (1|participant/leg), 
           data = predict_df, REML = FALSE) 
m2.2 <- lmer(mm_incr  ~ time +  sex +  intercept + slope + (1|participant/leg), 
           data = predict_df, REML = FALSE) 
m3 <- lmer(mm_incr  ~   intercept + slope + (1|participant/leg), 
           data = predict_df, REML = FALSE) 

summary(m2)


plot(m2.2)

# Apply the vif mer function to each model to check colinearity. 
# a threshold of 2 is considered meaningfull. 
lapply(list(m1 = m1, m2 = m2, m2.2 = m2.2, m3 =  m3), vif.mer)

# No model contains any problematic terms. 

anova(m1 ,m2, m2.2,  m3)

# Preliminary models suggests that sex/pre should be included in the modeling 

# Fitting with brms
summary(m1)






hyp.pre.m1 <- brm(bf(mm_incr  ~   sex + intercept + slope + (1|participant/leg)),
                  
               warmup = 4000, # number of samples before sampling
               iter = 10000,  # number of mcmc iterations 
               cores = 4, # number of cores used in sampling
               thin = 10, 
               chains = 4, # number of chains
               seed = 5, # seed for reproducible results
               control = list(adapt_delta = 0.95), 
               save_all_pars = TRUE, 
               data = predict_df)


summary(hyp.pre.m1)

hyp.pre.m2 <- brm(bf(mm_incr  ~  time + pre  + sex +  intercept + slope + (1|participant/leg)),
           
                  warmup = 4000, # number of samples before sampling
                  iter = 10000,  # number of mcmc iterations 
                  cores = 4, # number of cores used in sampling
                  thin = 10, 
                  chains = 4, # number of chains
                  seed = 5, # seed for reproducible results
                  control = list(adapt_delta = 0.95), 
                  save_all_pars = TRUE, 
                  data = predict_df)

summary(hyp.pre.m2)


hyp.pre.m3 <- brm(bf(mm_incr  ~   sex +  intercept + slope + (1|participant/leg)),
                  
                  warmup = 4000, # number of samples before sampling
                  iter = 10000,  # number of mcmc iterations 
                  cores = 4, # number of cores used in sampling
                  thin = 10, 
                  chains = 4, # number of chains
                  seed = 5, # seed for reproducible results
                  control = list(adapt_delta = 0.95), 
                  save_all_pars = TRUE, 
                  data = predict_df)


# Model checks:

# Posterior predictive checks
pp_check(hyp.pre.m1, type = "ecdf_overlay")
pp_check(hyp.pre.m1, type = "stat", stat = "median")
# The data is in agreement with what the model predicts, good.



# Leave one out statistics
loo(hyp.pre.m1, moment_match = TRUE)
loo(hyp.pre.m2, moment_match = TRUE)
loo(hyp.pre.m3, moment_match = TRUE)
# All Pareto estimates are OK


# Plot chains for convergance and RE estimates
plot(hyp.pre.m1)

# All effects show convergance and are interpretable (e.g. RE > 0).

summary(hyp.pre.m1)




# Save model 
saveRDS(hyp.pre.m1, "./data/derivedData/ubf-tot-rna-model/hyp_pre_m1.RDS")
saveRDS(hyp.pre.m2, "./data/derivedData/ubf-tot-rna-model/hyp_pre_m2.RDS")
saveRDS(hyp.pre.m3, "./data/derivedData/ubf-tot-rna-model/hyp_pre_m3.RDS")



summary(hyp.pre.m3)




#### Robust estimates of rna as a predictor ###################

# Above we have estimated the effect of total rna increase on muscle hypertrophy. 
# To see if these results are robust to ouliers, random exclusion of samples from 
# slopes estimates are used to estimate how many times we will find the observed result.

# The permutation is done on the level of RNA per mg over time estimate, i.e. the 
# estimate of the slope and intercept models for each participant. 

# To do this a single (or more) samples will be removed on each iteration and estimates
# of slope and intercept will be made together with estimates of effect of slope on muscle
# hypertrophy.

# A second leave one out iteration will be done on the level of muscle thickness.
# (Legs/Participant)



# List all samples for permutations


samples <- combined_df %>%
  filter(!(time %in% c( "post1w")), 
         cond != "ctrl_leg") %>%
  mutate(samples = paste(sep = "_", participant, leg, time)) %>%
  pull(samples)

# Save results
results <- list()

for(i in 1:length(samples)) {
  
  
  # Prepare data set
  
  predict_df <- combined_df %>%
    filter(!(time %in% c( "post1w"))) %>% # Removing the de-training estimate as does not represent training induced increase.
    
    # Change the intercept to represent average total-RNA. This gives more meaning to the intercept 
    # and it could be use for prediction. It will also remove some of the potential correlation between 
    # the intercept and slope.
    
    mutate(samples = paste(sep = "_", participant, leg, time), 
           time.c = time.c - 6 ) %>% # this centers number of sessions and the intercept becomes
    # the estimate total RNA (on log scale) in the mid of the training intervention.
    
    # Filter a specific sample in the i:th permutation
    filter(samples != samples[i]) %>%

    
    group_by(participant, leg) %>%
    
    # A simple linear model is used for each participant. 
    # A maximal and mean residual for each model is retained
    # for checking if this influences the model in any meaningfull 
    # way.
    
    do(intercept = coef(lm(log(rna.mg) ~ time.c, data = .))[1],
       slope =     coef(lm(log(rna.mg) ~ time.c, data = .))[2], 
       max.resid = max(resid(lm(log(rna.mg) ~ time.c, data = .))^2), 
       mean.resid =mean(resid(lm(log(rna.mg) ~ time.c, data = .))^2)) %>% 
    
    unnest(c(intercept, slope, max.resid, mean.resid)) %>%
    
    # Join with the ultra sound data 
    
    inner_join(us_data %>%
                 filter(cond != "ctrl_leg") %>%
                 dplyr::select(-time_pp, -detrain) %>%
                 pivot_wider(names_from = time, values_from = thickness) %>%
                 
                 # Calculates differences in muscle thickness.
                 
                 mutate(post = post - pre, 
                        post1w = post1w - pre) %>%  
                 group_by(participant) %>%
                 
                 
                 pivot_longer(names_to = "time", values_to = "mm_incr", cols = post:post1w) ) %>%
    # Scales slope to represent percentage increase per session
    # Scales the intercept to represent sd from the mean on the natural scale.
    mutate(slope =  100 * (exp(slope)-1), 
           intercept = scale(exp(intercept))) 
  
  
  m1 <- lmer(mm_incr  ~  slope + (1|participant/leg),
             data = predict_df, REML = FALSE) 
 
  m2 <- lmer(mm_incr  ~  sex + slope + (1|participant/leg),
             data = predict_df, REML = FALSE) 
  
  m3 <- lmer(mm_incr  ~  sex + slope + intercept + (1|participant/leg),
             data = predict_df, REML = FALSE) 
  
  
  
resdf <-   data.frame(model = c("m1", rep("m2", 2), rep("m3", 3)), 
             sample = samples[i], 
             coef = c("slope", "sex", "slope", "sex", "slope", "intercept"), 
             estimate = c(summary(m1)$coef[2, 1], 
                          summary(m2)$coef[2, 1], 
                          summary(m2)$coef[3, 1], 
                          summary(m3)$coef[2, 1], 
                          summary(m3)$coef[3, 1],
                          summary(m3)$coef[4, 1]), 
             t_val =  c(summary(m1)$coef[2, 3], 
                        summary(m2)$coef[2, 3], 
                        summary(m2)$coef[3, 3],
                        summary(m3)$coef[2, 3], 
                        summary(m3)$coef[3, 3],
                        summary(m3)$coef[4, 3]), 
             
             
             lwr =  c(confint.merMod(m1, method = "Wald")[5, 1], 
                      confint.merMod(m2, method = "Wald")[5, 1],
                      confint.merMod(m2, method = "Wald")[6, 1],
                      
                      confint.merMod(m3, method = "Wald")[5, 1],
                      confint.merMod(m3, method = "Wald")[6, 1],
                      confint.merMod(m3, method = "Wald")[7, 1]),
                      
             upr =  c(confint.merMod(m1, method = "Wald")[5, 2],
                      confint.merMod(m2, method = "Wald")[5, 2], 
                      confint.merMod(m2, method = "Wald")[6, 2],
                      confint.merMod(m3, method = "Wald")[5, 2],
                      confint.merMod(m3, method = "Wald")[6, 2],
                      confint.merMod(m3, method = "Wald")[7, 2]))

  results[[i]] <- resdf
  
  
  
}



# List all participants for permutations


participants <- unique(combined_df$participant)

# Save results
results2 <- list()

for(i in 1:length(participants)) {
  
  
  # Prepare data set
  
  predict_df <- combined_df %>%
    filter(!(time %in% c( "post1w"))) %>% # Removing the de-training estimate as does not represent training induced increase.
    
    # Change the intercept to represent average total-RNA. This gives more meaning to the intercept 
    # and it could be use for prediction. It will also remove some of the potential correlation between 
    # the intercept and slope.
    
    mutate(samples = paste(sep = "_", participant, leg, time), 
           time.c = time.c - 6 ) %>% # this centers number of sessions and the intercept becomes
    # the estimate total RNA (on log scale) in the mid of the training intervention.
    
    # Filter a specific sample in the i:th permutation
   filter(participant != participants[i]) %>%
    
    
    group_by(participant, leg) %>%
    
    # A simple linear model is used for each participant. 
    # A maximal and mean residual for each model is retained
    # for checking if this influences the model in any meaningfull 
    # way.
    
    do(intercept = coef(lm(log(rna.mg) ~ time.c, data = .))[1],
       slope =     coef(lm(log(rna.mg) ~ time.c, data = .))[2], 
       max.resid = max(resid(lm(log(rna.mg) ~ time.c, data = .))^2), 
       mean.resid =mean(resid(lm(log(rna.mg) ~ time.c, data = .))^2)) %>% 
    
    unnest(c(intercept, slope, max.resid, mean.resid)) %>%
    
    # Join with the ultra sound data 
    
    inner_join(us_data %>%
                 filter(cond != "ctrl_leg") %>%
                 dplyr::select(-time_pp, -detrain) %>%
                 pivot_wider(names_from = time, values_from = thickness) %>%
                 
                 # Calculates differences in muscle thickness.
                 
                 mutate(post = post - pre, 
                        post1w = post1w - pre) %>%  
                 group_by(participant) %>%
                 
                 
                 pivot_longer(names_to = "time", values_to = "mm_incr", cols = post:post1w) ) %>%
    # Scales slope to represent percentage increase per session
    # Scales the intercept to represent sd from the mean on the natural scale.
    mutate(slope =  100 * (exp(slope)-1), 
           intercept = scale(exp(intercept))) 
  
  
  m1 <- lmer(mm_incr  ~  slope + (1|participant/leg),
             data = predict_df, REML = FALSE) 
  
  m2 <- lmer(mm_incr  ~  sex + slope + (1|participant/leg),
             data = predict_df, REML = FALSE) 
  
  m3 <- lmer(mm_incr  ~  sex + slope + intercept + (1|participant/leg),
             data = predict_df, REML = FALSE) 
  
  
  
  resdf <-   data.frame(model = c("m1", rep("m2", 2), rep("m3", 3)), 
                        participant = participants[i], 
                        coef = c("slope", "sex", "slope", "sex", "slope", "intercept"), 
                        estimate = c(summary(m1)$coef[2, 1], 
                                     summary(m2)$coef[2, 1], 
                                     summary(m2)$coef[3, 1], 
                                     summary(m3)$coef[2, 1], 
                                     summary(m3)$coef[3, 1],
                                     summary(m3)$coef[4, 1]), 
                        t_val =  c(summary(m1)$coef[2, 3], 
                                   summary(m2)$coef[2, 3], 
                                   summary(m2)$coef[3, 3],
                                   summary(m3)$coef[2, 3], 
                                   summary(m3)$coef[3, 3],
                                   summary(m3)$coef[4, 3]), 
                        
                        
                        lwr =  c(confint.merMod(m1, method = "Wald")[5, 1], 
                                 confint.merMod(m2, method = "Wald")[5, 1],
                                 confint.merMod(m2, method = "Wald")[6, 1],
                                 
                                 confint.merMod(m3, method = "Wald")[5, 1],
                                 confint.merMod(m3, method = "Wald")[6, 1],
                                 confint.merMod(m3, method = "Wald")[7, 1]),
                        
                        upr =  c(confint.merMod(m1, method = "Wald")[5, 2],
                                 confint.merMod(m2, method = "Wald")[5, 2], 
                                 confint.merMod(m2, method = "Wald")[6, 2],
                                 confint.merMod(m3, method = "Wald")[5, 2],
                                 confint.merMod(m3, method = "Wald")[6, 2],
                                 confint.merMod(m3, method = "Wald")[7, 2]))
  
  results2[[i]] <- resdf
  
  
  
}


# The different model formulas were 

# m1 : mm_incr ~ slope
# m2 : mm_incr ~ sex + slope
# m3 : mm_incr ~ sex + slope + intercept


res <- bind_rows(results) %>%
  separate(sample, into = c("participant", "leg", "time")) 



res2 <- bind_rows(results2) 



saveRDS(list(loo.participant = res2 , 
             loo.sample = res), "./data/derivedData/ubf-tot-rna-model/leave-one-out.RDS")



























