# Ultrasound Bayesian/frequentist estimates ###############################

library(brms)
library(tidybayes)
library(readxl)
library(tidyverse)
library(lme4)



# Principally the same model is used for US data. Lesser data points are 
# available in this data set as pre testing is only done once. 


# Combine data sets 
us_data <- rbind(read_excel("./data/ultrasound/ultrasound_data.xlsx") %>%
                   inner_join(read_csv("./data/ultrasound/ultrasound_codekey.csv")) %>%
                   mutate(leg = gsub("VL", "", leg)) %>%
                   inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
                   dplyr::select(participant, time, leg, sex, cond, code, length) %>%
                   mutate(group = if_else(participant %in% paste("P", c(1:7,19:23), sep = ""), "experiment", "control")) %>%
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




## Model between conditions differences 
# 
# m1 <- lmer(thickness ~ time + time:cond + (time|participant), 
#            data = us_data[us_data$group == "experiment",])
# 
# summary(rePCA(m1))
# # Singular, all variation is explained at the intercept, removing correlation reducing slopes
# 
# m2 <- lmer(thickness ~ time + time:cond + (dummy(time, "post1w")||participant), 
#            data = us_data[us_data$group == "experiment",])
# 
# summary(rePCA(m2))
# 
# anova(m1, m2)
# # Does not improve fit, simplify further, removing random slopes
 m3 <-  lmer(thickness ~ time + time:cond + (1|participant), 
             data = us_data[us_data$group == "experiment",])
# 
# anova(m1, m2, m3)
# # Lowest AIC and non-significant improvement in models containing more complicated 
# # random effects structures. m3 is the preferred model.
# plot(m3)


summary(m3) 
confint(m3)

# There is no evidence for between condition differences in muscle thickness changes. This group is 
# subsequently analyzed as a whole against the control group. 

# This was confirmed using brms 
# m3.brms <- brm(thickness ~ time + time:cond + (1|participant), 
#                 
#                 warmup = 500, # number of samples before sampling
#                 iter = 3000,  # number of mcmc iterations 
#                 cores = 4, # number of cores used in sampling
#                 chains = 4, # number of chains
#                 seed = 123, # seed for reproducible results
#                 
#                 
#                 data = us_data[us_data$group == "experiment",])
# 
# 
# 
# summary(m3.brms)
# 



# Between group comparisons #########################

us_data2 <- us_data %>%

  mutate(time = as.character(time), 
         time = if_else(time == "pre", "baseline", time),
         time_pp = if_else(time == "post1w", "post", time),
       # The de-training period get its own coefficient

       group = if_else(group == "control", "con", "int"),
       
       detrain = if_else(time == "post1w" & group == "int", "detrain", "train"),
       # The effect of de training will be added to the model --within-- the intervention group 
       detrain = factor(detrain, levels = c("train", "detrain")), 

       tx = paste(time_pp, group, detrain, sep = "_"), 
       tx = factor(tx, levels = c("baseline_con_train", 
                                  "post_con_train", 
                                  "baseline_int_train", 
                                  "post_int_train", 
                                  "post_int_detrain"))) %>%
  # Removes duplicated data used for paralell modeling
  print()


# Modeling with brms
m0.us <- brm(thickness ~  tx + (1|participant),
             
             warmup = 1000, # number of samples before sampling
             iter = 4000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 123, # seed for reproducible results
             save_all_pars = TRUE,
             data = us_data2)

m1.us <- brm(thickness ~  tx + (1|participant/leg),
               
               warmup = 1000, # number of samples before sampling
               iter = 4000,  # number of mcmc iterations 
               cores = 4, # number of cores used in sampling
               chains = 4, # number of chains
               seed = 123, # seed for reproducible results
             save_all_pars = TRUE,
               data = us_data2)

# Save model 
saveRDS(m1.us, "./data/derivedData/us-freq-bayes/m1us.RDS")
saveRDS(m0.us, "./data/derivedData/us-freq-bayes/m0us.RDS")

m1.us <- readRDS("./data/derivedData/us-freq-bayes/m1us.RDS")
m0.us <- readRDS("./data/derivedData/us-freq-bayes/m0us.RDS")

loo(m0.us, m1.us, moment_match = TRUE, reloo = TRUE)
# Comparing these models indicates that the nested model (leg within participant)
# is better than the non nested (only participant). This is similar to the frequentis
# model selection. 

# Check if the model is ok

summary(m1.us)
pp_check(m1.us)
plot(m1.us)

# All looks ok, m1 from brms is selected. 

# Equivalent model in lme4

m0.us_freq <- lmer(thickness ~  tx + (1|participant),
     data = us_data2)

m1.us_freq <- lmer(thickness ~  tx + (1|participant/leg),
                   data = us_data2)

anova(m0.us_freq, m1.us_freq)
# Favors inclusion of random intercept per leg. 

em <- emmeans(m2.us_freq, specs = ~ time_pp * group  + time_pp:group:detrain)

# Combine data for plotting 

# Emmeans can be used as the contrats specified are combinations of factors. 
# emmeans::emmeans is used to create a emmGrid object:

em <- emmeans(m1.us_freq, specs =  ~ tx)

# To create contrasts, each levels can be specified by a indicator vector. 

baseline_con_train  =  c(1, 0, 0, 0, 0)
post_con_train      =  c(0, 1, 0, 0, 0) 
baseline_int_train  =  c(0, 0, 1, 0, 0) 
post_int_train      =  c(0, 0, 0, 1, 0) 
post_int_detrain    =  c(0, 0, 0, 0, 1)

# these can be combined in a named list to produce the contrasts of interest. 


comp_us <- contrast(em, 
                    method = list("post_con"         = post_con_train -  baseline_con_train, 
                                  "post_int"         = post_int_train - baseline_int_train, 
                                  "post1w_int"       = post_int_detrain - baseline_int_train, 
                                  "inter:post_int"   = (post_int_train - baseline_int_train) -  (post_con_train - baseline_con_train), 
                                  "inter:post1w_int" = (post_int_detrain - baseline_int_train) -  (post_con_train - baseline_con_train))) %>%
                         confint() %>%
                         data.frame() %>%
                   print()
                   


# Combine contrats for brms model

summary(m1.us)

h.us <- hypothesis(m1.us, c("txpost_con_train = 0", 
                                "txpost_int_train = txbaseline_int_train", 
                                "txpost_int_detrain = txbaseline_int_train", 
                                "txpost_int_train - txbaseline_int_train = txpost_con_train", 
                                "txpost_int_detrain - txbaseline_int_train = txpost_con_train"))



comp_us_bayes <- data.frame(time_group = c("post_con",        
                                                     "post_int",        
                                                     "post1w_int",      
                                                     "inter:post_int",  
                                                     "inter:post1w_int")) %>%
  cbind(
    data.frame(h.us$hypothesis)
  ) %>%
  mutate(estimate = Estimate, lower.CL = CI.Lower, upper.CL = CI.Upper) %>%
  dplyr::select(time_group, estimate, lower.CL, upper.CL) %>%
  print()



# Save data for plotting 
saveRDS(comp_us_bayes, "./data/derivedData/us-freq-bayes/comp_us_bayes.RDS")
saveRDS(comp_us, "./data/derivedData/us-freq-bayes/comp_us_freq.RDS")




