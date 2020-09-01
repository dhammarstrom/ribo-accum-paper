###### Strength and US analysis in a Bayesian framework

# This file produces analyses of strength and ultrasound in a Bayesian framework
# Output is equivalent to what is produced from strength-freq.R


# Use minimal packages as RStudio seems to chrash when RAM running low
library(brms)
library(tidybayes)
library(readxl)
library(tidyverse)





# Between-group comparisons (intervention vs. control) are done in two separate models 
# where post values in intervention (post/post1w) are compared to controls (post).

# Data is saved/plotted as estimated mean changes within groups (emmeans package) and
# interaction effects. 

strength_long <- read_excel("./data/tr010_humac.xlsx") %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  mutate(time = if_else(timepoint %in% c("B1", "B2", "fam"), "baseline", timepoint), 
         time = if_else(time == "post_ctrl", "post", time),
         cond = if_else(cond == "ctrl_leg", "ctrl", cond)) %>%
  
  dplyr::select(participant, 
                sex, 
                time,
                leg, 
                group,
                cond, 
                isok = isokinetic_torque,
                isom = isometric_torque) %>%
  print()




# Between group comparisons Strength ##########################
# Since the control group does not have two post-measurements a fully crossed model 
# will be rank deficient. Instead of a time * group interaction model, a new factor is 
# used to represent combinations in the data frame that we want to model. 
# Post-hoc comparisons are used to test interactions (i.e. change in intervention group 
# vs. change in the control group). 


strength_long2 <- strength_long %>%
  # Time coefficients are set as baseline/post
  mutate(time_pp = if_else(time == "post1w", "post", time),
         # The de-training period get its own coefficient
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



m1.isok <- brm(isok ~  tx + (1|participant),
               
               warmup = 500, # number of samples before sampling
               iter = 3000,  # number of mcmc iterations 
               cores = 4, # number of cores used in sampling
               chains = 4, # number of chains
               seed = 123, # seed for reproducible results
               
               data = strength_long2)


## Model diagnostics -- Posterior predictive Check

# Plots data and predictive values
pp_check(m1.isok)


pp_check(m1.isok, type = "stat_grouped", stat = "mean", group = "tx")
pp_check(m1.isok, type = "stat_2d", stat = c("max", "min"))

pp_check(m1.isok, type = "ribbon_grouped",  group = "tx") 
           


m1.isom <- brm(isom ~  tx + (1|participant),
               
               warmup = 500, # number of samples before sampling
               iter = 3000,  # number of mcmc iterations 
               cores = 4, # number of cores used in sampling
               chains = 4, # number of chains
               seed = 123, # seed for reproducible results
               
               data = strength_long2)

summary(m1.isom)

pp_check(m1.isom)


pp_check(m1.isom, type = "stat_grouped", stat = "mean", group = "tx")
pp_check(m1.isok, type = "stat_2d", stat = c("max", "min"))

pp_check(m1.isok, type = "ribbon_grouped",  group = "tx") 


plot(m1)


# "post-hoc" tests would be needed to test interactions. This can be done using 
# the hypothesis function. Hypotheses are written as character vectors 
# to create a credible interval around the parameter/estimate of interest -- these can be used for 
# inference. 

#### This is where IM at ########################

h <- hypothesis(m1.isok, "txpost_int_train - txbaseline_int_train  = txpost_con_train")

h


# Model diagnostics

# the ggs function transforms the brms output into a longformat tibble, 
# that we can use to make different types of plots.

m1_long <- ggs(m1.isok)



ggplot(filter(m1_long, Parameter %in% c("b_Intercept",
                                        "b_txpost_con_train",
                                        "b_txbaseline_int_train",
                                        "b_txpost_int_train",
                                        "b_txpost_int_detrain")),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line() +
  geom_vline(xintercept = 1000)+
  facet_grid(Parameter ~ . ,
             scale  = 'free_y',
             switch = 'y')+
  labs(title = "Caterpillar Plots", 
       col   = "Chains")






