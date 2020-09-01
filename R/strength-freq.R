####### Strength frequentist solution #########################

# This file produces strength analyses in a frequentis framework. 


# Statistical treatment of training effects (strength and muscle thickness) #
# Mixed effects models are fitted within experimental to evaluate condition effects
# and between intervention and control to evaluate effects of training. 
# Random effects are reduced from maximal models to the most parsimonious 
# based on lme4::rePCA and LRT tests of random effects.

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



# Between condition differences in strength increases ######

# Isokinetic data 


m1 <- lmer(isok ~ time + time:cond + (time|participant), 
           data = strength_long[strength_long$group == "int",])

summary(rePCA(m1))
# Singular, all variation is explained at the intercept, removing correlation reducing slopes

m2 <- lmer(isok ~ time + time:cond + (dummy(time, "post1w")||participant), 
           data = strength_long[strength_long$group == "int",])

summary(rePCA(m2))

anova(m1, m2)
# Does not improve fit, simplify further, removing random slopes
m3 <- lmer(isok ~  time + time:cond + (1|participant), 
           data = strength_long[strength_long$group == "int",])

anova(m1, m2, m3)
# Lowest AIC and non-significant improvement in models containing more complicated 
# random effects structures. m3 is the preferred model.
plot(m3)
confint(m3) 
# At least 1/3 of the 95% CI overlaps between conditions, little evidence for any between condition effects 
# in isokinetic strength. 

# Isometric data 

m1 <- lmer(log(isom) ~ time + time:cond + (time|participant), 
           data = strength_long[strength_long$group == "int",])

summary(rePCA(m1))
summary(m1)
# Singular, all variation is explained at the intercept, removing correlation reducing slopes

m2 <- lmer(log(isom) ~ time + time:cond + (dummy(time, "post1w")||participant), 
           data = strength_long[strength_long$group == "int",])

summary(rePCA(m2))

anova(m1, m2)
# Does not improve fit, simplify further, removing random slopes
m3 <- lmer(log(isom) ~  time + time:cond + (1|participant), 
           data = strength_long[strength_long$group == "int",])

anova(m1, m2, m3)
# Adding a random slope for time, uncorrelated to the intercept and only covering post1w improved the 
# fit. m2 for isom is the preferred model. 

#  confint(m2) 
# summary(m2)

# Symmetrical confidence intervals gives little evidence for condition effects on strength changes. 


# Between group comparisons Strength ##########################
# Since the control group does not have two post-measurements a fully crossed model 
# will be rank deficient. Insted of a time * group interaction model, a new factor is 
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


# Isokinetic models

m1 <- lmer(isok ~  tx + (1 + tx|participant), 
           data = strength_long2)

summary(rePCA(m1)) 
# the fit is not singular, but most of the variance is contained in the 
# intercept, simplify by removing correlation

m2 <- lmer(isok ~  tx + (dummy(tx, c("post_int_train", 
                                     "post_con_train"))||participant), 
           data = strength_long2)


summary(rePCA(m2)) # Sligthly more partition in the simpler model
# Fitting a simpler model

summary(m2)


m3 <- lmer(isok ~  tx + (1|participant), 
           data = strength_long2)

# The model is rank deficient (not all combinations of coefficients exists in the data)
# it is OK.

# Check if random effects structure matters
anova(m1, m2, m3) # No, the simplest model is preferred
plot(m3)

summary(m3)


## Isometric models 

m1.isom <- lmer(isom ~  tx + (1 + tx|participant), 
                data = strength_long2)

summary(rePCA(m1.isom)) 
plot(m1.isom)
# Not a singular fit, seeing if simplification may improve the model.

m2.isom <-  lmer(isom ~  tx + (dummy(tx, c("post_int_train", 
                                           "post_con_train"))||participant), 
                 data = strength_long2)


summary(rePCA(m2.isom)) # Some of the variance is contained in the uncorrelated random slope
# Fitting a simpler model
m3.isom <- lmer(isom ~  tx + (1|participant), 
                data = strength_long2)



# Check if random effects structure matters to model fit
anova(m1.isom, m2.isom, m3.isom) 

# The simplest model is preferred

plot(m3.isom)

summary(m3.isom)
# The random effects structure does not change the inference, but the simplest model has the lowest AIC
# and no improvement in fit with more elaborate models. The simplest model is preferred.



# Combine data for plotting 

# Emmeans can be used as the contrats specified are combinations of factors. 
# emmeans::emmeans is used to create a emmGrid object:

em.isok <- emmeans(m3, specs =  ~ tx)
em.isom <- emmeans(m3.isom, specs = ~ tx)


# The object contains no contrats, but simple means can be shown by calling the 
# object.  
em.isok

# To create contrats, each levels can be specified by a indicator vector. 

baseline_con_train  =  c(1, 0, 0, 0, 0)
post_con_train      =  c(0, 1, 0, 0, 0) 
baseline_int_train  =  c(0, 0, 1, 0, 0) 
post_int_train      =  c(0, 0, 0, 1, 0) 
post_int_detrain    =  c(0, 0, 0, 0, 1)

# these can be combined in a named list to produce the contrasts of interest. 


comp_strength <- rbind(contrast(em.isok, 
                                method = list("post_con"         = post_con_train -  baseline_con_train, 
                                              "post_int"         = post_int_train - baseline_int_train, 
                                              "post1w_int"       = post_int_detrain - baseline_int_train, 
                                              "inter:post_int"   = (post_int_train - baseline_int_train) -  (post_con_train - baseline_con_train), 
                                              "inter:post1w_int" = (post_int_detrain - baseline_int_train) -  (post_con_train - baseline_con_train))) %>%
                         confint() %>%
                         data.frame() %>%
                         mutate(type = "isok"), 
                       contrast(em.isom, 
                                method = list("post_con"         = post_con_train -  baseline_con_train, 
                                              "post_int"         = post_int_train - baseline_int_train, 
                                              "post1w_int"       = post_int_detrain - baseline_int_train, 
                                              "inter:post_int"   = (post_int_train - baseline_int_train) -  (post_con_train - baseline_con_train), 
                                              "inter:post1w_int" = (post_int_detrain - baseline_int_train) -  (post_con_train - baseline_con_train))) %>%
                         confint() %>%
                         data.frame() %>%
                         mutate(type = "isom")) %>%
  mutate(time_group = contrast) %>%
  print()



saveRDS(comp_strength, "./data/derivedData/strength-freq/comp_strength.RDS")






