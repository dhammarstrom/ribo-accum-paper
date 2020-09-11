### Western analaysis 

library(tidyverse)
library(brms)
library(bayesplot)
library(lme4)

western_data <- readRDS("./data/derivedData/western-compile/western_data.RDS")





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

  print()
  



s6.m1 <- brm(bf(log.expr ~ tx + (1|participant)),

          warmup = 1000, # number of samples before sampling
          iter = 4000,  # number of mcmc iterations 
          cores = 4, # number of cores used in sampling
          chains = 4, # number of chains
          seed = 5, # seed for reproducible results
      #    control = list(adapt_delta = 0.95), 
          data = west.dat[west.dat$target == "t-s6",])


ubf.m1 <- brm(bf(log.expr ~ tx + (1|participant)),
             
             warmup = 1000, # number of samples before sampling
             iter = 4000,  # number of mcmc iterations 
             cores = 4, # number of cores used in sampling
             chains = 4, # number of chains
             seed = 5, # seed for reproducible results
             #    control = list(adapt_delta = 0.95), 
             data = west.dat[west.dat$target == "t-UBF",])




# Save model 
saveRDS(s6.m1, "./data/derivedData/western-analysis/s6_tx_m1.RDS")
saveRDS(ubf.m1, "./data/derivedData/western-analysis/ubf_tx_m1.RDS")


s6.m1 <- readRDS( "./data/derivedData/western-analysis/s6_tx_m1.RDS")
loo(s6.m1)

# Model checks

# predictive posterior checks
pp_check(s6.m1, type = "stat", stat = "median")
pp_check(s6.m1, type = "ecdf_overlay")

pp_check(ubf.m1, type = "stat", stat = "median")
pp_check(ubf.m1, type = "ecdf_overlay")

# Leave one out diagnostics
loo(s6.m1)
loo(ubf.m1)


summary(s6.m1)
# combine factors for hypotheses

tx_contrasts <-  c("txS1_int_train -     txpre_int_train =      txS1_con_train",
                   # Interaction effects post-training
                   "txpost_int_train -     txpre_int_train =      txpost_con_train",
                   # Interaction effects detraining (post1w)
                   "txpost_int_detrain -     txpre_int_train =      txpost_con_train",
                   # Whitin group fold change
                   "txS1_int_train =     txpre_int_train",     
                   "txpost_int_train =     txpre_int_train",     
                   "txpost_int_detrain =     txpre_int_train",     
                   "txS1_con_train = 0",
                   "txpost_con_train = 0")



h.s6 <- hypothesis(s6.m1, tx_contrasts)
h.ubf <- hypothesis(ubf.m1, tx_contrasts)

# Create a data frame for the results 
tx_results_western <- data.frame(target = rep(c("rpS6","UBF"), each = 8), 
                               comparison = rep(c("inter:S1",
                                              "inter:post", 
                                              "inter:post1w", 
                                      
                                              "int_S1", 
                                              "int_post", 
                                              "int_post1w", 
                                      
                                              "con_S1",
                                              "con_post"), 2)) %>%
  cbind(data.frame(rbind(h.s6$hypothesis[,-1], 
                   h.ubf$hypothesis[,-1]) )) %>%
  print()


#### Fold change plot rpS6 and UBF 


tx_results_western %>%
  filter(comparison %in% c("int_S1", "int_post", "int_post1w", "con_S1", "con_post")) %>%
  mutate(Estimate = exp(Estimate), 
         CI.Lower= exp(CI.Lower), 
         CI.Upper = exp(CI.Upper)) %>%
  
  separate(comparison, into = c("group", "time")) %>%
  
  mutate(detrain = if_else(time == "post1w", "detrain", "train"),
         time = if_else(time == "post1w", "post", time), 
         time = factor(time, levels = c("S1", "post"), 
                       labels = c("Session 1", "Post-\ntraining"))) %>%
  
  ggplot(aes(time, Estimate, color = paste(group, detrain))) + 
  
  geom_errorbar(aes(ymin = CI.Lower, ymax = CI.Upper), 
                position = position_dodge(width = 0.2), 
                width = 0.1) +
  geom_point(position = position_dodge(width = 0.2)) + 
  
  facet_grid(target ~ .)


### Comparison plot 

tx_results_western %>%
  filter(!(comparison %in% c("int_S1", "int_post", "int_post1w", "con_S1", "con_post"))) %>%
  mutate(Estimate = exp(Estimate), 
         CI.Lower= exp(CI.Lower), 
         CI.Upper = exp(CI.Upper)) %>%
  mutate(comparison =  gsub("inter:", "", comparison), 
         comparison = factor(comparison, levels = c("S1", "post", "post1w"), 
                             labels = c("Session 1", "Post-training", 
                                        "Post-training\n+De-training")), 
         comparison = fct_rev(comparison)) %>%
  ggplot(aes(comparison, Estimate)) + 
  geom_errorbar(aes(ymin = CI.Lower, ymax = CI.Upper), 
                width = 0.1) +
    geom_point() + 
  facet_grid(.~ target) + coord_flip()




########## Continous data between conditions ##########

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
       time = factor(time, levels = c("S0", 
                                      "S1", 
                                      "S4", 
                                      "S5", 
                                      "S8", 
                                      "S9", 
                                      "S12")), 
       time.c = as.numeric(gsub("S", "", time))) %>%
  
  group_by(target, participant, time.c, cond, time) %>%
  summarise(expression = mean(expression)) %>%
  mutate(log.expr = log(expression)) %>%
  filter(!is.na(expression)) %>%
  
  print()


library(mgcv)

m1 <- gam(log.expr ~  cond + s(time.c, k = 7, by = cond), 
          random = list(participant = ~ 1), 
          data = west_cont)

m1 <- lmer(log.expr ~  cond * time + (1|participant), 
          data = west_cont[west_cont$target== "t-s6",])


summary(m1)
plot(m1)

west_cont %>%

  ggplot(aes(time.c,log.expr,
             color = cond,
             group = paste(participant, cond))) +
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap( ~ target)




plot(m1)


plot(m1)


west_conds <- western_data %>%
  filter(cond != "ctrl_leg", 
         time %in% c("S0", 
                    
                     "S12", 
                     "post1w")) %>%
  

  mutate(participant = factor(participant), 
         cond = factor(cond), 
         target = factor(target), 
         time = factor(time, levels = c("S12", "S0",  "post1w"))) %>%
  
  group_by(target, participant, cond, time) %>%
  summarise(expression = mean(expression)) %>%
  mutate(log.expr = log(expression)) %>%
  filter(!is.na(expression)) %>%
  
  print()


m1 <- lmer(log.expr ~ time + time:cond + (1|participant), 
     data = west_conds[west_conds$target == "t-UBF",])

summary(m1)





#### Graph of trends for all targets ### 
western_data %>%
  mutate(time = factor(time, levels = c("S0", "S1","S1c", "S4", "S5", "S6", 
                                        "S8", "S9", "S12", "post1w", "postctrl"))) %>%
  group_by(time, cond, target) %>%

summarise(expression = mean(log(expression), na.rm = TRUE)) %>%
  ggplot(aes(time, expression, group = cond, color = cond)) + geom_point() + geom_line() + 
  facet_wrap(~target, scales = "free") + 
  ylab("Mean log(expression)") +
  xlab("Time-point")


### Linear trends with session as continuous variable ###
western_data %>%
  mutate(time.c = gsub("S", "", time), 
         detrain = if_else(time.c == "post1w", "detrain", "train"),
         time.c = if_else(time.c == "post1w", "12", time.c),
         time.c = if_else(time.c == "postctrl", "12", time.c),
         time.c = gsub("c", "", time.c), 
         time.c = as.numeric(time.c), 
         time = factor(time, levels = c("S0", "S1","S1c", "S4", "S5", "S8", "S9", "S12", "post1w", "postctrl"))) %>%
  mutate(log.expr = log(expression)) %>%
  group_by(participant, leg, time.c, cond, detrain, target) %>%
  filter(detrain == "train") %>%
  summarise(log.expr = mean(log.expr, na.rm = TRUE),
            expression = mean(expression, na.rm = TRUE))%>%
  ggplot(aes(time.c, expression, color = cond)) + geom_smooth(method = "lm") + 
  facet_wrap(~target, scales = "free") +
  ylab("log(expression)") +
  xlab("Time-point")

### Effects of training-detraining in experimental group ###
western_data %>%
  filter(cond != "ctrl_leg", time %in% c("S0", "S12", "post1w")) %>%
  mutate(time = factor(time, levels = c("S0", "S12", "post1w"))) %>%
  group_by(participant, cond, time, target) %>%
  summarise(expression = mean(expression)) %>%
  ggplot(aes(time, log(expression), fill = cond)) + geom_boxplot() + 
  facet_wrap(~target, scales = "free") +
  ylab("log(expression)") +
  xlab("Time-point")


### Preliminary models ###  

# Treating time as a factor variable to assess effects on 
# training effect and de-training 

## t-UBF model  
western_data %>%
  filter(cond != "ctrl_leg", time %in% c("S0", "S12", "post1w")) %>%
  mutate(time = factor(time, levels = c("S12","S0",  "post1w"))) %>%
  filter(target == "t-UBF") %>%
  group_by(participant, leg, cond, time) %>%
  summarise(expression = mean(expression, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(log.expr = log(expression),
         leg = paste0(participant, leg)) %>%
  lmer(log.expr ~ time * cond + (1|participant/leg), data = .) %>%
  summary()

## t-s6 model
western_data %>%
  filter(cond != "ctrl_leg", time %in% c("S0", "S12", "post1w")) %>%
  mutate(time = factor(time, levels = c("S12","S0",  "post1w"))) %>%
  filter(target == "t-s6") %>%
  group_by(participant, leg, cond, time) %>%
  summarise(expression = mean(expression, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(log.expr = log(expression),
         leg = paste0(participant, leg)) %>%
  lmer(log.expr ~ time +  time:cond + (1|participant) + (1|leg), data = .) %>%
  summary()

##### Data for continuous models #### 
cont.dat <-  western_data %>%
  mutate(time.c = gsub("S", "", time), 
         detrain = if_else(time.c == "post1w", "detrain", "train"),
         time.c = if_else(time.c == "post1w", "12", time.c),
         time.c = if_else(time.c == "postctrl", "12", time.c),
         time.c = gsub("c", "", time.c), 
         time.c = as.numeric(time.c), 
         time = factor(time, levels = c("S0", "S1","S1c", "S4", "S5", "S8", "S9", "S12", "post1w", "postctrl")),
         detrain = factor(detrain, levels = c("train", "detrain")),
         repeated.biopsy = if_else(time %in% c("S1", "S5", "S9"), "repeated", "notrepeated")) %>%
  filter(detrain == "train", 
         cond != "ctrl_leg" ) %>%
  mutate(log.expr = log(expression),
         sample_gel = paste0(leg, "_g", gel), 
         cond = factor(cond),
         repeated.biopsy = factor(repeated.biopsy),
         sample_total = paste0(participant, leg, time),
         leg = paste0(participant, leg),
         # Use with piece wise models
         # piece-wise growth models or segmented analysis alternative 1
         # Interpretation alternative 1: each slope represents the actual 
         # regression line.
         time1 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), time.c, 8),
         time2 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8),
         time3 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8)) %>%
  # Alternative 2: Each slope except for slope "time1" represents 
  # deviation from the t1 slope
  # time1 = time.c,
  # time2 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8),
  # time3 = if_else(time %in% c("S0", "S1", "S4", "S5", "S8"), 0, time.c - 8)) %>%
  
  print()


cont.dat %>%
  ggplot(aes(time.c, log.expr, color = cond)) + geom_point() +
  geom_smooth() + facet_wrap(~target, scales = "free")




### Piece-wise models may be fitted with lmer or nlme

# Example UBF 

# Comments:
# There seems to be a repeated-biosy effect in samples from 48-h after biopsy. This is 
# handeled through modeling a repeated biopsy effect. This effect is significant in 
# UBF. There is also a deviation in the second segment in the var group compared to 
# const. 

# Model:

library(mgcv)

mx <- gam(log.expr ~ cond + s(time.c, k = 7, by = cond), 
          random = list(participant = ~ 1), 
          data = cont.dat[cont.dat$target == "t-s6",])


plot(mx)





m1 <- lmer(log.expr ~ time.c * cond +
             (1|participant), 
           data = cont.dat[cont.dat$target == "t-s6",])


summary(m1)


m2 <- lmer(log.expr ~ poly(time.c, 2, raw = TRUE) * cond + 
             (1|participant) +
             (1|leg) +
             (1|sample_gel), 
           data = cont.dat[cont.dat$target == "t-UBF",])


anova(m1, m2)


summary(m2)

library(effects)

effects(m1, "time:cond")

position = position_dodge(width = 0.2)
as.data.frame(allEffects(m1)) %>%
  data.frame() %>%
  mutate(time.cond.time = factor(time.cond.time, levels = c("S0", "S1", "S4", "S5", "S8", "S9", "S12"))) %>%
  ggplot(aes(time.cond.time, time.cond.fit, color = time.cond.cond)) + 
  geom_errorbar(aes(ymin = time.cond.lower, ymax = time.cond.upper), 
                width = 0.2,
                position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2))

confint(m1)
m2 <- lmer(log.expr ~ time1 + time2 + cond + time1:cond + time2:cond +  repeated.biopsy + 
             (1 +time1|participant) + 
             (1 + time2|sample_gel), 
           data = cont.dat[cont.dat$target == "t-UBF",])

m2 <- lmer(log.expr ~ time1 + time2 + cond + time1:cond + time2:cond +  repeated.biopsy + 
             (1 +time1|participant) + 
             (1|sample_gel), 
           data = cont.dat[cont.dat$target == "dgkz",])

anova(m2, m1) 
# No additional benefit of leg random intercept (nested within participant) 


# Summary of model   
summary(m1) 



#### DGKZ #### 

#### Graph of trends ### 
cont.dat %>%
  filter(target == "dgkz") %>%
  ggplot(aes(time.c, log(expression), color = cond)) + geom_point() +
  geom_smooth()
print()





western_data %>%
  mutate(time = factor(time, levels = c("S0", "S1","S1c", "S4", "S5", "S6", 
                                        "S8", "S9", "S12", "post1w", "postctrl"))) %>%
  filter(target == "dgkz") %>%
  group_by(participant, time, cond, target, round) %>%
  summarise(expression = mean(log(expression), na.rm = TRUE)) %>%
  
  
  ggplot(aes(time, expression, group = paste0(cond, participant, round), 
             color = round,
             shape = cond)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~ participant, scales = "free")

ylab("Mean log(expression)") +
  xlab("Time-point")


western_data %>% 
  filter(participant == "P4", target == "dgkz") %>%
  print()


western_data %>%
  mutate(time = factor(time, levels = c("S0", "S1","S1c", "S4", "S5", "S6", 
                                        "S8", "S9", "S12", "post1w", "postctrl"))) %>%
  filter(time  %in% c("S0", "S1", "S4", "S12")) %>%
  filter(target == "dgkz") %>%
  lmer(log(expression) ~ time  * cond + (1|participant),
       data = .) %>%
  summary()



group_by(time, participant) %>%
  summarise(expression = mean(log(expression), na.rm = TRUE)) %>%
  ggplot(aes(time, expression, group = participant)) + geom_line()



### Extract effects from models with effects package ###

# Scraps from here ---->

library(effects)

eff.m1 <- allEffects(m1, xlevels=list(x1 = 5, x2=5))
allEffects(m1) %>%
  as.data.frame() %>% 
  print()

eff.m1$`time1:cond`


ggplot(aes(time2, fit, group = cond)) + 
  geom_line(aes(color = cond)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  scale_x_continuous(breaks = seq(0:12))


