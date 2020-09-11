##### qpcr analysis using brms #########################

# qpcr analysis in a Bayesian framework.


# Use minimal packages as RStudio seems to chrash when RAM running low
library(brms)
library(tidybayes)
library(readxl)
library(tidyverse)




# qdat is the compiled data frame containing sample information and 
# qpcr estimates. Created in qpcr-compile.R. Use this for downstream analyses

qdat <- readRDS("./data/derivedData/qpcr/qpcr_compiled.RDS")



# Remove bad reactions prior to modeling. 

# Explorative plots

# Remove bad reaction prior to modelling. 
# Bad reactions are no amplification or estimatyed to cq < 5



qdat %>%
  ggplot(aes(cq, color = target)) + geom_histogram() + 
  facet_wrap(~ target)

qdat %>%
  group_by(target) %>%

  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 2 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 2 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  
  ggplot(aes(paste(participant, leg, time), cq, 
             color = target, 
             shape = outlier, 
             alpha = outlier)) + 
  geom_point() + theme_minimal()  + scale_alpha_manual(values = c(0.2, 0.8))

# 3 sd's away from target mean seems to capture the worst reactions.. 

  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 5 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 5 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  
  ggplot(aes(paste(participant, leg, time), cq, color = target, shape = outlier)) + 
  geom_point() + theme_minimal() 

# 5 sd's away from target mean seems to capture the worst reactions.. 




### Lambda normalization ####


nf <- qdat %>% group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 2 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 2 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  filter(target == "Lambda KIT") %>%
 filter(outlier == "in") %>% # Removes bad reactions

nf <- qdat %>%
  filter(target == "Lambda KIT") %>%
  filter(cq > 5) %>% # Removes bad reactions

  mutate(nf.w = (eff ^ -cq) * tissue_weight) %>%
  dplyr::select(participant, leg, time, cdna, nf.w) %>%
  group_by(participant, leg, time, cdna) %>%
  summarise(nf.w = mean(nf.w, na.rm = TRUE)) %>%

  
  ungroup() %>%
  # Scale the factor
  mutate(nf.w = nf.w/max(nf.w)) %>%

  print()
  

#### rRNA per tissue weight analysis #### 

qdat.rrna  <- qdat %>%
  group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 2 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 2 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
 filter(outlier == "in") %>%

  ungroup() %>%
  print()





#### rRNA per tissue weight analysis #### 

# THIW IS WHERE IM AT ###############

qdat.rrna  <- qdat %>%
  group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 5 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 5 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  filter(outlier == "in") %>%

  group_by(participant, leg, time, sex, cond, cdna, target) %>%
  summarise(cq = mean(cq, na.rm = TRUE), 
            eff = mean(eff, na.rm = TRUE), 
            tissue_weight = mean(tissue_weight, na.rm = TRUE)) %>%
  
  filter(target %in% c("rRNA18S F2R2",      "rRNA28S F2R2",     
                       "rRNA5.8S F2R2",       
                       "rRNA45S F5R5",      "rRNA45SITS F12R12", 

                       "rRNA47S F1R1",      "rRNA5S F3R3", 
                       "UBTF F4R4", "UBTF F6R6", 
                       "rpS6 F2R2")) %>%

                       "rRNA47S F1R1",      "rRNA5S F3R3")) %>%

  mutate(Ra = -cq * log(eff), # Relative abundance
         counts = as.integer(round(eff ^ (39-cq), 0)), # counts, ref Matz et al.
         technical = paste0(participant, leg, time, cdna), 
         biological = paste0(participant, leg, time), 
         leg.unique = paste0(participant, leg)) %>%
  mutate(time.c = gsub("S", "", time), 
         detrain = if_else(time.c == "post1w", "detrain", "train"),
         time.c = if_else(time.c == "post1w", "12", time.c),
         time.c = if_else(time.c == "postctrl", "12", time.c),
         time.c = gsub("c", "", time.c), 
         time.c = as.numeric(time.c), 
         time = factor(time, levels = c("S0", "S1","S1c", 
                                        "S4", "S5", "S8", 
                                        "S9", "S12", "post1w", "postctrl"))) %>%
  #print()
  ungroup() %>%
  inner_join(nf, by = c("participant", "leg", "time", "cdna")) %>%
  mutate(Ra.nf = (eff^-cq)/nf.w) %>%
  print()



## Control vs. intervention group #
# Due to rank deficiency in comparing ctrl vs. interv a new factor variable is formed 
# and comparisons are made post-hoc. This is similar to what is done in other analyses.


qdat.ctrl <- qdat.rrna %>%
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
                                    "post_int_detrain"))) %>%
  ungroup() %>%
         mutate(nf.w = nf.w / mean(nf.w, na.rm = TRUE))  %>%
  print()





m1 <- brm(bf(counts ~ 0 + target + target:tx + (1|participant) + (1|technical) + 
               offset(nf.w), 
             shape ~ target),
          family = negbinomial(),
          warmup = 1000, # number of samples before sampling
          iter = 8000,  # number of mcmc iterations 
          cores = 6, # number of cores used in sampling
          chains = 6, # number of chains
          seed = 123, # seed for reproducible results
          control = list(adapt_delta = 0.95), 
          data = qdat.ctrl)

# Save model 
saveRDS(m1, "./data/derivedData/qpcr-analysis-bayes/ctrl_vs_int.RDS")

m1 <- readRDS("./data/derivedData/qpcr-analysis-bayes/ctrl_vs_int.RDS")


pp_check(m1, type = "stat", stat = "median")
summary(m1)

m1 <- brm(bf(counts ~ 0 + target + target:tx + (1|participant) + (1|technical) + offset(nf.w), 
             shape ~ target),
          family = negbinomial(),
          warmup = 500, # number of samples before sampling
          iter = 4000,  # number of mcmc iterations 
          cores = 4, # number of cores used in sampling
          chains = 4, # number of chains
          seed = 123, # seed for reproducible results

          data = qdat.ctrl)




# combine factors for hypotheses


h <- hypothesis(m1, 
        # Interaction effects S1
        c("targetrRNA18SF2R2:txS1_int_train -     targetrRNA18SF2R2:txpre_int_train =      targetrRNA18SF2R2:txS1_con_train",
          "targetrRNA45SF5R5:txS1_int_train -     targetrRNA45SF5R5:txpre_int_train =      targetrRNA45SF5R5:txS1_con_train",
          "targetrRNA47SF1R1:txS1_int_train -     targetrRNA47SF1R1:txpre_int_train =      targetrRNA47SF1R1:txS1_con_train",
          "targetrRNA5SF3R3:txS1_int_train -      targetrRNA5SF3R3:txpre_int_train =       targetrRNA5SF3R3:txS1_con_train",
          "targetrRNA28SF2R2:txS1_int_train -     targetrRNA28SF2R2:txpre_int_train =      targetrRNA28SF2R2:txS1_con_train",
          "targetrRNA45SITSF12R12:txS1_int_train -targetrRNA45SITSF12R12:txpre_int_train = targetrRNA45SITSF12R12:txS1_con_train", 
          "targetrRNA5.8SF2R2:txS1_int_train -    targetrRNA5.8SF2R2:txpre_int_train =     targetrRNA5.8SF2R2:txS1_con_train", 

          
          "targetUBTFF4R4:txS1_int_train -     targetUBTFF4R4:txpre_int_train =      targetUBTFF4R4:txS1_con_train",
          "targetUBTFF6R6:txS1_int_train -targetUBTFF6R6:txpre_int_train = targetUBTFF6R6:txS1_con_train", 
          "targetrpS6F2R2:txS1_int_train -    targetrpS6F2R2:txpre_int_train =     targetrpS6F2R2:txS1_con_train", 
          

          # Interaction effects post-training
          "targetrRNA18SF2R2:txpost_int_train -     targetrRNA18SF2R2:txpre_int_train =      targetrRNA18SF2R2:txpost_con_train",
          "targetrRNA45SF5R5:txpost_int_train -     targetrRNA45SF5R5:txpre_int_train =      targetrRNA45SF5R5:txpost_con_train",
          "targetrRNA47SF1R1:txpost_int_train -     targetrRNA47SF1R1:txpre_int_train =      targetrRNA47SF1R1:txpost_con_train",
          "targetrRNA5SF3R3:txpost_int_train -      targetrRNA5SF3R3:txpre_int_train =       targetrRNA5SF3R3:txpost_con_train",
          "targetrRNA28SF2R2:txpost_int_train -     targetrRNA28SF2R2:txpre_int_train =      targetrRNA28SF2R2:txpost_con_train",
          "targetrRNA45SITSF12R12:txpost_int_train -targetrRNA45SITSF12R12:txpre_int_train = targetrRNA45SITSF12R12:txpost_con_train", 
          "targetrRNA5.8SF2R2:txpost_int_train -    targetrRNA5.8SF2R2:txpre_int_train =     targetrRNA5.8SF2R2:txpost_con_train", 

          
          
          "targetUBTFF4R4:txpost_int_train -     targetUBTFF4R4:txpre_int_train =      targetUBTFF4R4:txpost_con_train",
          "targetUBTFF6R6:txpost_int_train -targetUBTFF6R6:txpre_int_train = targetUBTFF6R6:txpost_con_train", 
          "targetrpS6F2R2:txpost_int_train -    targetrpS6F2R2:txpre_int_train =     targetrpS6F2R2:txpost_con_train", 
          
          
          

          # Interaction effects detraining (post1w)
          "targetrRNA18SF2R2:txpost_int_detrain -     targetrRNA18SF2R2:txpre_int_train =      targetrRNA18SF2R2:txpost_con_train",
          "targetrRNA45SF5R5:txpost_int_detrain -     targetrRNA45SF5R5:txpre_int_train =      targetrRNA45SF5R5:txpost_con_train",
          "targetrRNA47SF1R1:txpost_int_detrain -     targetrRNA47SF1R1:txpre_int_train =      targetrRNA47SF1R1:txpost_con_train",
          "targetrRNA5SF3R3:txpost_int_detrain -      targetrRNA5SF3R3:txpre_int_train =       targetrRNA5SF3R3:txpost_con_train",
          "targetrRNA28SF2R2:txpost_int_detrain -     targetrRNA28SF2R2:txpre_int_train =      targetrRNA28SF2R2:txpost_con_train",
          "targetrRNA45SITSF12R12:txpost_int_detrain -targetrRNA45SITSF12R12:txpre_int_train = targetrRNA45SITSF12R12:txpost_con_train",
          "targetrRNA5.8SF2R2:txpost_int_detrain -    targetrRNA5.8SF2R2:txpre_int_train =     targetrRNA5.8SF2R2:txpost_con_train", 
          

          "targetUBTFF4R4:txpost_int_detrain -     targetUBTFF4R4:txpre_int_train =         targetUBTFF4R4:txpost_con_train",
          "targetUBTFF6R6:txpost_int_detrain -     targetUBTFF6R6:txpre_int_train =         targetUBTFF6R6:txpost_con_train",
          "targetrpS6F2R2:txpost_int_detrain -     targetrpS6F2R2:txpre_int_train =         targetrpS6F2R2:txpost_con_train", 
          
          

          # Whitin group fold change
          "targetrRNA18SF2R2:txpost_int_detrain =     targetrRNA18SF2R2:txpre_int_train",     
          "targetrRNA45SF5R5:txpost_int_detrain =     targetrRNA45SF5R5:txpre_int_train",     
          "targetrRNA47SF1R1:txpost_int_detrain =     targetrRNA47SF1R1:txpre_int_train",     
          "targetrRNA5SF3R3:txpost_int_detrain =      targetrRNA5SF3R3:txpre_int_train",      
          "targetrRNA28SF2R2:txpost_int_detrain =     targetrRNA28SF2R2:txpre_int_train",    
          "targetrRNA45SITSF12R12:txpost_int_detrain = targetrRNA45SITSF12R12:txpre_int_train",
          "targetrRNA5.8SF2R2:txpost_int_detrain =    targetrRNA5.8SF2R2:txpre_int_train", 
          

          "targetUBTFF4R4:txpost_int_detrain =     targetUBTFF4R4:txpre_int_train",    
          "targetUBTFF6R6:txpost_int_detrain =     targetUBTFF6R6:txpre_int_train",
          "targetrpS6F2R2:txpost_int_detrain =     targetrpS6F2R2:txpre_int_train", 
          
          
          

          "targetrRNA18SF2R2:txpost_int_train =     targetrRNA18SF2R2:txpre_int_train",     
          "targetrRNA45SF5R5:txpost_int_train =     targetrRNA45SF5R5:txpre_int_train",     
          "targetrRNA47SF1R1:txpost_int_train =     targetrRNA47SF1R1:txpre_int_train",     
          "targetrRNA5SF3R3:txpost_int_train =      targetrRNA5SF3R3:txpre_int_train",      
          "targetrRNA28SF2R2:txpost_int_train =     targetrRNA28SF2R2:txpre_int_train",    
          "targetrRNA45SITSF12R12:txpost_int_train = targetrRNA45SITSF12R12:txpre_int_train",
          "targetrRNA5.8SF2R2:txpost_int_train =    targetrRNA5.8SF2R2:txpre_int_train", 
          

          "targetUBTFF4R4:txpost_int_train =      targetUBTFF4R4:txpre_int_train",    
          "targetUBTFF6R6:txpost_int_train =      targetUBTFF6R6:txpre_int_train",
          "targetrpS6F2R2:txpost_int_train =      targetrpS6F2R2:txpre_int_train", 
          

          "targetrRNA18SF2R2:txS1_int_train =     targetrRNA18SF2R2:txpre_int_train",     
          "targetrRNA45SF5R5:txS1_int_train =     targetrRNA45SF5R5:txpre_int_train",     
          "targetrRNA47SF1R1:txS1_int_train =     targetrRNA47SF1R1:txpre_int_train",     
          "targetrRNA5SF3R3:txS1_int_train =      targetrRNA5SF3R3:txpre_int_train",      
          "targetrRNA28SF2R2:txS1_int_train =     targetrRNA28SF2R2:txpre_int_train",    
          "targetrRNA45SITSF12R12:txS1_int_train = targetrRNA45SITSF12R12:txpre_int_train",
          "targetrRNA5.8SF2R2:txS1_int_train =    targetrRNA5.8SF2R2:txpre_int_train", 

          "targetUBTFF4R4:txS1_int_train =    targetUBTFF4R4:txpre_int_train",    
          "targetUBTFF6R6:txS1_int_train =    targetUBTFF6R6:txpre_int_train",
          "targetrpS6F2R2:txS1_int_train =    targetrpS6F2R2:txpre_int_train", 
          
          
          

          "targetrRNA18SF2R2:txS1_con_train = 0",
          "targetrRNA45SF5R5:txS1_con_train = 0",
          "targetrRNA47SF1R1:txS1_con_train = 0",
          "targetrRNA5SF3R3:txS1_con_train = 0",
          "targetrRNA28SF2R2:txS1_con_train = 0",
          "targetrRNA45SITSF12R12:txS1_con_train = 0", 
          "targetrRNA5.8SF2R2:txS1_con_train = 0", 
          

          "targetUBTFF4R4:txS1_con_train = 0",
          "targetUBTFF6R6:txS1_con_train = 0", 
          "targetrpS6F2R2:txS1_con_train = 0", 
          
          

          "targetrRNA18SF2R2:txpost_con_train = 0",
          "targetrRNA45SF5R5:txpost_con_train = 0",
          "targetrRNA47SF1R1:txpost_con_train = 0",
          "targetrRNA5SF3R3:txpost_con_train = 0",
          "targetrRNA28SF2R2:txpost_con_train = 0",
           "targetrRNA45SITSF12R12:txpost_con_train = 0",

           "targetrRNA5.8SF2R2:txpost_con_train = 0", 
          
          "targetUBTFF4R4:txpost_con_train = 0",
          "targetUBTFF6R6:txpost_con_train = 0", 
          "targetrpS6F2R2:txpost_con_train = 0" 
          
          

           "targetrRNA5.8SF2R2:txpost_con_train = 0"

          ))
h

# Create a data frame for the results 
qpcr_res_int_con <- data.frame(target = rep(c("rRNA18SF2R2"      ,      
                                          "rRNA45SF5R5"      ,
                                          "rRNA47SF1R1"      ,
                                          "rRNA5SF3R3"       ,
                                          "rRNA28SF2R2"      ,    
                                          "rRNA45SITSF12R12" ,

                                          "rRNA5.8SF2R2", 
                                          "UBTFF4R4",
                                          "UBTFF6R6", 
                                          "rpS6F2R2" ), 8), 
                               comparison = c(rep("inter:S1", 10),
                                              rep("inter:post", 10),
                                              rep("inter:post1w", 10), 
                                              
                                              rep("int_post1w", 10), 
                                              rep("int_post", 10),
                                              rep("int_S1", 10),
                                              
                                              rep("con_S1", 10),
                                              rep("con_post", 10))) %>%

                                          "rRNA5.8SF2R2"), 8), 
                               comparison = c(rep("inter:S1", 7),
                                              rep("inter:post", 7),
                                              rep("inter:post1w", 7), 
                                              
                                              rep("int_S1", 7), 
                                              rep("int_post", 7),
                                              rep("int_post1w", 7),
                                              
                                              rep("con_S1", 7),
                                              rep("con_post", 7))) %>%

  cbind(data.frame(h$hypothesis)[,-1]) %>%
  print()



saveRDS(qpcr_res_int_con, "./data/derivedData/qpcr-analysis-bayes/qpcr_res_int_con.RDS")





# pre- to mature-rRNA analysis ###################

# The aim of the normalization procedure is to map pre-rRNA fluctuations per 
# mature rRNA to explore de novo transcription of rRNA. 

# Creating the normalization factor. 

nf <- qdat.rrna %>%
  filter(cond != "ctrl_leg") %>%
  filter(time != "post1w") %>%
  filter(target %in% c("rRNA28S F2R2",
                       "rRNA5.8S F2R2", 
                       "rRNA18S F2R2")) %>%
  group_by(participant, leg, time, cond, time.c, cdna) %>%
  summarise(nf = sum(counts, na.rm = TRUE)) %>%
  ungroup() %>%
  print()

  
nf %>%
  ggplot(aes(time.c, log(nf))) + geom_point() +
  facet_grid(leg ~  participant)


#### Model 47S/45S per mature rRNA ##############################





qdat_pre <- qdat.rrna %>%
  filter(cond != "ctrl_leg") %>%
  filter(time != "post1w") %>%
  filter(target %in% c("rRNA45S F5R5", 
                       "rRNA47S F1R1")) %>%
  


  inner_join(nf) %>%
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
        time.c = as.numeric(gsub("S", "", time)), 
        norm = Ra - log(nf)) %>% 
  filter(!(participant == "P4" & time == "S5" & cond == "const" & cdna == "cDNA1")) %>%

qdat_rrna2 <- qdat.rrna %>%
  filter(cond != "ctrl_leg") %>%
  dplyr::select(participant, leg, time, time.c, cond, technical, biological, target,Ra, counts) %>%
  
  filter(target %in% c("rRNA18S F2R2", 
                       "rRNA28S F2R2", 
                       "rRNA5.8S F2R2", 
                       "rRNA45S F5R5", 
                       "rRNA47S F1R1")) %>%
  
  mutate(type = if_else(target %in% c("rRNA18S F2R2", 
                                      "rRNA28S F2R2", 
                                      "rRNA5.8S F2R2"), 
                        "mature", "pre")) %>%

  filter(time != "post1w") %>%


  print()



qdat_pre %>%
  ggplot(aes(nf, norm, color = participant)) + geom_point() +
  geom_text(data = subset(qdat_pre, norm > -32), 
            aes(label = paste(participant, time, cond, cdna)), 
            hjust = 0)


qdat_pre %>%
  ggplot(aes(time.c, log(counts)-log(nf), color = cond)) + geom_point() +
  geom_smooth() + facet_grid(.~target)

qdat_pre %>%
  filter(time == "S0") %>%
  ggplot(aes(participant, cq, color = cond)) + 
  geom_point() + 
  facet_grid(.~ target)



library(nlme); library(emmeans); library(mgcv); library(ggeffects)



m <- lme(norm ~ target + time * cond, 
         random = list(participant = ~ 1, technical = ~ 1),
         weights = varIdent(form = ~ 1|target),
             data = qdat_pre)

plot(m)

summary(m)

m <- gam(counts ~ target + cond + s(time.c, by = cond, k = 7) + offset(log(nf)), 
         data = qdat_pre, 
         family = nb)
plot(m)




m2x <- brm(bf(counts ~  target + target:time + target:cond + target:time:cond + (1|participant) + (1|participant:cdna) + offset(log(nf)), 
              shape ~ target),
           family = negbinomial(),
           warmup = 1000, # number of samples before sampling
           iter = 3000,  # number of mcmc iterations 
           cores = 6, # number of cores used in sampling
           chains = 6, # number of chains
           seed = 123, # seed for reproducible results
           control = list(adapt_delta = 0.8),
           data = qdat_pre)
saveRDS(m2x, "./data/derivedData/temp/temp.RDS")

summary(m2x)
conditional_effects(m2x)

hypothesis(m1x, c("timeS12:condvar = 0"))

pp_check(m2x, type = "ecdf_overlay")


ggeffect(mx)

??ggpredict

qdat_pre$resid <- resid(m)
qdat_pre$fitted <- fitted(m)


qdat_pre %>%
  
    ggplot(aes(fitted, resid )) + geom_point(shape = 21) + 
  geom_text(data = subset(qdat_pre, resid > 2.5|resid < -1.5), aes(label = paste(participant, leg, time)))





em <- emmeans(m, specs = ~ time | cond)

em %>%
  data.frame() %>%
  mutate(time.c = as.numeric(gsub("S", "", time))) %>% 
  
  ggplot(aes(time.c, emmean, color = cond, group = cond)) +
  
  geom_line() 





summary(m)


qdat_pre %>%
  group_by(participant, time.c, time, leg, cond, target) %>%
  
  summarise(nf = mean(nf), 
            counts = mean(counts)) %>%
  
  ggplot(aes(time.c, nf, color = cond,
             group = paste(participant, leg))) + 
  geom_line() +
  geom_point() + 
  facet_grid(target ~ participant)









log(nf$nf)

get_prior(bf(counts ~ cond + s(time.c, by = cond, k = 7) + (1|participant) + offset(log(nf))), 
           data = qdat_pre[qdat_pre$target == "rRNA47S F1R1",])


m1x <- brm(bf(log(norm) ~ cond + s(time.c, by = cond, k = 7) + 
               (1|participant)),

    
  
          warmup = 1000, # number of samples before sampling
          iter = 4000,  # number of mcmc iterations 
          cores = 6, # number of cores used in sampling
          chains = 6, # number of chains
          seed = 123, # seed for reproducible results
          control = list(adapt_delta = 0.95),
          data = qdat_pre)




pp_check(m1)
summary(m1)


plot(m1)


conditional_smooths(m1x)


library(mgcv)

mx <- gam(counts ~ cond + s(time.c, by = cond, k = 6) + 
            
            s(participant, bs = "re"), 
          offset = log(nf),
          family = poisson,

          data = qdat_pre[qdat_pre$target == "rRNA47S F1R1",])


plot(mx)


library(mgcv)

qdat_rrna2 %>%
  filter(rel.counts < 0.001) %>%
  ggplot(aes(time.c, rel.counts)) + geom_point()


m1 <- brm(bf(counts ~ s(time.c, by = target, k = 7) + (1|participant) + (1|technical), 
             sigma ~ target),

          warmup = 1000, # number of samples before sampling
          iter = 6000,  # number of mcmc iterations 
          cores = 4, # number of cores used in sampling
          chains = 4, # number of chains
          seed = 123, # seed for reproducible results
          thin = 5,
          data = qdat_rrna2)

saveRDS(m1, "./data/derivedData/temp/qpcr_smooth.RDS")
m1 <- readRDS("./data/derivedData/temp/qpcr_smooth.RDS")

conditional_effects(m1)




?brm








