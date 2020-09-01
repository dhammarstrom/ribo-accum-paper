##### qpcr analysis using brms #########################

# qpcr analysis in a Bayesian framework.


# Use minimal packages as RStudio seems to chrash when RAM running low
library(brms)
library(tidybayes)
library(readxl)
library(tidyverse)




# brms offers several 






# qdat is the compiled data frame containing sample information and 
# qpcr estimates. Created in qpcr-compile.R. Use this for downstream analyses

qdat <- readRDS("./data/derivedData/qpcr/qpcr_compiled.RDS")


# Remove bad reaction prior to modelling. 
# Bad reactions are no amplification or estimatyed to cq < 5


qdat %>%
  ggplot(aes(cq, color = target)) + geom_histogram() + 
  facet_wrap(~ target)

qdat %>%
  group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 5 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 5 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  
  ggplot(aes(paste(participant, leg, time), cq, color = target, shape = outlier)) + 
  geom_point() + theme_minimal() 

# 5 sd's away from target mean seems to capture the worst reactions.. 





### Lambda normalization ####

nf <- qdat %>%
  filter(target == "Lambda KIT") %>%
  filter(cq > 5) %>% # Removes bad reactions
  mutate(nf.w = (eff ^ -cq) * tissue_weight) %>%
  dplyr::select(participant, leg, time, cdna, nf.w) %>%
  group_by(participant, leg, time, cdna) %>%
  summarise(nf.w = mean(nf.w, na.rm = TRUE)) %>%
  ungroup() %>%
  print()





#### rRNA per tissue weight analysis #### 


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
          # Interaction effects post-training
          "targetrRNA18SF2R2:txpost_int_train -     targetrRNA18SF2R2:txpre_int_train =      targetrRNA18SF2R2:txpost_con_train",
          "targetrRNA45SF5R5:txpost_int_train -     targetrRNA45SF5R5:txpre_int_train =      targetrRNA45SF5R5:txpost_con_train",
          "targetrRNA47SF1R1:txpost_int_train -     targetrRNA47SF1R1:txpre_int_train =      targetrRNA47SF1R1:txpost_con_train",
          "targetrRNA5SF3R3:txpost_int_train -      targetrRNA5SF3R3:txpre_int_train =       targetrRNA5SF3R3:txpost_con_train",
          "targetrRNA28SF2R2:txpost_int_train -     targetrRNA28SF2R2:txpre_int_train =      targetrRNA28SF2R2:txpost_con_train",
          "targetrRNA45SITSF12R12:txpost_int_train -targetrRNA45SITSF12R12:txpre_int_train = targetrRNA45SITSF12R12:txpost_con_train", 
          "targetrRNA5.8SF2R2:txpost_int_train -    targetrRNA5.8SF2R2:txpre_int_train =     targetrRNA5.8SF2R2:txpost_con_train", 
          # Interaction effects detraining (post1w)
          "targetrRNA18SF2R2:txpost_int_detrain -     targetrRNA18SF2R2:txpre_int_train =      targetrRNA18SF2R2:txpost_con_train",
          "targetrRNA45SF5R5:txpost_int_detrain -     targetrRNA45SF5R5:txpre_int_train =      targetrRNA45SF5R5:txpost_con_train",
          "targetrRNA47SF1R1:txpost_int_detrain -     targetrRNA47SF1R1:txpre_int_train =      targetrRNA47SF1R1:txpost_con_train",
          "targetrRNA5SF3R3:txpost_int_detrain -      targetrRNA5SF3R3:txpre_int_train =       targetrRNA5SF3R3:txpost_con_train",
          "targetrRNA28SF2R2:txpost_int_detrain -     targetrRNA28SF2R2:txpre_int_train =      targetrRNA28SF2R2:txpost_con_train",
          "targetrRNA45SITSF12R12:txpost_int_detrain -targetrRNA45SITSF12R12:txpre_int_train = targetrRNA45SITSF12R12:txpost_con_train",
          "targetrRNA5.8SF2R2:txpost_int_detrain -    targetrRNA5.8SF2R2:txpre_int_train =     targetrRNA5.8SF2R2:txpost_con_train", 
          
          # Whitin group fold change
          "targetrRNA18SF2R2:txpost_int_detrain =     targetrRNA18SF2R2:txpre_int_train",     
          "targetrRNA45SF5R5:txpost_int_detrain =     targetrRNA45SF5R5:txpre_int_train",     
          "targetrRNA47SF1R1:txpost_int_detrain =     targetrRNA47SF1R1:txpre_int_train",     
          "targetrRNA5SF3R3:txpost_int_detrain =      targetrRNA5SF3R3:txpre_int_train",      
          "targetrRNA28SF2R2:txpost_int_detrain =     targetrRNA28SF2R2:txpre_int_train",    
          "targetrRNA45SITSF12R12:txpost_int_detrain = targetrRNA45SITSF12R12:txpre_int_train",
          "targetrRNA5.8SF2R2:txpost_int_detrain =    targetrRNA5.8SF2R2:txpre_int_train", 
          
          "targetrRNA18SF2R2:txpost_int_train =     targetrRNA18SF2R2:txpre_int_train",     
          "targetrRNA45SF5R5:txpost_int_train =     targetrRNA45SF5R5:txpre_int_train",     
          "targetrRNA47SF1R1:txpost_int_train =     targetrRNA47SF1R1:txpre_int_train",     
          "targetrRNA5SF3R3:txpost_int_train =      targetrRNA5SF3R3:txpre_int_train",      
          "targetrRNA28SF2R2:txpost_int_train =     targetrRNA28SF2R2:txpre_int_train",    
          "targetrRNA45SITSF12R12:txpost_int_train = targetrRNA45SITSF12R12:txpre_int_train",
          "targetrRNA5.8SF2R2:txpost_int_train =    targetrRNA5.8SF2R2:txpre_int_train", 
          
          "targetrRNA18SF2R2:txS1_int_train =     targetrRNA18SF2R2:txpre_int_train",     
          "targetrRNA45SF5R5:txS1_int_train =     targetrRNA45SF5R5:txpre_int_train",     
          "targetrRNA47SF1R1:txS1_int_train =     targetrRNA47SF1R1:txpre_int_train",     
          "targetrRNA5SF3R3:txS1_int_train =      targetrRNA5SF3R3:txpre_int_train",      
          "targetrRNA28SF2R2:txS1_int_train =     targetrRNA28SF2R2:txpre_int_train",    
          "targetrRNA45SITSF12R12:txS1_int_train = targetrRNA45SITSF12R12:txpre_int_train",
          "targetrRNA5.8SF2R2:txS1_int_train =    targetrRNA5.8SF2R2:txpre_int_train", 
          
          "targetrRNA18SF2R2:txS1_con_train = 0",
          "targetrRNA45SF5R5:txS1_con_train = 0",
          "targetrRNA47SF1R1:txS1_con_train = 0",
          "targetrRNA5SF3R3:txS1_con_train = 0",
          "targetrRNA28SF2R2:txS1_con_train = 0",
          "targetrRNA45SITSF12R12:txS1_con_train = 0", 
          "targetrRNA5.8SF2R2:txS1_con_train = 0", 
          
          "targetrRNA18SF2R2:txpost_con_train = 0",
          "targetrRNA45SF5R5:txpost_con_train = 0",
          "targetrRNA47SF1R1:txpost_con_train = 0",
          "targetrRNA5SF3R3:txpost_con_train = 0",
          "targetrRNA28SF2R2:txpost_con_train = 0",
           "targetrRNA45SITSF12R12:txpost_con_train = 0",
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


























