
# mcmcglmm modeling of qpcr
#
#
#


# Notes:

# Matz suggested using a poisson mixed-effects model to model qpcr-data 
# without the use of reference genes as random effects on technical duplictaes
# could be used to control for technical variation (Matz el al. 2013).

# An implementation of these models is done in the mcmc.qpcr package. Althogh 
# this is a nice implementation we would like to have more control in the 
# fitting process and downstream comparisons. To this end we use the 
# mcmcglmm package directly. 

# Two models are used, in model 1 estimates are retrieved from an "unadjusted model"
# modeling transcript counts per total RNA (as this is the biological input). Model 2 
# uses the external reference gene as a normalization factor modeling counts as a rate
# per tissue weight. The normalization factor is entered as a fixed effect with strong priors
# with the interpretation being counts per normalization factor. 
# See discussion on offset in MCMCglmm:
# https://hannahdugdale.wordpress.com/2016/03/16/why-you-shouldnt-use-the-offset-function-in-mcmcglmm/
# https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q3/016817.html


# Load packages 
library(MCMCglmm)
library(tidybayes)
library(readxl)
library(tidyverse)
library(parallel)
library(coda)



# qdat is the compiled data frame containing sample information and 
# qpcr estimates. Created in qpcr-compile.R. Use this for downstream analyses

qdat <- readRDS("./data/derivedData/qpcr/qpcr_compiled.RDS")



# Remove bad reactions prior to modeling. 

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

# 2 sd's away from target mean seems to capture the worst reactions.. 


### Lambda normalization factor ####


nf <- qdat %>% group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 2 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 2 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  filter(target == "Lambda KIT") %>%
  filter(outlier == "in") %>% # Removes bad reactions

  # The raw data is trasnformed to counts and 
  # thereafter multiplied with the input muscle 
  # weight to create the normalization factor
  # this is proportional to the amount of muscle 
  # used in each cDNA synthesis.
  
  mutate(nf.w = (eff ^ (37-cq)) * tissue_weight) %>%
  dplyr::select(participant, leg, time, cdna, nf.w) %>%
  group_by(participant, leg, time, cdna) %>%
  # Sets the norm factor on the same scale as the other targets
  summarise(nf.w = as.integer(round(mean(nf.w, na.rm = TRUE),0))) %>%
  ungroup() %>%
  # Scale the factor
  print()


  



#### rRNA per tissue weight analysis #### 

qdat.rrna  <- qdat %>%
  group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 2 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 2 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  filter(outlier == "in") %>%
  
  group_by(participant, leg, time, sex, cond, cdna, target) %>%
  summarise(cq = mean(cq, na.rm = TRUE), 
            eff = mean(eff, na.rm = TRUE), 
            tissue_weight = mean(tissue_weight, na.rm = TRUE)) %>%
  # To calibrate the models we use all avaliable transcripts. 
  filter(target %in% c("rRNA18S F2R2",      
                       "rRNA28S F2R2",     
                       "rRNA5.8S F2R2",       
                       "rRNA45S F5R5",      
                       "rRNA45SITS F12R12", 
                       "MyHC1 F1R1", 
                       "MyHC2A F5R5",
                       "MyHC2X F5R5",  
                       "rRNA47S F1R1",      
                       "rRNA5S F3R3", 
                       "UBTF F4R4", 
                       "UBTF F6R6", 
                       "rpS6 F2R2")) %>%

  
  mutate(Ra = -cq * log(eff), # Relative abundance
         counts = as.integer(round(eff ^ (37-cq), 0)), # counts, ref Matz et al.
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
  print()




### Summarise counts over all transcripts
qdat.rrna %>%
  group_by(target) %>%
  summarise(m = mean(counts),
            min = min(counts),
            max = max(counts)) %>%
  print()

# Looks good. Big differences as expected...

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
  data.frame() %>%
  print()


# Notes on modeling:
# A random intercept for the participant and for the technical replicate. 
#
# The mcmc.qpcr approach adds slopes by random factor for the gene effect. 

# The priors for the variance components are set with the limit (V) set to 1 and the 
# belief parameter (nu) set to 0.002. The course notes (p. 13) says this is a common
# use for variance components.
# A  This prior is specified in the R part. 
# We set this to nu = 4.002 to make as similar to as possible mcmc.qpcr 

# Random effects allow for random intercepts by subject, technical sample and biological sample
# and random gene slopes by subject.

# Library parallel needed for parallel fitting of models

prior1 <- list(R = list(V = diag(13), # n genes/targets 
                        nu = 4.002),
               # the G part are the variance components for the 
               # random effects
               G = list(G1 = list(nu = 0.002, V = 1), 
                        G2 = list(nu = 0.002, V = 1),
                        G3 = list(nu = 0.002, V = diag(13))))


# use only 80% of the machine's total 
# logical processors. 

setCores <- 4 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

# load the MCMCglmm package within the cluster
cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) 
# import each object that's necessary to run the function
clusterExport(cl, "prior1")
clusterExport(cl, "qdat.ctrl")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
model1 <- parLapply(cl=cl, 1:4, function(i) {
  
  MCMCglmm(fixed = counts ~ 0 + target + target:tx,
           random = ~participant + technical + idh(target):participant, 
           prior = prior1, 
           family = "poisson", 
           data = qdat.ctrl, 
           rcov = ~idh(target):units,
           thin = 10,
           burnin = 5000,
           nitt = 50000)
  
  }
)

# once it's finished, use stopCluster() to stop running the parallel cluster
stopCluster(cl)





## In model 2 we will use the normalization factor as an "offset"
# specified with strong priors as suggested in 
# https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q3/016817.html

# New priors are specified

prior2 <- list(R = list(V = diag(13), # n genes/targets 
                        nu = 4.002),
               B = list(V = diag(92) * 1e5, 
                        mu = c(1, rep(0, 91))),
               # the G part are the variance components for the 
               # random effects
               G = list(G1 = list(nu = 0.002, V = 1), 
                        G2 = list(nu = 0.002, V = 1),
                        G3 = list(nu = 0.002, V = diag(13))))


prior2$B$V[1,1] <- 1e-5 # replace the offset variance



setCores <- 4 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

# load the MCMCglmm package within the cluster
cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) 
# import each object that's necessary to run the function
clusterExport(cl, "prior2")
clusterExport(cl, "qdat.ctrl")


# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
model2 <- parLapply(cl=cl, 1:4, function(i) {
  
  MCMCglmm(fixed = counts ~ 0 + log(nf.w) + target + target:tx,
           random = ~participant + technical + idh(target):participant, 
           prior = prior2, 
           family = "poisson", 
           data = qdat.ctrl, 
           rcov = ~idh(target):units,
           thin = 10,
           burnin = 5000,
           nitt = 50000)
  
}
)

summary(model2[[1]])


# once it's finished, use stopCluster() to stop running the parallel cluster
stopCluster(cl)



# Combine all chains for diagnostics

m1 <- lapply(model1, function(m) m$Sol)

m1 <- do.call(mcmc.list, m1)


m2 <- lapply(model2, function(m) m$Sol)

m2 <- do.call(mcmc.list, m2)


# The gelman plot from 
gelman.diag(m1)
gelman.diag(m2)

plot(m1)






#### Extract contrasts #######################

# Defines a function to get model/target specific estimates 

contr_extr <- function(model, target) {
  
  # Define the model
  mc1 <- model
  
  # Comparisons are based on these factors
  comps <- c("txS1_con_train", 
             "txpost_con_train",
             
             "txpre_int_train", 
             "txS1_int_train", 
             "txpost_int_train", 
             "txpost_int_detrain")
  
  # Make new vector containg target info
  comps2 <- paste0("target", target,  ":",  comps)
  
  inter_S1          <- as.numeric( (mc1$Sol[, comps2[4]] - mc1$Sol[, comps2[3]]) -  mc1$Sol[, comps2[1]])
  inter_post        <- as.numeric( (mc1$Sol[, comps2[5]] - mc1$Sol[, comps2[3]]) -  mc1$Sol[, comps2[2]])
  inter_post1w      <- as.numeric( (mc1$Sol[, comps2[6]] - mc1$Sol[, comps2[3]]) -  mc1$Sol[, comps2[2]])
  
  within_int_S1     <-  as.numeric( (mc1$Sol[, comps2[4]] - mc1$Sol[, comps2[3]]))
  within_int_post   <-  as.numeric( (mc1$Sol[, comps2[5]] - mc1$Sol[, comps2[3]]))
  within_int_post1w <-  as.numeric( (mc1$Sol[, comps2[6]] - mc1$Sol[, comps2[3]]))
  
  within_con_S1     <-  as.numeric( mc1$Sol[, comps2[1]] )
  within_con_post   <-  as.numeric( mc1$Sol[, comps2[2]] )
  
  # Combine all values in a df
  contr <- data.frame(cbind(inter_S1,         
                            inter_post,       
                            inter_post1w,     
                            within_int_S1,    
                            within_int_post,  
                            within_int_post1w,
                            within_con_S1,    
                            within_con_post))
  
 results <-  contr %>%
    pivot_longer(names_to = "coef", values_to = "val", cols = 1:8) %>%
    group_by(coef) %>%
    summarise(est = median(val), 
              lwr = quantile(val, 0.025), 
              upr = quantile(val, 0.975))  %>%
   mutate(target = target) %>%
   dplyr::select(target, coef, est, lwr, upr)
  
  return(results)
}


targets <- unique(qdat.ctrl$target)

results <- list()

for(i in 1:length(targets)) {
  
 m1.res  <- contr_extr(model1[[1]], targets[i]) %>%
   mutate(model = "total_rna")
 m2.res  <- contr_extr(model2[[1]], targets[i]) %>%
   mutate(model = "tissue")
  
 results[[i]] <- rbind(m1.res, m2.res)
 
}

# Combine all results
ctrl_vs_int <- bind_rows(results) %>%
  mutate(contrast = if_else(coef == "inter_S1", 
                            "inter:S1", 
                            if_else(coef == "inter_post", 
                                    "inter:post", 
                                    if_else(coef == "inter_post1w", 
                                            "inter:post1w", coef))), 
         contrast = gsub("within_", "", contrast)) %>%
  dplyr::select(contrast, target, estimate = est, lower.CL = lwr, upper.CL = upr, model) %>%
  print()
                                       


### Save results 


saveRDS(ctrl_vs_int, "./data/derivedData/qpcr-analysis-bayes2/qpcr_res_int_con.RDS")


############## Within-training-group models ##############################


qdat_int <- qdat.rrna %>%
  filter(cond != "ctrl_leg") %>%
  
  mutate(time = factor(time, 
                       levels = c("S0", 
                                  "S1", 
                                  "S4", 
                                  "S5", 
                                  "S8", 
                                  "S9", 
                                  "S12", 
                                  "post1w"))) %>%
  print()




prior3 <- list(R = list(V = diag(13), # n genes/targets 
                        nu = 4.002),
               # the G part are the variance components for the 
               # random effects
               G = list(G1 = list(nu = 0.002, V = 1), 
                        G2 = list(nu = 0.002, V = 1),
                        G3 = list(nu = 0.002, V = diag(13))))


setCores <- 4 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

# load the MCMCglmm package within the cluster
cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) 
# import each object that's necessary to run the function
clusterExport(cl, "prior3")
clusterExport(cl, "qdat_int")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
model3 <- parLapply(cl=cl, 1:4, function(i) {
  
  MCMCglmm(fixed = counts ~ 0 + target + target:time + target:time:cond,
           random = ~participant + technical + idh(target):participant, 
           prior = prior3, 
           family = "poisson", 
           data = qdat_int, 
           rcov = ~idh(target):units,
           thin = 10,
           burnin = 10000,
           nitt = 65000)
  
}
)

# once it's finished, use stopCluster() to stop running the parallel cluster
stopCluster(cl)




# New priors are specified

prior4 <- list(R = list(V = diag(13), # n genes/targets 
                        nu = 4.002),
               B = list(V = diag(209) * 1e5, 
                        mu = c(1, rep(0, 208))),
               # the G part are the variance components for the 
               # random effects
               G = list(G1 = list(nu = 0.002, V = 1), 
                        G2 = list(nu = 0.002, V = 1),
                        G3 = list(nu = 0.002, V = diag(13))))


prior4$B$V[1,1] <- 1e-5 # replace the offset variance


setCores <- 4 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

# load the MCMCglmm package within the cluster
cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) 
# import each object that's necessary to run the function
clusterExport(cl, "prior4")
clusterExport(cl, "qdat_int")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
model4 <- parLapply(cl=cl, 1:4, function(i) {
  
  MCMCglmm(fixed = counts ~ 0 + log(nf.w) + target + target:time + target:time:cond,
           random = ~participant + technical + idh(target):participant, 
           prior = prior4, 
           family = "poisson", 
           data = qdat_int, 
           rcov = ~idh(target):units,
           thin = 10,
           burnin = 10000,
           nitt = 65000)
  
}
)

# once it's finished, use stopCluster() to stop running the parallel cluster
stopCluster(cl)


##### Mature vs pre-rna ##########################
# An attempt to model pre-rrna and mature rna together 
# Warning, experimental!


qdat_int2 <- qdat_int %>%
  filter(target != "rRNA5.8S F2R2") %>%
  mutate(Target = if_else(target %in% c("rRNA45S F5R5", 
                                        "rRNA47S F1R1"), "pre", 
                          if_else(target %in% c("rRNA18S F2R2", 
                                                "rRNA5.8S F2R2", 
                                                "rRNA28S F2R2"), "mature", target))) %>%
  print()

qdat_int2 %>%
  filter(target == "rRNA18S F2R2") %>%
  dplyr::select(target, Target) %>%
  print()


prior5 <- list(R = list(V = diag(12), # n genes/targets 
                        nu = 4.002),
               B = list(V = diag(161) * 1e5, 
                        mu = c(1, rep(0, 160))),
               # the G part are the variance components for the 
               # random effects
               G = list(G1 = list(nu = 0.002, V = 1), 
                        G2 = list(nu = 0.002, V = 1),
                        G3 = list(nu = 0.002, V = diag(12))))



prior5$B$V[1,1] <- 1e-5 # replace the offset variance


setCores <- 4 # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores", setCores))

# load the MCMCglmm package within the cluster
cl.pkg <- clusterEvalQ(cl, library(MCMCglmm)) 
# import each object that's necessary to run the function
clusterExport(cl, "prior5")
clusterExport(cl, "qdat_int2")

# use parLapply() to execute 10 runs of MCMCglmm(), each with nitt=100000
model5 <- parLapply(cl=cl, 1:4, function(i) {
  
  MCMCglmm(fixed = counts ~ 0 + log(nf.w) + Target + Target:time + Target:time:cond,
           random = ~participant + technical + idh(target):participant, 
           prior = prior5, 
           family = "poisson", 
           data = qdat_int2, 
           rcov = ~idh(target):units,
           thin = 10,
           burnin = 10000,
           nitt = 65000)
  
}
)

# once it's finished, use stopCluster() to stop running the parallel cluster
stopCluster(cl)

















# Combine all chains for diagnostics

m3 <- lapply(model3, function(m) m$Sol)

m3 <- do.call(mcmc.list, m3)
m4 <- lapply(model4, function(m) m$Sol)

m4 <- do.call(mcmc.list, m4)

m5 <- lapply(model5, function(m) m$Sol)

m5 <- do.call(mcmc.list, m5)

# The gelman plot from 
gelman.diag(m3)
gelman.diag(m4)
gelman.diag(m5)



# Combine data per target from regression coefficients #


reg_cof <- rbind(cbind(summary(m4)$statistics, summary(m4)$quantiles[,c(1,5)]) %>%
   data.frame() %>%
   mutate(coef = rownames(.), 
          model = "tissue"), 
   cbind(summary(m3)$statistics, summary(m3)$quantiles[,c(1,5)]) %>%
     data.frame() %>%
     mutate(coef = rownames(.), 
            model = "tot_rna"), 
   cbind(summary(m5)$statistics, summary(m5)$quantiles[,c(1,5)]) %>%
     data.frame() %>%
     mutate(coef = rownames(.), 
            model = "combine_rrna")) %>%
   print()



interaction_rrna <- reg_cof %>%
  separate(coef, into = c("target", "time", "cond"), sep = ":") %>%
  mutate(target = gsub("target", "", target), 
         target = gsub("Target", "", target)) %>%
  dplyr::select(model, target, time, cond, Mean, lwr =  X2.5., upr = X97.5.) %>%
  filter(target %in% c("rRNA47S F1R1",
                       "rRNA45S F5R5",
                       "rRNA45SITS F12R12",
                       "rRNA18S F2R2", 
                       "rRNA5.8S F2R2", 
                       "rRNA28S F2R2", 
                       "rRNA5S F3R3",
                       "pre", "mature")) %>%

  mutate(time = gsub("time", "", time), 
         time = factor(time, levels = c("S0", 
                                          "S1", 
                                          "S4", 
                                          "S5", 
                                          "S8", 
                                          "S9", 
                                          "S12", 
                                          "post1w")), 
         robust = if_else(lwr > 0 | upr < 0, "robust", "notrobust")) %>%
  
  
  filter(cond == "condvar", 
         model %in% c("tissue", "combine_rrna"))%>%
  print()
  



estimated_means_rrna <- rbind(emmeans(model3[[2]], specs = ~ time * cond|target, data = qdat_int) %>%
  data.frame() %>%
    mutate(model = "total_rna"), 
  emmeans(model4[[2]], specs = ~ time * cond|target, data = qdat_int) %>%
    data.frame() %>%
    mutate(model = "tissue"),
  emmeans(model5[[2]], specs = ~ time * cond|Target, data = qdat_int2) %>%
  data.frame() %>%
  mutate(model = "combine_rrna") %>% 
  dplyr::select(time, cond, target = Target, emmean:model) %>%
    filter(target %in% c("pre", "mature"))) %>%
  mutate(time = factor(time, levels = c("S0", 
                                        "S1", 
                                        "S4",
                                        "S5", 
                                        "S8", 
                                        "S9", 
                                        "S12", 
                                        "post1w"))) %>%
  filter(model %in% c("tissue","combine_rrna"),  
         target %in% c("rRNA47S F1R1",
                       "rRNA45S F5R5",
                       "rRNA45SITS F12R12",
                       "rRNA18S F2R2", 
                       "rRNA5.8S F2R2", 
                       "rRNA28S F2R2",
                       "rRNA5S F3R3",
                       "pre", 
                       "mature")) %>%
  print()
  
  
  

## Save data 

rrna_timecourse <- list(estimated_means = estimated_means_rrna, 
                        estimated_diff = interaction_rrna)


saveRDS(rrna_timecourse, "./data/derivedData/qpcr-analysis-bayes2/time_course_qpcr.RDS")



  
  
  