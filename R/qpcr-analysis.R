######################### qPCR analysis #############################


source("./R/libs.R")

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
# In comparing control vs. intervention there is an issue with rank deficiency
# nlme:lme does not handle rank deficiency but do handle variance functions 
# which would allow for multi-target modeling further accounting for 
# technical variation through random effects with each technical replicate. 

# lme4::lmer does not handle variance componenets per strata (gene/target)
# mcmc.qpcr is a possibility which uses mcmcglmm, however post-hoc comparisons
# are not straight-forward. 

# A possible solution is 

qdat.ctrl <- qdat.rrna %>%
  filter(time %in% c("S0", "S1c", "S1", "S12", "postctrl", "post1w")) %>%
  mutate(group = if_else(cond == "ctrl_leg", "con", "int"), 
         detrain = if_else(time == "post1w", "detrain", "train"), 
         time = if_else(time == "S0", "pre", 
                        if_else(time == "S1c", "S1", 
                                if_else(time %in% c("postctrl", "post1w", "S12"), "post", time))), 
         time = factor(time, levels = c("pre", "S1", "post")))  %>%
  
  
  
  print()


# brms solution ########################

# THIS IS WHERE IM AT ! #################


library(brms)


m1 <- brm(counts ~ 0 + target + target:time + target:group + 
            target:time:group + target:time:group:detrain +
            (target||participant) + (1|technical) + (1|biological), 
          data = qdat.ctrl, 
          family = poisson, 
          iter = 50000, 
          cores = 16)


h <- hypothesis(m1, "targetrRNA18SF2R2:timepost:groupint - targetrRNA18SF2R2:groupint= targetrRNA18SF2R2:timepost - targetrRNA18SF2R2")


h$hypothesis
summary(m1)


### mcmcglmm solution #############

library(MCMCglmm); library(parallel)

### The first model is specified with random effects as in the basic mixed linear-model implementation. 
# A random intercept for the participant and for the technical replicate. 
#
# The mcmc.qpcr approach adds slopes by random factor for the gene effect. First we ignore this to get a model formulation 
# more like the one in the lme-approach.


# The priors for the variance components are set with the limit (V) set to 1 and the 
# belief parameter (nu) set to 0.002. The course notes (p. 13) says this is a common
# use for variance components.
# A  This prior is specified in the R part. 
# We set this to nu = 4.002 to make as similar to as possible mcmc.qpcr 

# Random effects allow for random intercepts by subject, technical sample and biological sample
# and random gene slopes by subject.

# Library parallel needed for parallel fitting of models
library(parallel)


prior1 <- list(R = list(V = diag(7), nu = 4.002),
               G = list(G1 = list(nu = 0.002, V = 1), 
                        G2 = list(nu = 0.002, V = 1),
                        G3 = list(nu = 0.002, V = 1),
                        G4 = list(nu = 0.002, V = diag(7))))


setCores<-round(detectCores()*1) # use detectCores() by itself if you want all CPUs

# make the cluster
cl <- makeCluster(getOption("cl.cores",setCores))


# load the MCMCglmm package within the cluster
cl.pkg <- clusterEvalQ(cl,library(MCMCglmm)) 

clusterExport(cl,"prior1")
clusterExport(cl,"qdat.ctrl")


model2_10runs <- parLapply(cl=cl,1:10, function(i) {
  MCMCglmm(fixed = counts ~ 0 + target + target:time + target:group + target:time:group + target:time:group:detrain ,
           random = ~participant + technical + biological + idh(target):participant, 
           prior = prior1, 
           family = "poisson", 
           data = qdat.ctrl, 
           rcov = ~idh(target):units,
           thin = 20,
           burnin = 5000,
           nitt = 55000)}
)

# once it's finished, use stopCluster() to stop running the parallel cluster
stopCluster(cl)





mc1.models.diagn <- do.call(mcmc.list, lapply(model2_10runs, function(m) m$VCV))

par(mfrow=c(4,2), mar=c(2,2,1,2))
gelman.plot(mc1.models.diagn, auto.layout=F)


summary(model2_10runs[[1]])

## Examining the gelman plots and scale reduction factors suggests that the model is ok. 
# The model converges. 

mc1 <- model2_10runs[[1]]


mc1$

model2_10runs[[1]]


con.post <- matrix(c(1, 1, 0, 0, 0),1) - matrix(c(1, 0, 0, 0, 0),1) 
int.post <- (matrix(c(1, 1, 1, 1, 0),1) - matrix(c(1, 0, 1, 0, 0),1))
int.post.detrain <- (matrix(c(1, 1, 1, 1, 1),1) - matrix(c(1, 0, 1, 0, 0),1)) 

k <-  rbind(con.post, 
            int.post, 
            int.post.detrain,
            int.post - con.post, 
            int.post.detrain - con.post) 
# Set rownames
rownames(k) <- c("post_con", "post_int", "post1w_int", "inter:post_int", "inter:post1w_int") 

## These confidence intervals are not adjusted. 
ci.isok <- confint(glht(m3, linfct = k), calpha = univariate_calpha())

















## MCMC diagnostics 
# First we plot traces and density-plots for all parameters


  mcmc_trace(mc1$VCV)
  mcmc_hist(mc1$VCV)
  
  mcmc_trace(mc1$Sol)
  mcmc_hist(mc1$Sol)
  
  # These looks good, traces stationary and seems to converge. How about autocorrelations?
  
  mcmc_acf(mc1$VCV)
  mcmc_acf(mc1$Sol)




# No autocorrelation at lag >~10... A good sign!

# How about the variance components?
vcv_sum(mc1$VCV)

# The variance components are well above zero (95% credible interval), this is indication 
# to keep them. 

# In the second model we replicate the "full" model 5 in the nlme approach. We include random slopes on the 
# subject level and biological sample level.

prior2 <- list(R = list(V = diag(4), nu = 4.002),
               G = list(G1 = list(nu = 0.002, V = 1), 
                        G2 = list(nu = 0.002, V = 1),
                        G3 = list(V = diag(4), nu = 0.002),
                        G4 = list(V = diag(4), nu = 0.002),
                        G5 = list(V = 1, nu = 0.002)))


mc2.models <- mclapply(1:4, function(i) {
  MCMCglmm(fixed = count ~ 0 + gene + gene:timepoint + gene:sets + gene:timepoint:sets,
           random = ~subject + technical + idh(gene):subject + idh(gene):biological + biological, 
           prior = prior2, 
           family = "poisson", 
           data = qs, 
           rcov = ~idh(gene):units,
           thin = 20,
           burnin = 5000,
           nitt = 55000)
}, mc.cores=1)


mc2.models.diagn <- do.call(mcmc.list, lapply(mc2.models, function(m) m$Sol))

par(mfrow=c(4,2), mar=c(2,2,1,2))
# gelman.plot(mc2.models.diagn, auto.layout=F)
if(USE_mcmcPLOTS == TRUE){
  gelman.diag(mc2.models.diagn)
}
## Examining the gelman plots and scale reduction factors suggests that the model is ok. 
# The model converges. 
mc2 <-mc2.models[[1]]

# Assessing the model fit by DIC
DIC(list(mc1 = mc1, mc2 = mc2))

# Suggests that model mc 2 is better

## Looking at the diagnostics for the second model
# Traces

if(USE_mcmcPLOTS == TRUE){
  mcmc_trace(mc2$VCV) + geom_smooth(color = "black")
  mcmc_hist(mc2$VCV)
  
  ## The variance components are small, but all traces seems ok.
  
  # We can also look at effective samples
  effectiveSize(mc2$VCV);effectiveSize(mc2$Sol)
  
  # Effective size looks alright even though some VCV-components seems to be sampled less effective. Increse
  # iterations in the next modeling (nitt = 55000, burning = 5000, thin = 20)
  
  # Looking at autocorrelations
  mcmc_acf(mc2$Sol)
  mcmc_acf(mc2$VCV)}
# Also autocorr looks ok.











m1 <- lme(log(Ra.nf) ~ 0+target + target:time + target:group + target:time:group, 
          random = list(participant = ~ 1, 
                        technical = ~ 1), 
          weights = varIdent(form = ~ 1|target), 
          data = qdat.ctrl[qdat.ctrl$detrain == "train",],
          control = list(msMaxIter = 120, 
                         opt = "nloptwrap", 
                         msVerbose = TRUE), 
          method = "REML", na.action = na.exclude) 
          
# In model 1, residual variances are modeled per target. To test if more ellaborate 
# modelling is needed.


# varfun3 <- varPower(form = ~ cq|target)
# varfun4 <- varExp(form = ~ cq|target)
# 
# 
# m3 <- update(m1, weights = varfun3)
# m4 <- update(m1, weights = varfun4)

# m5 <- update(m1, random = list(participant = pdDiag(~ 1 + target), technical = ~ 1))

# by testing intervals(m3) shows non pos def vcov. Too complex, stick with m1.

# Residual plots looks good enough
plot(m1, resid(., type = "n") ~ fitted(.)|target)
qqnorm(m1, ~ resid(., type = "n")|target, abline = c(0,1))


# The same model is used to compare control pre-post to intervention pre-detrain

m1.detrain <- lme(log(Ra.nf) ~ 0+ target + target:time + target:group + target:time:group, 
                  random = list(participant = ~ 1, 
                                technical = ~ 1), 
                  weights = varIdent(form = ~ 1|target), 
                  data = qdat.ctrl[!(qdat.ctrl$detrain == "train" & 
                                       qdat.ctrl$time == "post" &
                                       qdat.ctrl$group == "int"),],
                  control = list(msMaxIter = 120, 
                                 opt = "nloptwrap", 
                                 msVerbose = TRUE), 
                  method = "REML", na.action = na.exclude) 


### Get coefficients from the model using custom contrasts 

data.frame(pairs(emmeans(m1, specs = ~ time*group|target), reverse = TRUE)) %>%
 filter(contrast %in% c()) 

intervals(m1)$fixed


colnames(model.matrix(m1, data = qdat.ctrl[qdat.ctrl$detrain == "train",]))



con.post <- matrix(c(1, 1, 0, 0, 0),1) - matrix(c(1, 0, 0, 0, 0),1) 
int.post <- (matrix(c(1, 1, 1, 1, 0),1) - matrix(c(1, 0, 1, 0, 0),1))
int.post.detrain <- (matrix(c(1, 1, 1, 1, 1),1) - matrix(c(1, 0, 1, 0, 0),1)) 

k <-  rbind(con.post, 
            int.post, 
            int.post.detrain,
            int.post - con.post, 
            int.post.detrain - con.post) 
# Set rownames
rownames(k) <- c("post_con", "post_int", "post1w_int", "inter:post_int", "inter:post1w_int") 

## These confidence intervals are not adjusted. 
ci.isok <- confint(glht(m3, linfct = k), calpha = univariate_calpha())








temp <- qdat.rrna %>%
  group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 5 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 5 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%

  filter(outlier == "in", 
         cond != "ctrl_leg", 
         time != "post1w") %>%
  mutate(target = factor(target), 
         participant = factor(participant), 
         cond = factor(cond), 
         technical = factor(technical))
 

temp %>%
  mutate(ra.tissue = Ra - log(nf.w)) %>%
  ggplot(aes(time.c, ra.tissue, group = paste(target, cond), color = cond)) + 
  geom_point() + facet_wrap(~ target, scales = "free") + 
  geom_smooth()



temp <- qdat.rrna %>%
  group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 5 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 5 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  filter(outlier == "in", 
         time != "post1w") %>%
  mutate(target = factor(target), 
         participant = factor(participant), 
         cond = factor(cond), 
         technical = factor(technical)) %>%

  mutate(time = if_else(time == "S1c", "S1",
                        if_else(time == "postctrl", "S12", time))) %>%
  filter(time %in% c("S0", "S1", "S12")) %>%
  ungroup() %>%
  mutate(nf.w.centered = nf.w - mean(nf.w, na.rm = TRUE), 
         ra.tissue = Ra - log(nf.w)) %>%
  print()        




fxd <- ra.tissue ~ 0 + target + target:time + target:cond + target:time:cond
rand <- list(participant = ~ 1, technical = ~ 1)  


m1 <- lme(fxd, random = rand, data = temp, 
          control = list(msMaxIter = 120, 
                         opt = "nloptwrap", 
                         msVerbose = TRUE), 
          method = "REML", na.action = na.exclude)

varfun2 <- varIdent(form = ~ 1|target)
varfun3 <- varPower(form = ~ cq|target)
varfun4 <- varExp(form = ~ cq|target)

m2 <- update(m1, weights = varfun2)
m3 <- update(m1, weights = varfun3)
m4 <- update(m1, weights = varfun4)


anova(m1,m2, m3, m4)

m5 <- update(m4, random = list(participant = pdDiag(~ 1 + target), technical = ~ 1))

anova(m4, m5)
intervals(m3)

plot(m3, resid(., type = "n") ~ fitted(.)|target)
qqnorm(m3, ~ resid(., type = "n")|target, abline = c(0,1))

intervals(m3)


rna.interaction <- data.frame(coef(summary(m3))) %>%
  mutate(coef = rownames(.)) %>%

  separate(coef, into = c("target", "time", "cond"), sep = ":") %>%
  print()
  


emmeans(m1, specs = ~ "target|time+cond") %>%
  data.frame() %>%
  mutate(time = factor(time, levels = c("S0", "S1", "S12"))) %>%
  ggplot(aes(time, emmean, fill = cond)) + 
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                  position = position_dodge(width = 0.1), 
                  width = 0.05) +
    geom_point(shape = 21, position = position_dodge(width = 0.1), size = 2.5) +
    facet_wrap(~ target, scales = "free")






