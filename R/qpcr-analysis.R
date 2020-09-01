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
         tx = paste(time, group, detrain, sep = "_"))  %>%
  
  print()






m1 <- lme(log(Ra.nf) ~ 0 + target + target:tx, 
          random = list(participant = ~ 1, 
                        technical = ~ 1), 
          weights = varIdent(form = ~ 1|target), 
          control = list(msMaxIter = 120, 
                         opt = "nloptwrap", 
                         msVerbose = TRUE), 
          method = "REML", 
          na.action = na.exclude,
          data = qdat.ctrl)

# In model 1, residual variances are modeled per target. To test if more ellaborate 
# modelling is needed.


varfun2 <- varPower(form = ~ cq|target)
varfun3 <- varExp(form = ~ cq|target)
# 
# 
 m2 <- update(m1, weights = varfun2)
 m3 <- update(m1, weights = varfun3)

 # Check all models for pos def vcov using intervals. 
 intervals(m3)
 
 # m2 does not give a pos def vcov. Comapre m1 and m3
 
 anova(m1, m3)
 
 # m3 gives a better fit, try more elaborate random effects
 
 m4 <- update(m3, random = list(participant = pdDiag(~ 1 + target), technical = ~ 1))
 m5 <- update(m3, random = list(participant = ~ 1, biological = ~ 1, technical = ~ 1))
 intervals(m5)

 # m3 keeps beeing stabel, check model diagnostics
 

# Residual plots looks good enough
plot(m3, resid(., type = "n") ~ fitted(.)|target)
qqnorm(m3, ~ resid(., type = "n")|target, abline = c(0,1))

# Looks good!


# Using emmeans to get coefficients of interest.

em.m3 <- emmeans(m3, specs = ~ tx|target)
em <- em.m3

# A function to produce contrasts 

contr.fun <- function(em, tx, target) {
  # convert to data frame
  em <- data.frame(em)
  # count number of rows
  rows <- nrow(em)
  # create a vector of zeros
  vec <- rep(0, rows)
  # Specify what indicator should be 1
  vec[which(em$tx == tx & em$target == target)] <- 1
  return(vec)
  
}

c1 <- (contr.fun(em.m3, "post_int_train", "rRNA18S F2R2") - contr.fun(em.m3, "pre_int_train", "rRNA18S F2R2") ) -
  (contr.fun(em.m3, "post_con_train", "rRNA18S F2R2") - contr.fun(em.m3, "pre_con_train", "rRNA18S F2R2") )

contrast(em.m3, 
         method = list("rRNA18S:inter:post_int" = c1)) %>%
  print()



###########################






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






