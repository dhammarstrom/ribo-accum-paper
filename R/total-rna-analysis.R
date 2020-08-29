######## Total RNA analysis ########################

# Load libraries and data 
source("./R/libs.R")

rna <- readRDS("./data/derivedData/tot-rna/tot-rna.RDS")



### 


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
                                        "S9", "S12", "post1w", "postctrl"))) %>%
  print()






# Time course data in experimental group #

  
time_course <- rna_complete %>%
  filter(cond != "ctrl_leg") %>%
  filter(time != "post1w") %>%
  
  group_by(participant,leg, time,time.c, cond, ) %>%
  summarise(rna = mean(rna, na.rm = TRUE), 
            tissue_weight = mean(tissue_weight, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(participant = factor(participant), 
         cond = factor(cond), 
         rna.tissue = log(rna/tissue_weight), 
         tw = tissue_weight - mean(tissue_weight), 
         time.c.cent = time.c - 4) %>%
  print()
 
# Exploratory plot for each participant

time_course %>%
  ggplot(aes(time.c, rna.tissue, 
             color = cond,
             group = paste0(participant, leg))) + geom_line() +
  facet_wrap(~ participant)


# Polynomial effects -- is the relation curve linear? 
f0 <- lme(rna ~ tw + time.c * cond, random = list(participant = ~ 1, leg = ~ 1), 
          data = time_course, method = "ML")

f1 <- lme(rna ~ tw + poly(time.c, 2) * cond, random = list(participant = ~ 1, leg = ~ 1), 
    data = time_course, method = "ML")

f2 <- lme(rna ~ tw + poly(time.c, 3) * cond, random = list(participant = ~ 1, leg = ~ 1), 
          data = time_course, method = "ML")


anova(f0, f1, f2)
# Evidence for curvelinear relationship but 3rd degree does not fit better than 2nd. 

# Predicted values divided by mean tissue weight to get RNA / tissue weight (ng/mg)
pr <- ggpredict(f2, c("time.c [all]", "cond"), type = "fe") %>%
  data.frame() %>%
  mutate(predicted = predicted / mean(time_course$tissue_weight), 
         conf.low = conf.low / mean(time_course$tissue_weight), 
        conf.high = conf.high / mean(time_course$tissue_weight)) %>%
  print()



# For plotting th e
pr %>%
  ggplot(aes(x, predicted, color = group)) + geom_line() 
  

# More accurate predictions with gam models

# Is there an effect of splines (curvature)?
# Effect of condition in curvature?

# No splines
m0 <- gam(rna ~ tw + cond + time.c +  s(participant, bs = "re"), 
          data = time_course, method = "ML")
# Splines but not per group
m1 <- gam(rna ~ tw + cond +  s(time.c,  k = 7) +  s(participant, bs = "re"), 
          data = time_course, method = "ML")
# Splines per group
m2 <- gam(rna ~ tw + cond + s(time.c, by = cond,  k = 7) +  s(participant, bs = "re"), 
          data = time_course, method = "ML")

# Splines improves fit (as with polynomials above)
anova(m0, m1, m2, test = "Chisq")

# AIC does not improve when splines are fitted per group -- indication of no "interaction"
AIC(m1, m2)

# REML estimation full model
m3 <- gam(rna ~ tw +  s(time.c, k = 7,  by = cond) + cond +  s(participant, bs = "re"), 
          data = time_course, method = "REML", select = TRUE)


library(ggeffects)
## Prediction from the full model 
pr <- ggpredict(m3, c("time.c", "cond"), type = "re")

str(pr)

# Prediction from the reduced model 
pr <- ggpredict(m2, c("time.c"), type = "re")

plot(pr)

summary(m3)

dev.off()
?s
anova(m0, m1, test = "F")

summary(m2)


predict(m1, newdata = data.frame(time.c = c(0, seq(1:12), 0, seq(1:12)), 
                                 cond = rep(c("var", "const"), each = 13)), 
        type = "response")



plot(m1, pages = 1)
gam.check(m1)




m0.brm <- brm(rna ~ tw +  s(time.c, k = 7) +  s(participant, bs = "re"),
              control = list(adapt_delta = 0.99),
              cores = 16,
              chains = 4,
              data = temp)
m1.brm <- brm(rna ~ tw +  s(time.c, by = cond, k = 7) +  s(participant, bs = "re"),
              control = list(adapt_delta = 0.99),
              cores = 16,
              chains = 4,
              data = temp)




AIC(m0, m1)



summary(m1)
summary(m1.brm)

pp_check(m1.brm)

plot(conditional_effects(m1.brm), points = TRUE)


mx <- brm(rna ~ tw + time*cond +  (1|participant), 
          chains = 2, silent = FALSE,
          cores = 12, data = temp)


loo(m0.brm, m1.brm)


library(brms)
bprior1 <- prior(student_t(5,0,10), class = b) +
  prior(cauchy(0,2), class = sd)
fit1 <- brm(count ~ zAge + zBase * Trt + (1|patient),
            chains = 1,
            data = epilepsy, family = poisson(), prior = bprior1)

summary(fit1)


anova(m1, m2, test = "Chisq")
par(mfrow = c(1, 1))
gam.check(m2)



plot(m1, pages = 1, unconditional = FALSE)

pr <- ggpredict(m2, c("time.c", "cond"), type = "re")

plot(pr)


data.frame( ggpredict(m1, c("time.c", "cond"), type = "fe")) %>%
    ggplot(aes(x, predicted, color = group)) + geom_line() + geom_point()
  
data.frame( ggpredict(m2, c("time.c", "cond"), type = "re")) %>%
  ggplot(aes(x, predicted, color = group)) + geom_line() + geom_point() 



?ggpredict

anova(m0, m1, test = "Chisq")

summary(m1)



test <- gamm4.grptest(rna.tissue ~ s(time.c, k = 7, by = cond), 
                    random =  ~ (1|participant), 
                    test = ~ cond,
                    data = temp)


anova(m0.gamm$mer, m1.gamm$mer)
plot(m1.gamm$mer)


temp$pred.gamm <- predict(m1.gamm$gam)


temp %>%
  ggplot(aes(time.c, pred.gamm, color = cond, group = paste(participant, cond))) + 

  geom_line()


temp %>%
  filter(time != "post1w", 
         rna/tissue_weight < 1200) %>%
  ggplot(aes(time.c, rna/tissue_weight, color = cond)) +
  geom_point() + 
  geom_smooth()
  geom_text_repel(data = . %>%
                    filter(rna/tissue_weight > 1200), 
                  aes(label = paste(participant, sample)))
  
  
  ############### Time as factor expr vs. control #################

  temp <- rna_complete %>%
    filter(time %in% c("S0", "S1", "S1c",  "S12", "postctrl", "post1w")) %>%
    mutate(Time = if_else(time == "S1c", "S1", 
                          if_else(time %in% c("post1w", "postctrl", "S12"), 
                                  "post", as.character(time))), 
           tw = tissue_weight - mean(tissue_weight), 
           Time = factor(Time, levels = c("S0", "S1", "post"))) %>%
   # group_by(participant, Time, detrain, cond) %>%
   # summarise(rna = mean(rna), 
   #           tw = mean(tw)) %>%
    mutate(Cond = paste(Time, cond, detrain, sep = ":"), 
           participant = factor(participant),
           detrain = factor(detrain, levels = c("train", "detrain")), 
           sample = paste0(participant, leg, time)) %>%
    filter(time != "S12") %>%
    mutate(grp = if_else(cond == "ctrl_leg", "ctrl", "expr")) %>%

    print()
    
    
  m1 <- lme(log(rna) ~ tw + Time * grp, 
            random = list(participant = ~ 1, leg = ~ 1, sample = ~ 1), 
            data = temp)  
  
  
  summary(m1)
  
  
  mx <- gam(rna ~ tw + Time * cond + detrain + detrain:cond + s(participant, bs = "re"), 
            method = "REML", data = temp)
  
  summary(mx)
  
  plot(m1)
  summary(m1)
  

  temp %>%
    filter(cond == "ctrl_leg") %>%
    group_by(participant, time) %>%
    summarise(rna = mean(rna), 
              tissue_weight = mean(tissue_weight)) %>%
    ggplot(aes(time, rna/tissue_weight, group = paste(participant), 
               label = participant)) + geom_line() + geom_text() + 
    geom_point()
  
  
  
  summary(m1)
  plot(m1)
  
  emmeans(m1, specs = ~ "Time|grp") %>%
    data.frame() %>%

    mutate(time = factor(Time, levels = c("S0", "S1", "post"))) %>%
    ggplot(aes(time, exp(emmean)/ mean(temp$tissue_weight), shape = grp)) + 
    
    
    geom_errorbar(aes(ymin = exp(lower.CL)/ mean(temp$tissue_weight),  ymax = exp(upper.CL)/ mean(temp$tissue_weight)), 
                  width = 0.15, 
                  position = position_dodge(width = 0.1)) +
    geom_point(size = 3, position = position_dodge(width = 0.1)) + 
    
    scale_shape_manual(values = c(21,22))
  
  
  filter(cond != "ctrl_leg") %>%
    filter(time != "post1w") %>%
    
    group_by(participant,leg, time,time.c, cond, ) %>%
    summarise(rna = mean(rna, na.rm = TRUE), 
              tissue_weight = mean(tissue_weight, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(participant = factor(participant), 
           cond = factor(cond), 
           rna.tissue = log(rna/tissue_weight), 
           tw = tissue_weight - mean(tissue_weight), 
           time.c.cent = time.c - 4) %>%
    print()
  
  
  
  
  
