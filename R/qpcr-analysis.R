######################### qPCR analysis #############################


source("./R/libs.R")

# qdat is the compiled data frame containing sample information and 
# qpcr estimates. Created in qpcr-compile.R. Use this for downstream analyses

qdat <- readRDS("./data/derivedData/qpcr/qpcr_compiled.RDS")



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
  print()
  

# Remove bad reaction prior to modelling. 
# Bad reactions are no amplification or estimatyed to cq < 5
 
unique(qdat.rrna$target)


qdat.rrna %>%
     ggplot(aes(cq, color = target)) + geom_histogram() + 
  facet_wrap(~ target)

qdat.rrna %>%
  group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 5 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 5 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  
  ggplot(aes(technical, cq, color = target, shape = outlier)) + 
  geom_point() + theme_minimal() 

# 5 sd's away from target mean seems to capture the worst reactions.. 


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
intervals(m2)

plot(m2, resid(., type = "n") ~ fitted(.)|target)
qqnorm(m2, ~ resid(., type = "n")|target, abline = c(0,1))

intervals(m1)


rna.interaction <- data.frame(coef(summary(m2))) %>%
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






