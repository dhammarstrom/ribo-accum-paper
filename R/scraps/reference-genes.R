#### Reference gene stability in time-course data set ##########

# library(devtools)
# install_github("dhammarstrom/generefer")

library(tidyverse)
library(generefer)


# Load data
qdat <- readRDS("./data/derivedData/qpcr/qpcr_compiled.RDS")



qdat2 <- qdat %>%
  group_by(target) %>%
  mutate(outlier = if_else(cq > mean(cq, na.rm = TRUE) + 3 * sd(cq, na.rm = TRUE)|
                             cq < mean(cq, na.rm = TRUE) - 3 * sd(cq, na.rm = TRUE), 
                           "outlier", "in")) %>%
  filter(outlier == "in", target != "Lambda KIT") %>%
  group_by(participant, leg, time, sex, cond, cdna, target) %>%
  summarise(cq = mean(cq, na.rm = TRUE), 
            eff = mean(eff, na.rm = TRUE), 
            tissue_weight = mean(tissue_weight, na.rm = TRUE)) %>%
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

  print()


qdat2$target

combinations <- combn(unique(data.frame(qdat2)[,7 ]), 2)







temp <- qdat2 %>%
  filter(cond != "ctrl_leg", 
         time != "post1w", 
         !(target %in% c("rRNA45S F1R1", 
                         "rRNA45S F5R5", 
                         "rRNA45SITS F12R12", 
                         "rRNA47S F1R1"))) %>%
  mixedmodel_stability(target = "target", 
                       response = "Ra", 
                       fixed.effects = "target + time * cond", 
                       random.effect = "(1|participant)",
                       form = "Ra ~ target + time * cond + (1|participant)", 
                       reduced.model = "Ra ~ target + (1|participant)", 
                       icc.model = NULL,
                       n.genes = c(2, 3),
                       hypothesis.test = "LRT")



df <- qdat2 %>%
  filter(target %in% c("rRNA28S F2R2",
                       "rRNA5.8S F2R2", 
                       "rRNA5 F3R3"),
         cond != "ctrl_leg", 
         time != "post1w") %>%
  print()


m <- lme4::lmer(Ra ~ target + time * cond + (1|participant), 
     data = df)

summary(m)
m0 <- lme4::lmer(Ra ~ target + (1|participant), 
                 data = df)
anova(m, m0)
summary(m)


emmeans::emmeans(m, specs = ~ cond|time) %>%
  data.frame() %>%
  ggplot(aes(time, emmean, color = cond, group = cond)) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL)) +
  geom_line()



icc.calc <- function(model){
  icc <- as.numeric(data.frame(lme4::VarCorr(model))[1,4]) / sum(as.numeric(data.frame(lme4::VarCorr(model))[,4]))
  return(icc)
}

icc.calc(m0)

lme4::VarCorr(m0)


temp %>%
  filter(!is.na(icc)) %>%
    mutate(icc.l = round(icc.l, 2)) %>%
  arrange(-icc.l) %>%
  print()


