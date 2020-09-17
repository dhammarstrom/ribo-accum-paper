# Training load data 


library(tidyverse)
library(readxl)
library(brms)
library(lme4)



load.stats <- read_excel("./data/tr010_training.xlsx", sheet = 1, na = "NA") %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  filter(exercise == "legext") %>%
  mutate(set.load = repetitions * load, 
         week = if_else(session %in% c(1:4), "W1", if_else(session %in% c(5:8), "W2", "W3"))) %>%
  group_by(participant,leg,sex, cond, session, week) %>%
  summarise(load = mean(load, na.rm = TRUE)) %>%
  mutate(week = factor(week, levels = c("W1", "W2", "W3")))  %>%
  print()

  

load.stats %>%
  ggplot(aes(week, load, color = sex, shape = cond)) + geom_point(position = position_jitter(width = 0.2))


# Is there differences between conditions?



m1 <- brm(bf(load ~ week * cond + (1|participant)),
               
               warmup = 1000, # number of samples before sampling
               iter = 4000,  # number of mcmc iterations 
               cores = 4, # number of cores used in sampling
               chains = 4, # number of chains
               seed = 5, # seed for reproducible results
               #    control = list(adapt_delta = 0.95), 
               data = load.stats)

pp_check(m1, type = "stat", stat = "median")
# The data is in agreement with what the model predicts, good.
# Leave one out statistics
loo(m1)

# No effects of condition, these are disregarded to produce 

m2 <- brm(bf(load ~ week + (1|participant/leg)),
          
          warmup = 1000, # number of samples before sampling
          iter = 4000,  # number of mcmc iterations 
          cores = 4, # number of cores used in sampling
          chains = 4, # number of chains
          seed = 5, # seed for reproducible results
          #    control = list(adapt_delta = 0.95), 
          data = load.stats)

pp_check(m2, type = "stat", stat = "median")
# The data is in agreement with what the model predicts, good.
# Leave one out statistics
loo(m1)

summary(m2)
pp_check(m2, type = "ecdf_overlay")




saveRDS(m2, "./data/derivedData/training-data/m2.RDS")
saveRDS(m1, "./data/derivedData/training-data/m1.RDS")




# All Pareto estimates are OK
