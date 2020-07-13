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




  
temp <- rna_complete %>%
  filter(cond != "ctrl_leg") %>%
  filter(time != "post1w", 
         rna/tissue_weight < 1200) %>%
  
  group_by(participant,leg, time,time.c, cond, tissue_weight) %>%
  summarise(rna = mean(rna, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(participant = factor(participant), 
         cond = factor(cond), 
         rna.tissue = log(rna/tissue_weight), 
         tw = tissue_weight - mean(tissue_weight), 
         time.c.cent = time.c - 4) %>%
  print()
 
library(gamm4); library(gamm4.test); 
library(ggeffects); library(brms)









m0.gamm <- gamm4(rna.tissue ~  s(time.c, k = 7), 
                 random = ~ (1|participant), 
                 data = temp)

m1.gamm <- gamm4(rna.tissue ~ s(time.c.cent, k = 7, by = cond), 
            random = ~ (1|participant), 
            data = temp)



f1 <- lme(rna ~ tw + time * cond, random = list(participant = ~ 1, leg = ~ 1), 
    data = temp)


plot(f1)



m1 <- gam(rna ~  tw + s(time.c, by = cond, k = 7) +  s(participant, bs = "re"), 
          data = temp)

m2 <- gam(rna ~ tw +  s(time.c, by = cond, k = 7) +  s(participant, bs = "re"), 
          data = temp, method = "REML")


m1.brm <- brm(bf(rna ~ tw +  s(time.c, by = cond, k = 7) +  s(participant, bs = "re")), 
              data = temp)




anova(m1, m2, test = "Chisq")
par(mfrow = c(1, 1))
gam.check(m2)



plot(m2, pages = 1, unconditional = FALSE)

pr <- ggpredict(m2, c("time.c", "cond"), type = "fe")

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
  
  
  
  
