
### Leg lean mass and muscle thickness ##############





# Load libraries
source("R/libs.R")



# Leg randomization data ##############

leg <- read_excel("./data/leg_randomization.xlsx")


# DXA data ################

# This file has a modified ROI using only the thigh muscle
dxa_modified_roi <- read_excel("./data/tr010_dxa.xlsx") %>%
  dplyr::select(participant, group, time, leg, lean_g) %>%
  print()



# Full data set
dxa <- read_excel("./data/tr010_dxa_standarddefinitions.xlsx") %>%
  dplyr::select(Etternavn, Kjønn,Målingsdato,Alder, Vekt:`HeleVev%fett`) %>%
  pivot_longer(names_to = "variable", values_to = "value", cols = Vekt:`HeleVev%fett`) %>%
  mutate(value = gsub("kg", "", value), 
         value = gsub("[[:space:]]", "", value, fixed = FALSE), 
         value = as.numeric(value)) %>%
  filter(variable %in% c("BeinhøyreMagermasse", "BeinvenstreMagermasse")) %>%
  separate(Etternavn, into = c("Study", "participant")) %>%
  dplyr::select(participant, Date = Målingsdato, variable, value) %>%
  filter(!(participant %in% c("P08", "P20", "P17"))) %>%
  group_by(participant) %>%
  mutate(time = if_else(Date == min(Date), "S0", 
                        if_else(Date == max(Date), "post1w", 
                                "S12")), 
         leg = if_else(variable == "BeinhøyreMagermasse", "R", "L")) %>%
  dplyr::select(participant, Date, time, leg, lean.mass = value) %>%
  inner_join(leg) %>%
  mutate(cond = if_else(cond == "ctrl_leg", "ctrl", cond), 
         time = factor(time, levels = c("S0", "S12", "post1w")), 
         group = if_else(cond == "ctrl", "ctrl", "exp")) %>%
  dplyr::select(-Date) %>%
  inner_join(dxa_modified_roi) %>%
  mutate(lean.mean = (lean.mass + lean_g)/2, 
         group = factor(group, levels = c("ctrl", "exp")), 
         time = factor(time, levels = c("S12", "S0", "post1w"))) %>%
#  filter( group != "ctrl") %>%
  print()








# Ultra sound data ########################

### Load data ####

### Fix this

us_data <- read_excel("./data/ultrasound/ultrasound_data.xlsx") %>%
  inner_join(read_csv("./data/ultrasound/ultrasound_codekey.csv")) %>%
  mutate(leg = gsub("VL", "", leg)) %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  dplyr::select(participant, time, leg, sex, cond, code, length) %>%
  mutate(group = if_else(participant %in% paste("P", 1:7, sep = ""), "experiment", "control")) %>%
  group_by(participant, time, leg, sex, cond, group) %>%
  summarise(thickness = mean(length, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = c("pre", "post", "post1w")),
         time.pp = factor(gsub("1w", "", time), levels = c("pre", "post")),
         detrain = if_else(time == "post1w", "detrain", "train")) %>%
  print()





us_data_2019 <- read_excel("./data/ultrasound/ultrasound_data_2019.xlsx") %>%
  inner_join(read_csv("./data/ultrasound/ultrasound_codekey_2019.csv")) %>%
  mutate(leg = gsub("VL", "", leg)) %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  dplyr::select(participant, time, leg, sex, cond, code, length) %>%
  mutate(group = if_else(participant %in% paste("P", c(1:7, 19:23), sep = ""), "experiment", "control")) %>%
  group_by(participant, time, leg, sex, cond, group) %>%
  summarise(thickness = mean(length, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = c("pre", "post", "post1w")),
         time.pp = factor(gsub("1w", "", time), levels = c("pre", "post")),
         detrain = if_else(time == "post1w", "detrain", "train")) %>%
  print()




#### Mean center pre 

us_temp <- us_data %>%
  rbind(us_data_2019) %>%
  dplyr::select(participant, time, leg, sex, cond, group, thickness) %>%

  pivot_wider(names_from = time, values_from = thickness) %>%
  
  mutate(post   = post -  pre, 
         post1w = post1w -  pre, 
         pre    = pre - mean(pre)) %>%
  
  pivot_longer(names_to = "time", values_to = "change", post:post1w) %>%
  filter(!is.na(change)) %>%
  mutate(cond = if_else(cond == "ctrl_leg", "ctrl", cond),
         group_time = paste0(group, "_", time), 
         group_time = factor(group_time, levels = c("control_post", 
                                                    "experiment_post", 
                                                    "experiment_post1w"))) %>%
   #      group_time = factor(group_time, levels = c("ctrl_post", 
   #                                                 "var_post", 
   #                                                 "const_post", 
   #                                                 "var_post1w", 
   #                                                 "const_post1w"))) %>%
  print()
 


  




## CONST vs. VAR 

m1 <- lme(change ~  sex + pre + group_time,
          random = list(participant = ~ 1),
          na.action = na.omit,
          control = list(msMaxIter = 120, 
                         opt = "nloptwrap", 
                         msVerbose = TRUE), 
          data = us_temp)


marginal.means <- emmeans(m1, specs = ~group_time) %>%
  data.frame()




raw_scores <- us_temp %>%
  ggplot(aes(group_time, change)) + 
  
  geom_hline(yintercept = 0, color = "gray90") +
  
  geom_jitter(width = 0.05) + 
  geom_errorbar(data = marginal.means, 
                aes(group_time, emmean, ymin = lower.CL, ymax = upper.CL), 
                position = position_nudge(x = 0.1), width = 0) +
  geom_point(data = marginal.means, 
             aes(group_time, emmean), 
             position = position_nudge(x = 0.1), color = "white", size = 0.5, shape = 15) + 
  
  labs(y = "\U0394 Muscle thickness (mm)") +
  
  plot_theme() +
  theme(axis.line.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank())
  
# scale_y_continuous(limits = c(-8, 8), 
#                    breaks = c(-2.5, 0, 2.5, 5), 
#                    labels = c(-2.5,   0,  2.5, 5), 
#                    expand = c(0, 0))



comp_panel <- intervals(m1)$fixed %>%
  data.frame() %>%
  mutate(coef = rownames(.)) %>%
  filter(coef %in% c("group_timeexperiment_post", "group_timeexperiment_post1w")) %>%

  mutate(group_time = gsub("group_time", "", coef)) %>%
  dplyr::select(-coef) %>%
  rbind(data.frame(lower = NA, est. = NA, upper = NA, group_time = "control_post")) %>%
  
  mutate(group_time = factor(group_time,  
                             levels = c("control_post", 
                                       "experiment_post", 
                                       "experiment_post1w"), 
                             labels = c("CON", "EXP S12", "EXP De-train"))) %>%

  ggplot(aes(group_time, est.)) + 
  
  geom_hline(yintercept = 0, color = "gray90") +
  
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) + 
  geom_point(size = 2) + 

  scale_y_continuous(limits = c(-0.5, 2.5), 
                     breaks = c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5), 
                     labels = c("",   0,  "", 1,  "", 2,  ""), 
                     expand = c(0, 0)) +
  
  labs(y = "\U0394 EXP \U2212 \U0394 CON") +
  
  plot_theme() +
  
  theme(axis.title.x = element_blank())
  
  




plot_grid(raw_scores, comp_panel, nrow = 2, align = "v")



 


emmeans(m1, specs = ~group_time) %>%
  data.frame() %>%
  ggplot(aes(group_time, emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, 
                position = position_dodge(width = 0.2)) +
  geom_point(shape = 21, size = 3, position = position_dodge(width = 0.2))




spread(time, thickness) %>%
  mutate(post = ((post/pre)-1)*100,
         post1w = ((post1w/pre)-1)*100) %>%
  # group_by(sex) %>%
  mutate(pre = pre - mean(pre, na.rm = TRUE)) %>%
  gather(timepoint, change, post:post1w) %>%
  mutate(cond_time = paste(group, timepoint)) %>%
  filter(!is.na(change)) %>%
  mutate(cond_time = factor(cond_time, 
                            levels = c("control post", 
                                       "experiment post", 
                                       
                                       "experiment post1w"))) %>%
  print()


us_temp %>%
  ggplot(aes(participant, change, color = cond, shape = timepoint)) + geom_point(size = 2) +
  xlab("Participant") + ylab("Percentage change")

m1 <- lme(change ~ cond_time,
          random = list(participant = ~1, leg = ~1),
          na.action = na.exclude,
          data = us_temp)


plot(m1)


us_em <- emmeans(ref_grid(m1), specs = ~cond_time)


summary(m1)

#### Fix from here #### 

pos <- position_dodge(width = 0.2)



us_em %>%
  data.frame() %>%
  separate(cond_time, into = c("cond", "time"), sep = " ") %>%
  ggplot(aes(time, emmean, color = cond)) + 
  geom_point(position = pos) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                position = pos, 
                width = 0.2) 
geom_point(data = us_temp, aes())
print()


summary(m1)

us_temp %>%
  group_by(timepoint, cond) %>%
  
  ggplot(aes(timepoint, change, color = cond, group = cond)) + 
  geom_point(position = position_dodge(width = 0.2)) + 
  geom_label(aes(label = participant), position = position_jitter())
print()




filter(group == "experiment") %>%
  ggplot(aes(time, thickness, group = paste(leg, participant), color = participant)) + geom_line() 


us_data_exp <- us_data %>%
  mutate()
print()


m1 <- lmer(thickness ~ time *  cond  + (1|participant/leg), data = us_data_exp)
summary(m1)

drop1(m1, test ="Chi")



us_data_change <- us_data %>%
  dplyr::select(-detrain, -time.pp) %>%
  spread(time, thickness) %>%
  mutate(post = post - pre,
         post1w = post1w - pre) %>%
  gather(time, change, post:post1w) %>%
  print()

m1 <- lmer(change ~ pre + time * group + (1|participant), data = us_data_change)

summary(m1)

### Experimental group ###


us_data_change <- us_data %>%
  filter(group == "experiment") %>%
  dplyr::select(-detrain, -time.pp) %>%
  spread(time, thickness) %>%
  mutate(post = post - pre,
         post1w = post1w - pre) %>%
  gather(time, change, post:post1w) %>%
  print()

m2 <- lmer(change ~ pre + time * cond + (1|participant), data = us_data_change)   

drop1(m2, test = "Chisq")


summary(m2)

















