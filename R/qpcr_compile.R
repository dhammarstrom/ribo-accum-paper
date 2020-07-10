#### qPCR import and data cleaning #######


source("./R/libs.R")


#### 2018-1 import #####


batch <- prepare_batch("./data/wetlab/qPCR/exports/run2018-1", equipment = "quant", skip = 45)


# Filter away targets not used in the present analysis 


targets <- c("UBTF F4R4", 
             "UBTF F6R6",
             "rRNA5.8S F2R2",
             "rRNA28S F2R2",
             "rRNA18S F2R2",
             "rRNA45S F1R1",
             "RPL32 F1R1", 
             "rpS6 F2R2",
             "MyHC1 F1R1",
             "MyHC2A F5R5",
             "MyHC2X F5R5",
             "Lambda KIT")  

batch.filtered <- batch %>%
  filter(target %in% targets) %>%
  print()

# Preliminary models 

models <- model_qpcr(batch.filtered)

# Model tests
model.tests <- test_models(models)


# Plot best model per target
data.frame(table(model.tests$best.model, model.tests$target)) %>%
  ggplot(aes(Var2, Freq, fill = Var1)) + geom_bar(stat="identity") + coord_flip()

# Best model per target are stored for modeling with best model
best.models <- data.frame(target = names(apply(table(model.tests$best.model, model.tests$target), 2, which.max)),
                          model = as.character(row.names(table(model.tests$best.model, model.tests$target))[apply(table(model.tests$best.model, model.tests$target), 2, which.max)]))




## load data with best model
results <- list()

# Loop through all targets in best.models data frame
for(i in 1:nrow(best.models)){
  
  results[[i]] <- batch.filtered %>%
    filter(target == best.models[i,1]) %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>% # use the best model in each model_qpcr
    analyze_models() # analyze models for cpD2
  
}

# combine all results and str split id variables
qpcrdat <- rbind_all(results) 

id.var <- str_split_fixed(qpcrdat$ID, "_", 5) 
colnames(id.var) <- c("subject", "sample", "x", "target", "cDNA")  
qpcrdat <- cbind(id.var, qpcrdat[,-1])



## estimate efficiencies ##
efficiencies <- list()

# use the same loop to analyze efficiencies
for(i in 1:nrow(best.models)){
  
  efficiencies[[i]] <- batch.filtered %>%
    filter(target == best.models[i,1]) %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2])))) %>%
    analyze_efficiency(method = "cpD2", model = "linexp")
  
}

# combine results and use str split to extract id variables
efficiencies <- rbind_all(efficiencies) 
id.var <- str_split_fixed(efficiencies$ID, "_", 4) 
colnames(id.var) <- c("subject", "sample", "x", "target")  
efficiencies <- cbind(id.var, efficiencies[,-1])


efficiencies %>%
  filter(eff > 1.5 & eff < 2.5)%>% # remove outliers from efficiency estimation
  group_by(target)%>%
  summarise(efficiency = mean(eff, na.rm = TRUE),
            max.eff = max(eff, na.rm = TRUE),
            min.eff = min(eff, na.rm = TRUE),
            sd.eff = sd(eff, na.rm = TRUE))%>%
  ggplot(aes(target, efficiency)) + geom_point() + coord_flip()


effs <- efficiencies %>%
  filter(eff > 1.5 & eff < 2.5)%>% # remove outliers from efficiency estimation
  group_by(target)%>%
  summarise(mean.eff = mean(eff, na.rm = TRUE))


### Combine and save data 


qpcrdat1 <- qpcrdat %>%
  dplyr::select(-eff) %>% # remove prelimin efficiencies
  inner_join(effs, by = "target") %>%
  dplyr::select(participant = subject, sample, target, cDNA, cq = cpD2, eff) %>%
  print()



saveRDS(qpcrdat1, file = "./data/derivedData/qpcr/qpcr_run1.RDS")


############### Round 2020-1 ####################





batch <- prepare_batch("./data/wetlab/qPCR/exports/run2020-1", equipment = "quant", skip = 47)


unique(batch$target)

### Change names to correspond to last run ...
batch.filtered <- batch %>%
  mutate(target = if_else(target %in% c("47S F1R1",      "45S F5R5",      "45SITS F12R12", "5.8S F2R2",    
                                        "28S F2R2",      "18S F2R2",      "5S F3R3"), paste0("rRNA", target),
                          if_else(target == "MyHC2a", "MyHC2A", 
                                  if_else(target == "MyHC2x", "MyHC2X", 
                                          if_else(target == "RPS6", "rpS6", 
                                                  if_else(target == "Lambda Kit", "Lambda KIT", target)))))) %>%
  print()




# Preliminary models 

models <- model_qpcr(batch.filtered)

# Model tests
model.tests <- test_models(models)


# Plot best model per target
data.frame(table(model.tests$best.model, model.tests$target)) %>%
  ggplot(aes(Var2, Freq, fill = Var1)) + geom_bar(stat="identity") + coord_flip()

# Best model per target are stored for modeling with best model
best.models <- data.frame(target = names(apply(table(model.tests$best.model, model.tests$target), 2, which.max)),
                          model = as.character(row.names(table(model.tests$best.model, model.tests$target))[apply(table(model.tests$best.model, model.tests$target), 2, which.max)]))




## load data with best model
results <- list()
i <- 1
# Loop through all targets in best.models data frame
for(i in 1:nrow(best.models)){
  
  results[[i]] <- batch.filtered %>%
    filter(target == best.models[i,1]) %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>% # use the best model in each model_qpcr
    analyze_models() # analyze models for cpD2
  
}

# combine all results and str split id variables
qpcrdat <- rbind_all(results) 

id.var <- str_split_fixed(qpcrdat$ID, "_", 5) 
colnames(id.var) <- c("subject", "sample", "x", "target", "cDNA")  
qpcrdat <- cbind(id.var, qpcrdat[,-1])



## estimate efficiencies ##
efficiencies <- list()

# use the same loop to analyze efficiencies
for(i in 1:nrow(best.models)){
  
  efficiencies[[i]] <- batch.filtered %>%
    filter(target == best.models[i,1]) %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>%
    analyze_efficiency(method = "cpD2", model = "linexp")
  
}

temp <- batch.filtered %>%
  filter(target == best.models[i,1]) %>%

  model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>%
  print()
  
  
temp[1:2]


  analyze_efficiency(temp[1:2], method = "outlier", model = "linexp")


temp[[1]]$MODEL


# combine results and use str split to extract id variables

efficiencies.df <- bind_rows(efficiencies)
id.var <- str_split_fixed(efficiencies.df$ID, "_", 5) 
colnames(id.var) <- c("subject", "sample", "x", "target", "cDNA")  
efficiencies.df <- cbind(id.var, efficiencies.df[,-1])


unique(efficiencies.df$target)

efficiencies.df %>%
  filter(eff > 1.5 & eff < 2.5)%>% # remove outliers from efficiency estimation
  group_by(target)%>%
  summarise(efficiency = mean(eff, na.rm = TRUE),
            max.eff = max(eff, na.rm = TRUE),
            min.eff = min(eff, na.rm = TRUE),
            sd.eff = sd(eff, na.rm = TRUE))%>%
  ggplot(aes(target, efficiency)) + geom_point() + coord_flip()


effs <- efficiencies.df %>%
  filter(eff > 1.5 & eff < 2.5)%>% # remove outliers from efficiency estimation
  group_by(target)%>%
  summarise(mean.eff = mean(eff, na.rm = TRUE))


### Combine and save data 


qpcrdat2 <- qpcrdat %>%
  dplyr::select(-eff) %>% # remove prelimin efficiencies
  inner_join(effs, by = "target") %>%
  dplyr::select(participant = subject, sample, target, cDNA, cq = cpD2, eff) %>%
  print()



saveRDS(qpcrdat1, file = "./data/derivedData/qpcr/qpcr_run1.RDS")


