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
qpcrdat <- bind_rows(results) 

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
efficiencies <- bind_rows(efficiencies) 
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
  dplyr::select(participant = subject, sample, target, cDNA, cq = cpD2, eff = mean.eff) %>%
  print()



saveRDS(qpcrdat1, file = "./data/derivedData/qpcr/qpcr_run1.RDS")


############### Round 2020-1 ####################



batch <- prepare_batch("./data/wetlab/qPCR/exports/run2020-1", equipment = "quant", skip = 47)


### Change names to correspond to last run ...
batch.filtered <- batch %>%
  mutate(target = if_else(target %in% c("47S F1R1",      "45S F5R5",      "45SITS F12R12", "5.8S F2R2",    
                                        "28S F2R2",      "18S F2R2",      "5S F3R3"), paste0("rRNA", target),
                          if_else(target == "MyHC2a F5R5", "MyHC2A F5R5", 
                                  if_else(target == "MyHC2x F5R5", "MyHC2X F5R5", 
                                          if_else(target == "MyHC2a F5RF", "MyHC2A F5R5", 
                                          if_else(target == "RPS6 F2R2", "rpS6 F2R2", 
                                                  if_else(target == "Lambda Kit", "Lambda KIT", target))))))) %>%
  filter(replicate %in% c("cDNA1", "cDNA2")) %>%
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
qpcrdat <- bind_rows(results) 

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


# combine results and use str split to extract id variables

efficiencies.df <- bind_rows(efficiencies)
id.var <- str_split_fixed(efficiencies.df$ID, "_", 5) 
colnames(id.var) <- c("subject", "sample", "x", "target", "cDNA")  
efficiencies.df <- cbind(id.var, efficiencies.df[,-1])


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

# Borrow efficiencies from prevous runs and avarege...
effs <- qpcrdat1 %>%
  group_by(target) %>%
  summarise(mean.eff = mean(eff, na.rm = TRUE )) %>%
  rbind(effs) %>%
  group_by(target) %>%
  summarise(mean.eff = mean(mean.eff, na.rm = TRUE )) %>%
  print()
 

### Combine and save data 


qpcrdat2 <- qpcrdat %>%
  dplyr::select(-eff) %>% # remove prelimin efficiencies
  inner_join(effs, by = "target") %>%
  dplyr::select(participant = subject, sample, target, cDNA, cq = cpD2, eff = mean.eff) %>%
  print()



saveRDS(qpcrdat2, file = "./data/derivedData/qpcr/qpcr_run2.RDS")


###### Clean qpcr data and insert sample information ################

# Notes:
# Two qPCR rounds have been performed and amplifications estimates are
# stored in "./data/derivedData/qpcr/qpcr_run1.RDS" and
# "./data/derivedData/qpcr/qpcr_run2.RDS". To get muscle weights for normalization
# use "./data/tr010_mRNASamples.xlsx" for round 1.




qpcrdat1 <- readRDS(file = "./data/derivedData/qpcr/qpcr_run1.RDS")
qpcrdat2 <- readRDS(file = "./data/derivedData/qpcr/qpcr_run2.RDS")




### Samples are loaded from extractions setup
samples <-  read_excel("./data/tr010_mRNASamples.xlsx", na = "NA") %>%
  filter(IncludeTotalRNA == "YES") %>%
  mutate(sample = paste("S", ExtractionNR, sep = "")) %>%
  dplyr::select(participant, leg, time, sample, tissue_weight) %>%
  print()



### Combine data sets, qpcr, samples and leg randomizations
qdat1 <- qpcrdat1 %>%
  inner_join(samples) %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  dplyr::select(participant, leg, time, sex, cond, target, cq, eff, cdna = cDNA, tissue_weight) %>%
  print()


################## qPCR dat 2 ###############



qdat2 <- qpcrdat2 %>%
  filter(cDNA %in% c("cDNA1", "cDNA2")) %>%
  mutate(sample = gsub("S", "", sample), 
         series = gsub("cDNA", "", cDNA)) %>%
  unite(col = "sample_series", c(sample, series)) %>%
  
  mutate(sample_series = if_else(participant == "P19" & sample_series == "1_1", 
                                 "7_1", 
                                 if_else(participant == "P22" & sample_series == "1_2", 
                                         "3_1", sample_series))) %>%
  separate(sample_series, into = c("sample", "series"), convert = TRUE) %>% 
  dplyr::select(participant, sample, series, target, cq, eff) 



  ### Join sample information 
  # Subject 19, 21-23  
 sample_info <-   rbind(read.csv("./data/written/extraction_numbers_round2.csv", sep = ";") %>%
      separate(extraction_nr,
               into = c("series", "sample"),
               convert = TRUE) %>%
      mutate(series = if_else(series == "I", 1, 2)) %>%
      dplyr::select(participant, leg, time, tissue_weight, series, sample) %>%
      filter(!(is.na(tissue_weight))),
      # Add information on subject 18 (page 181 in lab note book for sample setup)
        data.frame(
          Sample = c("P18-1", "P18-2", "P18-3", "P18-4", "P18-5", "P18-6"),
          cDNA = c("mRNA1", "mRNA1", "mRNA2", "mRNA2", "mRNA2", "mRNA1"),
          time = c("postctrl", "S0", "S0", "postctrl", "S0", "S0"),
          leg = c("L", "R", "L", "L", "R", "L")
        ) %>%
          inner_join(
            read_excel("./data/tissue/tr010_tissue.xlsx", na = "NA") %>%
              filter(participant == "P18", sample != "prot") %>%
              
              dplyr::select(participant, leg, time, cDNA = sample, tissue_weight)
          ) %>%
          mutate(sample = gsub("P18-", "", Sample),
                 series = 1) %>%
          dplyr::select(participant, leg, time, tissue_weight, series, sample)
      # Add participants from round 1 (duplicate series to get cDNA1 and 2)
      ,samples %>%
        mutate(sample = gsub("S", "", sample), 
               series = 1) %>%
        dplyr::select(participant, leg, time, tissue_weight, series, sample),
      # Participants from round 1, cDNA 2 (duplicate information)
      samples %>%
        mutate(sample = gsub("S", "", sample), 
               series = 2) %>%
        dplyr::select(participant, leg, time, tissue_weight, series, sample)) %>%
      mutate(series = as.character(series)) %>%  # fix data type for compability when inner_join
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>% # leg randomization file
  dplyr::select(participant, leg, time, sample, series, sex, cond, tissue_weight) %>%
  print()

str(sample_info)
str(qdat2)

qdat2_complete <-  qdat2 %>%
  mutate(sample = as.character(sample), 
         series = as.character(series)) %>%
   inner_join(sample_info) %>%
  mutate(cdna = paste0("cDNA", series)) %>%
   dplyr::select(participant, leg, time, sex, cond, target, cq, eff, cdna, tissue_weight) %>%
  print()



qdat <- rbind(qdat1, qdat2_complete) %>%
  print()

############## Save qdat for analysis ##############################
saveRDS(qdat, "./data/derivedData/qpcr/qpcr_compiled.RDS")
















