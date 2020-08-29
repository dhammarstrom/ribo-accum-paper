########################  Total RNA from replicates ############################



# Total RNA estimates are loaded from replicates. Sample weights, participant, leg, and time data 
# are retrieved from different data sources and a single file saved containing rna amounts 
# per sample. Use this for modeling. 

# Round 2 data: participants P18 (control), P19, P21-23.
# Round 1 data all other participants.


source("./R/libs.R")



######## Round 2 data ################################################

files <- list.files("./data/tot_rna/total-rna-replicates/")



results <- list()

for(i in 1:length(files)){
  
  
  results[[i]] <-  read_excel(paste0("./data/tot_rna/total-rna-replicates/", files[i]), skip = 0, na = "NaN")
  
  
}


tot_rna <- bind_rows(results) %>%
  dplyr::select(Sample, dilution = `Dilution factor`, concentration = ...7) %>%
  mutate(concentration = concentration * dilution) %>%
  dplyr::select(-dilution)



# Extraction numbers round 2

# Filter P18 and code samples
# Elution volume is 20 ul in P18.
p18_samples <- tot_rna %>%
  filter(Sample %in% paste0("P18-", rep(1:6))) %>%
  inner_join(data.frame(Sample = c("P18-1", "P18-2", "P18-3", "P18-4", "P18-5", "P18-6"), 
                        sample = c("mRNA1","mRNA1", "mRNA2", "mRNA2", "mRNA2", "mRNA1" ),
                        time = c("postctrl", "S0", "S0", "postctrl", "S0", "S0"), 
                        leg = c("L", "R", "L", "L", "R", "L"))) %>%
  mutate(concentration = ((concentration/3)/(50*(10/0.51))) * (40 *(10/0.51)) ,
         rna = concentration * 20) %>%
  inner_join(read_excel("./data/tissue/tr010_tissue.xlsx", na = "NA") %>%
               filter(participant == "P18", sample != "prot")) %>%
  mutate(series = 1) %>%
  dplyr::select(participant, series, sample, leg, time, tissue_weight, rna) %>%
  print()
 


# Elution volume for P19, 21, 22, 23: 25 ul. 
tot_rna_round2 <- tot_rna %>%
  filter(!(Sample %in% paste0("P18-", rep(1:6)))) %>%
  separate(Sample, into = c("participant", "series", "sample"), convert = TRUE) %>%
  inner_join(read.csv("./data/written/extraction_numbers_round2.csv", sep = ";") %>%
               separate(extraction_nr, into = c("series", "sample"), convert = TRUE) %>%
               mutate(series = if_else(series == "I", 1, 2))) %>%
  
  # Change concentration as this was estimated with the wrong factor in nanodrop
  # Calculations:
  # RNA = Absorbance * 40 * (10/0.51) * dilutionfactor (= 3)
  # [The above was mistakenly calculated as DNA in raw data:
  # DNA = Absorbance * 50 * (10/0.51)]
  mutate(concentration = ((concentration/3)/(50*(10/0.51))) * (40 *(10/0.51)) ,
         rna = concentration * 25) %>%
  dplyr::select(participant, series, sample, leg, time, tissue_weight, rna) %>%
  print()


tot_rna_round2 <- rbind(tot_rna_round2, p18_samples)



#################### Round 1 data ###############################################


sample_setup <- read_excel("./data/tr010_mRNASamples.xlsx", na = "NA") %>%
  inner_join(read_excel("./data/tissue/tr010_tissue.xlsx", na = "NA") %>%
               filter(sample == "mRNA") %>%
               mutate(date = as.Date(as.numeric(date), origin = "1899-12-30")) %>%
               dplyr::select(participant, leg, time, samplenr, date)) %>%
  filter(IncludeTotalRNA == "YES") %>%
  dplyr::select(participant, leg, time, ExtractionNR, tissue_weight, elution) %>%
  print()



### Read files from Total RNA measurements
files <- list.files("./data/tot_rna/total-rna-raw-round1/")

results <- list()

for(i in 1:length(files)) {
  
  results[[i]] <- read_excel(paste0("./data/tot_rna/total-rna-raw-round1/", files[i], sep = ""), na = "NA")  
  
}

tot_rna_round1 <- bind_rows(results) %>%
  filter(Sample != "Blank1") %>%
  separate(Sample, c("participant", "ExtractionNR")) %>%
  dplyr::select(participant, ExtractionNR, TotalRNA) %>%
  mutate(ExtractionNR = as.numeric(ExtractionNR)) %>%
  inner_join(sample_setup) %>%
  mutate(rna = TotalRNA * elution, 
         series = 1) %>%
  dplyr::select(participant, series, sample = ExtractionNR, leg, time, tissue_weight, rna) %>%
  print()
  
  


############### Combined data set ################################################


rna <- rbind(tot_rna_round1, tot_rna_round2) %>%
  filter(!is.na(tissue_weight)) %>%
  print()




#### Remove outliers based on weight to rna relationship ###
# model for RNA concentration to weight relationship
l.mod <- lm(log(rna) ~ log(tissue_weight), data = rna) 

# Calculates standardized residuals
rna$resid<- resid(l.mod)/sd(resid(l.mod))
# store predicted values for diagnostic plotting
rna$pred<- predict(l.mod)

# plot predicted vs. standardized residuals
rna %>%
  mutate(outlier = if_else(resid < -2| resid > 2, "out", "in")) %>%
  # filter(resid > 3) %>%
  ggplot(aes(log(tissue_weight), log(rna), color = time)) + geom_point() 



### Save data 

rna.save <- rna %>%
#  mutate(outlier = if_else(resid < -2.5| resid > 2.5, "out", "in")) %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  print()
  

saveRDS(rna.save, "./data/derivedData/tot-rna/tot-rna.RDS")



