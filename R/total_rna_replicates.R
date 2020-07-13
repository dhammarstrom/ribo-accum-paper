##### Total RNA from replicates #####



library(tidyverse); library(readxl)



files <- list.files("./data/tot_rna/total-rna-replicates/")



results <- list()

for(i in 1:length(files)){
  
  
  results[[i]] <-  read_excel(paste0("./data/tot_rna/total-rna-replicates/", files[i]), skip = 0, na = "NaN")
  
  
}


tot_rna <- bind_rows(results) %>%
  filter(!is.na(`Dilution factor`)) %>%
  dplyr::select(Sample, concentration = `...7`) %>%
  print()



# Extraction numbers round 2

# Filter P18
p18_samples <- tot_rna %>%
  filter(Sample %in% paste0("P18-", rep(1:6))) %>%
  print()


# Elution volume for P19, 21, 22, 23: 25 ul. 
tot_rna %>%
  filter(!(Sample %in% paste0("P18-", rep(1:6)))) %>%
  separate(Sample, into = c("participant", "series", "sample"), convert = TRUE) %>%
  inner_join(read.csv("./data/written/extraction_numbers_round2.csv", sep = ";") %>%
               separate(extraction_nr, into = c("series", "sample"), convert = TRUE) %>%
               mutate(series = if_else(series == "I", 1, 2))) %>%
  mutate(rna = concentration * 25) %>%
  dplyr::select(participant, series, sample, leg, time, tissue_weight, rna) %>%
  print()




  
  
  ggplot(aes(tissue_weight, rna)) + geom_point() +
  theme(legend.position = "none")
  
  
  print()

















