
#### Total RNA round 2 ####


library(tidyverse); library(readxl)



files <- list.files("./data/tot_rna/")

results <- list()

for(i in 1:length(files)){
  
  
results[[i]] <-  read_excel(paste0("./data/tot_rna/", files[i]), skip = 4)
  
  
}


tot_rna <- bind_rows(results) 


tot_rna %>%
  dplyr::select(Sample, conc = AVG) %>%
  separate(Sample, into = c("Subject", "Rep", "Sample")) %>%
  mutate(Stock = round(400/conc,2), 
         H2O = round(11 - Stock, 2)) %>%
  write_csv(path = "./data/written/cDNA_synthesis.csv") %>%
  print()
  
 

## Check correlation with tissue weight

tw <- read.csv("./data/written/extraction_numbers_round2.csv", sep = ";") %>%
  dplyr::select(participant, extraction_nr, tissue_weight) %>%
  separate(extraction_nr, into = c("rep", "sample")) %>%
  mutate(rep = as.character(if_else(rep == "I", 1, 2))) %>%
  
  print()
 
tw2 <- tot_rna %>%
  dplyr::select(Sample, conc = AVG) %>%
  separate(Sample, into = c("participant", "rep", "sample"), sep = "_") %>%
  inner_join(tw) %>%
  print()
  
  tw2 %>%
    
    mutate(check = if_else(tissue_weight > 7.5 & conc < 150, "check", "ok")) %>%
    
    ggplot(aes(tissue_weight, conc, color = participant, shape = check)) + geom_point()

  
  
  tw2 %>%
    
    mutate(check = if_else(tissue_weight > 7.5 & conc < 150, "check", "ok")) %>%
    
  filter(check == "check") %>%
    print()
  
  
  
  
  
  
  
    
  





