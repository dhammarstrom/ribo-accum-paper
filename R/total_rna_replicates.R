##### Total RNA from replicates #####




library(tidyverse); library(readxl)



files <- list.files("./data/total-rna-replicates/")



results <- list()

for(i in 1:length(files)){
  
  
  results[[i]] <-  read_excel(paste0("./data/total-rna-replicates/", files[i]), skip = 0, na = "NaN")
  
  
}


tot_rna <- bind_rows(results) %>%
  filter(!is.na(`Dilution factor`)) %>%
  dplyr::select(Sample, concentration = `...7`) %>%
  print()



